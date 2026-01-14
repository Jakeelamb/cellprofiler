#!/usr/bin/env python3
"""
Calculate true IOD (Integrated Optical Density) following Hardie et al. (2002)
"From Pixels to Picograms" methodology.

IOD = Σ log₁₀(I_background / I_pixel) for all pixels in nucleus

This script post-processes CellProfiler output to calculate proper IOD values
with background correction.

Usage:
    python calculate_true_iod.py <image_path> <nuclei_csv> <output_csv>
"""

import argparse
import csv
import math
import numpy as np
import tifffile
from pathlib import Path


def estimate_background(image: np.ndarray, percentile: float = 95) -> float:
    """
    Estimate background intensity (I₀) from the image.

    Uses the high percentile of pixel values as background estimate,
    assuming most of the bright areas are unstained background.

    Args:
        image: Grayscale image array
        percentile: Percentile to use for background (default 95)

    Returns:
        Background intensity value
    """
    return np.percentile(image, percentile)


def calculate_pixel_od(pixel_value: float, background: float, epsilon: float = 1.0) -> float:
    """
    Calculate optical density for a single pixel.

    OD = log₁₀(I_background / I_pixel)

    Args:
        pixel_value: Pixel intensity (I)
        background: Background intensity (I₀)
        epsilon: Small value to prevent log(0)

    Returns:
        Optical density value
    """
    # Prevent division by zero and log of zero
    pixel_value = max(pixel_value, epsilon)
    if pixel_value >= background:
        return 0.0  # No absorption (bright as background)

    return math.log10(background / pixel_value)


def calculate_nucleus_iod(image: np.ndarray, center_x: float, center_y: float,
                          area: float, background: float) -> dict:
    """
    Calculate true IOD for a nucleus using circular approximation.

    Args:
        image: Grayscale image
        center_x, center_y: Nucleus centroid
        area: Nucleus area in pixels
        background: Background intensity

    Returns:
        Dictionary with IOD and related measurements
    """
    radius = int(math.sqrt(area / math.pi))
    cx, cy = int(center_x), int(center_y)

    # Get pixels within nucleus (circular mask)
    y_indices, x_indices = np.ogrid[:image.shape[0], :image.shape[1]]
    mask = (x_indices - cx)**2 + (y_indices - cy)**2 <= radius**2

    # Ensure within image bounds
    nucleus_pixels = image[mask]

    if len(nucleus_pixels) == 0:
        return {'IOD': 0, 'MeanOD': 0, 'PixelCount': 0}

    # Calculate OD for each pixel and sum (IOD)
    od_values = [calculate_pixel_od(float(p), background) for p in nucleus_pixels]
    iod = sum(od_values)
    mean_od = iod / len(od_values) if od_values else 0

    return {
        'IOD': iod,
        'MeanOD': mean_od,
        'PixelCount': len(nucleus_pixels),
        'MaxOD': max(od_values) if od_values else 0,
        'MinOD': min(od_values) if od_values else 0
    }


def process_image(image_path: str, nuclei_csv: str, output_csv: str,
                  background_percentile: float = 95,
                  trim_percent: float = 5.0):
    """
    Process an image and calculate true IOD for all nuclei.

    Args:
        image_path: Path to the image file
        nuclei_csv: Path to CellProfiler nuclei measurements CSV
        output_csv: Path for output CSV with true IOD
        background_percentile: Percentile for background estimation
        trim_percent: Percent to trim from top/bottom of IOD values (QC)
    """
    # Load image
    print(f"Loading image: {image_path}")
    image = tifffile.imread(image_path)

    if image.ndim == 3:
        # If RGB, use green channel
        image = image[:, :, 1]

    # Estimate background
    background = estimate_background(image, background_percentile)
    print(f"Background intensity (I₀): {background:.1f}")

    # Load nuclei measurements
    with open(nuclei_csv, 'r') as f:
        reader = csv.DictReader(f)
        nuclei = list(reader)

    print(f"Processing {len(nuclei)} nuclei...")

    # Calculate true IOD for each nucleus
    results = []
    for n in nuclei:
        cx = float(n['AreaShape_Center_X'])
        cy = float(n['AreaShape_Center_Y'])
        area = float(n['AreaShape_Area'])

        iod_data = calculate_nucleus_iod(image, cx, cy, area, background)

        results.append({
            'ObjectNumber': n['ObjectNumber'],
            'Center_X': cx,
            'Center_Y': cy,
            'Area_pixels': area,
            'True_IOD': iod_data['IOD'],
            'Mean_OD': iod_data['MeanOD'],
            'Max_OD': iod_data['MaxOD'],
            'Background_I0': background,
            'CellProfiler_IntegratedIntensity': n.get('Intensity_IntegratedIntensity_OrigGreen', 'N/A')
        })

    # Sort by IOD for QC filtering
    results.sort(key=lambda x: x['True_IOD'])
    n_total = len(results)

    # Apply QC: remove top and bottom trim_percent
    if trim_percent > 0 and n_total > 10:
        trim_n = int(n_total * trim_percent / 100)
        results_qc = results[trim_n:-trim_n] if trim_n > 0 else results
        print(f"QC: Removed {trim_n} from top and bottom ({trim_percent}%)")
        print(f"Retained {len(results_qc)} of {n_total} nuclei")
    else:
        results_qc = results

    # Calculate summary statistics
    iod_values = [r['True_IOD'] for r in results_qc]
    if iod_values:
        mean_iod = sum(iod_values) / len(iod_values)
        variance = sum((x - mean_iod)**2 for x in iod_values) / len(iod_values)
        std_iod = math.sqrt(variance)
        cv = (std_iod / mean_iod * 100) if mean_iod > 0 else 0

        print(f"\n=== IOD Summary (after QC) ===")
        print(f"N nuclei: {len(iod_values)}")
        print(f"Mean IOD: {mean_iod:.2f}")
        print(f"Std IOD: {std_iod:.2f}")
        print(f"CV: {cv:.1f}%")
        print(f"Range: {min(iod_values):.2f} - {max(iod_values):.2f}")

    # Write output CSV
    with open(output_csv, 'w', newline='') as f:
        fieldnames = list(results[0].keys()) + ['QC_Passed']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        qc_set = set(r['ObjectNumber'] for r in results_qc)
        for r in results:
            r['QC_Passed'] = 'Yes' if r['ObjectNumber'] in qc_set else 'No'
            writer.writerow(r)

    print(f"\nResults saved to: {output_csv}")

    return results_qc


def main():
    parser = argparse.ArgumentParser(
        description="Calculate true IOD following Hardie et al. (2002) methodology"
    )
    parser.add_argument("image_path", help="Path to image file")
    parser.add_argument("nuclei_csv", help="Path to CellProfiler nuclei CSV")
    parser.add_argument("output_csv", help="Output CSV path")
    parser.add_argument("--background-percentile", type=float, default=95,
                        help="Percentile for background estimation (default: 95)")
    parser.add_argument("--trim-percent", type=float, default=5.0,
                        help="Percent to trim from top/bottom for QC (default: 5)")

    args = parser.parse_args()

    process_image(
        args.image_path,
        args.nuclei_csv,
        args.output_csv,
        args.background_percentile,
        args.trim_percent
    )


if __name__ == "__main__":
    main()
