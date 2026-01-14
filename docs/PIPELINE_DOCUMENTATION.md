# CellProfiler Pipeline Documentation

## Overview

This project contains two CellProfiler pipelines for analyzing salamander blood cell images to measure:
- **Nucleus size** (area in pixels)
- **Cell size** (area in pixels) - for brightfield images only
- **IOD (Integrated Optical Density)** - correlates with DNA content/genome size

Based on methodology from: [Hardie, Gregory & Hebert (2002) "From Pixels to Picograms: A Beginners' Guide to Genome Quantification by Feulgen Image Analysis Densitometry"](https://journals.sagepub.com/doi/10.1177/002215540205000601)

---

## Important: True IOD Calculation

CellProfiler's `IntegratedIntensity` is **NOT** the same as true IOD from the Hardie et al. paper.

### The Difference

| Measurement | Formula | Use |
|-------------|---------|-----|
| CellProfiler IntegratedIntensity | Σ pixel_values | Raw sum of intensities |
| **True IOD** | Σ log₁₀(I₀/I) | Proper optical density with background correction |

### True IOD Formula

For each pixel in the nucleus:
```
OD_pixel = log₁₀(I_background / I_pixel)
```

Then sum all pixels:
```
IOD = Σ OD_pixel (for all pixels in nucleus)
```

Where:
- `I_background` (I₀) = intensity of clear/unstained area (incident light)
- `I_pixel` (I) = intensity of each pixel in the nucleus

### Post-Processing Script

Use `calculate_true_iod.py` to compute proper IOD values:

```bash
python calculate_true_iod.py <image> <cellprofiler_nuclei.csv> <output.csv>
```

This script:
1. Estimates background from image (95th percentile)
2. Calculates true OD = log₁₀(I₀/I) for each pixel
3. Sums to get IOD per nucleus
4. Applies QC: removes top/bottom 5% of values (per Hardie et al.)
5. Reports mean IOD with coefficient of variation (CV)

### Genome Size Calculation

To convert IOD to genome size (picograms), you need an **internal standard** with known genome size:

```
Genome_size_unknown = (IOD_unknown / IOD_standard) × Genome_size_standard
```

Common standards:
- Chicken erythrocytes: 1C = 1.25 pg
- Rainbow trout erythrocytes: 1C = 2.60 pg

### Quality Control (per Hardie et al.)

1. **Minimum nuclei**: Measure ≥20 nuclei per sample
2. **Trim outliers**: Remove top and bottom 5% of IOD values
3. **Eliminate poor slides**: Discard slides with faint nuclei or pink backgrounds
4. **CV threshold**: CV should typically be <20% for good preparations

---

## Pipeline 1: Oil Immersion (`pipeline_oil_immersion.cppipe`)

**Use for:** Oil immersion microscopy images where only the nucleus is visible (no cytoplasm)

### Processing Steps

#### Step 1: Load Images (`Images` module)
- Filters input to only process image files
- Excludes hidden directories

#### Step 2: Extract Metadata (`Metadata` module)
- Currently disabled (no metadata extraction)
- Can be configured to extract sample info from filenames

#### Step 3: Name and Type Images (`NamesAndTypes` module)
- Assigns loaded images the name `OrigGreen`
- Treats as grayscale images
- Sets intensity range from image metadata (0-255 for uint8)

#### Step 4: Group Images (`Groups` module)
- Currently disabled (processes all images individually)
- Can group by metadata for batch analysis

#### Step 5: Invert Image (`ImageMath` module)
- **Operation:** Invert
- **Input:** `OrigGreen`
- **Output:** `GreenChannel`
- **Purpose:** Converts dark-on-bright (Feulgen staining) to bright-on-dark for proper segmentation
- Dark nuclei become bright objects that CellProfiler can detect

#### Step 6: Identify Nuclei (`IdentifyPrimaryObjects` module)
- **Input:** `GreenChannel` (inverted image)
- **Output:** `Nuclei` objects
- **Parameters:**
  - Typical diameter: 30-200 pixels (sized for large salamander nuclei)
  - Discard objects outside diameter range: Yes
  - Discard border objects: Yes (removes partial cells at image edges)
  - Declumping method: Intensity-based (separates touching nuclei)
  - Thresholding: Global Otsu, three-class
  - Fill holes: Yes (after thresholding and declumping)

#### Step 7: Measure Size and Shape (`MeasureObjectSizeShape` module)
- **Objects measured:** `Nuclei`
- **Key outputs:**
  - `AreaShape_Area` - nucleus area in pixels
  - `AreaShape_Center_X/Y` - centroid coordinates
  - `AreaShape_Perimeter` - boundary length
  - `AreaShape_Eccentricity` - shape elongation
  - `AreaShape_FormFactor` - circularity measure

#### Step 8: Measure Intensity (`MeasureObjectIntensity` module)
- **Image:** `OrigGreen` (original, non-inverted image)
- **Objects:** `Nuclei`
- **Key outputs:**
  - `Intensity_IntegratedIntensity_OrigGreen` - **IOD (Integrated Optical Density)**
  - `Intensity_MeanIntensity_OrigGreen` - average pixel intensity
  - `Intensity_StdIntensity_OrigGreen` - intensity variation

> **Note:** IOD is measured on the ORIGINAL image (not inverted) to get correct optical density values for DNA quantification.

#### Step 9: Export Results (`ExportToSpreadsheet` module)
- **Output format:** CSV files
- **Files generated:**
  - `Oil_Nuclei.csv` - per-nucleus measurements
  - `Oil_Image.csv` - per-image summary statistics
- **Includes:** Mean, median, and standard deviation per image

---

## Pipeline 2: Brightfield (`pipeline_brightfield.cppipe`)

**Use for:** Brightfield microscopy images where both nucleus and cytoplasm are visible

### Processing Steps

Steps 1-6 are identical to the oil immersion pipeline.

#### Step 7: Identify Cells (`IdentifySecondaryObjects` module)
- **Input objects:** `Nuclei` (as seeds)
- **Input image:** `GreenChannel` (inverted)
- **Output:** `Cells` objects
- **Method:** Propagation (expands from nuclei to find cell boundaries)
- **Parameters:**
  - Maximum expansion: 100 pixels
  - Regularization factor: 0.01 (allows flexible boundaries)
  - Thresholding: Adaptive Otsu, three-class
  - Adaptive window: 100 pixels
  - Discard border cells: Yes
- **Also outputs:** `FilteredNuclei` (nuclei with valid cells)

#### Step 8: Measure Size and Shape (`MeasureObjectSizeShape` module)
- **Objects measured:** `Nuclei` AND `Cells`
- Enables calculation of nucleus-to-cell ratio

#### Step 9: Measure Intensity (`MeasureObjectIntensity` module)
- **Image:** `OrigGreen` (original image)
- **Objects:** `Nuclei`
- Measures IOD on original image for accurate DNA quantification

#### Step 10: Export Results (`ExportToSpreadsheet` module)
- **Files generated:**
  - `Measurements_Nuclei.csv` - per-nucleus measurements including IOD
  - `Measurements_Cells.csv` - per-cell measurements
  - `Measurements_Image.csv` - per-image summary statistics

---

## Key Measurements Explained

### Nucleus/Cell Size
- **AreaShape_Area**: Total pixels within the object boundary
- Convert to physical units using: `Area_µm² = Area_pixels × (pixel_size_µm)²`

### IOD (Integrated Optical Density)
- **Intensity_IntegratedIntensity_OrigGreen**: Sum of all pixel intensities within the nucleus
- Proportional to total stain amount (DNA content for Feulgen staining)
- Used for genome size estimation when compared to a known standard

### Nucleus/Cell Ratio
- Calculate from: `N/C_ratio = Nucleus_Area / Cell_Area`
- Only applicable for brightfield images

---

## Usage

### Command Line (Headless)
```bash
# Activate environment
source ~/miniforge3/bin/activate
conda activate cellprofiler

# Oil immersion images
cellprofiler -p pipeline_oil_immersion.cppipe -c -r -o output -i /path/to/images

# Brightfield images
cellprofiler -p pipeline_brightfield.cppipe -c -r -o output -i /path/to/images

# Single file or file list
cellprofiler -p pipeline_oil_immersion.cppipe -c -r -o output --file-list filelist.txt
```

### Parameters to Tune

If detection is not accurate, adjust these parameters:

1. **Nucleus size** (`IdentifyPrimaryObjects`):
   - `Typical diameter (Min,Max)`: Currently 30-200 pixels
   - Increase for larger cells, decrease for smaller

2. **Threshold sensitivity** (`IdentifyPrimaryObjects`):
   - `Threshold correction factor`: Default 1.0
   - Increase (>1) to detect fainter nuclei
   - Decrease (<1) to be more selective

3. **Cell expansion** (`IdentifySecondaryObjects` - brightfield only):
   - `Number of pixels to expand`: Currently 100
   - Increase if cells are larger than detected

---

## Output Files

### Per-Object CSVs
Each row = one detected nucleus/cell with columns:
- `ImageNumber` - source image index
- `ObjectNumber` - object ID within image
- `AreaShape_Area` - size in pixels
- `Intensity_IntegratedIntensity_OrigGreen` - IOD value
- Plus 20+ additional shape and intensity measurements

### Per-Image CSVs
Summary statistics across all objects in each image:
- Mean, median, standard deviation of all measurements
- Object counts
- Image metadata

---

## Quality Control

### Visual Verification
Generate overlay images to check segmentation:
```python
# See output/oil_segmentation_zoom.png or output/segmentation_check.png
```

### Common Issues

| Issue | Solution |
|-------|----------|
| Too few objects detected | Lower threshold correction factor or increase diameter range |
| Too many false detections | Increase minimum diameter or raise threshold |
| Merged/clumped nuclei | Adjust declumping settings or smoothing filter |
| Cells not expanding | Increase expansion pixels or lower secondary threshold |

---

## File Organization

```
cellprofiler_test/
├── pipeline_oil_immersion.cppipe    # For oil immersion images
├── pipeline_brightfield.cppipe       # For brightfield images
├── data/                             # Converted green channel images
│   ├── *_green.ome.tiff             # Input images
│   └── test_crop.tiff               # Test image (2000x2000 crop)
├── output/                           # CellProfiler results
│   ├── Oil_Nuclei.csv               # Oil immersion measurements
│   ├── Measurements_Nuclei.csv      # Brightfield nucleus measurements
│   ├── Measurements_Cells.csv       # Brightfield cell measurements
│   └── *_segmentation*.png          # Visual overlays
└── PIPELINE_DOCUMENTATION.md         # This file
```
