import cv2
import numpy as np
import os
from tqdm import tqdm

file_path = 'Process_374_raw.ome.btf'
tile_size = 2048  # Reasonable tile size for memory efficiency
green_channel = 1  # Green in RGB (0=Red, 1=Green, 2=Blue)
output_dir = 'tiles'

os.makedirs(output_dir, exist_ok=True)

print("Processing image with OpenCV (memory-efficient tiling)...")

# Read image with OpenCV
print("Loading image with OpenCV...")
img = cv2.imread(file_path, cv2.IMREAD_UNCHANGED)

if img is None:
    print("Error: Could not load image with OpenCV")
    exit(1)

print(f"Image loaded: {img.shape}")
print(f"Image dtype: {img.dtype}")

# Check if image has multiple channels
if len(img.shape) == 3:
    height, width, channels = img.shape
    print(f"Processing {height}x{width} image with {channels} channels")
    
    # Extract green channel (OpenCV uses BGR order, so green is channel 1)
    green_image = img[:, :, green_channel]
else:
    height, width = img.shape
    channels = 1
    print(f"Processing {height}x{width} grayscale image")
    green_image = img

print(f"Green channel extracted: {green_image.shape}")

# Calculate total number of tiles for progress bar
total_tiles = ((height + tile_size - 1) // tile_size) * ((width + tile_size - 1) // tile_size)
print(f"Will create {total_tiles} tiles")

tile_idx = 0
with tqdm(total=total_tiles, desc="Creating tiles") as pbar:
    for y in range(0, height, tile_size):
        for x in range(0, width, tile_size):
            tile_h = min(tile_size, height - y)
            tile_w = min(tile_size, width - x)

            try:
                # Extract tile from green channel
                tile_data = green_image[y:y+tile_h, x:x+tile_w]

                # Save as TIFF using OpenCV
                output_file = os.path.join(output_dir, f'Process_374_green_tile_{tile_idx:04d}_{y}_{x}.tif')
                cv2.imwrite(output_file, tile_data)
                
                pbar.set_postfix({
                    'Current': f'{tile_h}x{tile_w}',
                    'Position': f'({y},{x})'
                })
                tile_idx += 1
                
            except Exception as e:
                print(f"Error processing tile at ({y},{x}): {e}")
                continue
            
            pbar.update(1)

print(f"Tiling complete! Created {tile_idx} tiles in '{output_dir}' directory.")
