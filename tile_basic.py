import tifffile
import numpy as np
import os
from tqdm import tqdm

file_path = 'Process_374_raw.ome.btf'
tile_size = 512  # Very small tiles to minimize memory usage
green_channel = 1  # Green in RGB (0=Red, 1=Green, 2=Blue)
output_dir = 'tiles'

os.makedirs(output_dir, exist_ok=True)

print("Basic memory-safe tiling approach...")

# Get image info first
with tifffile.TiffFile(file_path) as tif:
    # Get basic info without loading data
    page = tif.pages[0]
    print(f"Page shape: {page.shape}")
    print(f"Page dtype: {page.dtype}")
    print(f"Page axes: {page.axes}")
    
    # Try to read just the first page to understand structure
    try:
        first_page = page.asarray()
        print(f"First page shape: {first_page.shape}")
        
        # If it's a multi-page TIFF, we need to handle it differently
        if len(tif.pages) > 1:
            print(f"Multi-page TIFF with {len(tif.pages)} pages")
            # For now, let's just process the first page
            img_data = first_page
        else:
            img_data = first_page
            
    except Exception as e:
        print(f"Error reading first page: {e}")
        print("Trying alternative approach...")
        
        # Try reading with different parameters
        try:
            img_data = tif.asarray()
            print(f"Image data shape: {img_data.shape}")
        except Exception as e2:
            print(f"Alternative approach failed: {e2}")
            exit(1)

# Determine dimensions
if len(img_data.shape) == 3:
    if img_data.shape[2] == 3:  # RGB
        height, width, channels = img_data.shape
        print(f"RGB image: {height}x{width}x{channels}")
    else:
        height, width = img_data.shape[:2]
        channels = 1
        print(f"Multi-channel image: {height}x{width}x{channels}")
        green_channel = 0
else:
    height, width = img_data.shape
    channels = 1
    green_channel = 0
    print(f"Grayscale image: {height}x{width}")

print(f'Processing {height}x{width} image with {channels} channels.')

# Extract green channel if needed
if channels > 1:
    green_image = img_data[:, :, green_channel]
else:
    green_image = img_data

print(f"Green channel shape: {green_image.shape}")

# Calculate number of tiles
total_tiles = ((height + tile_size - 1) // tile_size) * ((width + tile_size - 1) // tile_size)
print(f'Will create {total_tiles} tiles')

tile_idx = 0
with tqdm(total=total_tiles, desc="Creating tiles") as pbar:
    for y in range(0, height, tile_size):
        for x in range(0, width, tile_size):
            tile_h = min(tile_size, height - y)
            tile_w = min(tile_size, width - x)

            try:
                # Extract tile from green channel
                tile_data = green_image[y:y+tile_h, x:x+tile_w]

                # Save as TIFF
                output_file = os.path.join(output_dir, f'Process_374_green_tile_{tile_idx:04d}_{y}_{x}.tif')
                tifffile.imwrite(output_file, tile_data, photometric='minisblack', compression='zlib')
                
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

