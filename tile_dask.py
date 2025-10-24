import dask.array as da
import tifffile
import numpy as np
import os
from tqdm import tqdm

file_path = 'Process_374_raw.ome.btf'
tile_size = 2048  # Reasonable tile size
green_channel = 1  # Green in RGB (0=Red, 1=Green, 2=Blue)
output_dir = 'tiles'

os.makedirs(output_dir, exist_ok=True)

print("Loading image with Dask for memory-efficient processing...")

# Load image as Dask array for memory-efficient processing
with tifffile.TiffFile(file_path) as tif:
    img = tif.series[0]
    if 'S' in img.axes:  # RGB (YXS axes)
        num_channels = img.shape[img.axes.index('S')]
        height, width = img.shape[img.axes.index('Y')], img.shape[img.axes.index('X')]
    else:  # Grayscale
        num_channels = 1
        height, width = img.shape[-2], img.shape[-1]
        green_channel = 0

print(f'Processing {height}x{width} image with {num_channels} channels.')

# Create Dask array from the image
# This creates a lazy array that doesn't load data until needed
dask_img = da.from_array(img.asarray(), chunks=(tile_size, tile_size, num_channels))

print(f"Dask array created with chunks: {dask_img.chunks}")

# Extract green channel
if num_channels > 1:
    green_dask = dask_img[:, :, green_channel]
else:
    green_dask = dask_img

print(f"Green channel extracted: {green_dask.shape}")

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
                # Extract tile using Dask (memory efficient)
                tile_data = green_dask[y:y+tile_h, x:x+tile_w].compute()
                
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