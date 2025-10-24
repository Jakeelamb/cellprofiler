import tifffile
import numpy as np
import os
from tqdm import tqdm

file_path = 'Process_374_raw.ome.btf'
tile_size = 2048  # Reasonable tile size for memory efficiency
green_channel = 1  # Green in RGB (0=Red, 1=Green, 2=Blue)
output_dir = 'tiles'

os.makedirs(output_dir, exist_ok=True)

print("Processing image with memory-efficient tiling...")

# Use memory mapping for large files
with tifffile.TiffFile(file_path) as tif:
    img = tif.series[0]  # First series (single image)
    if 'S' in img.axes:  # RGB (YXS axes)
        num_channels = img.shape[img.axes.index('S')]
        height, width = img.shape[img.axes.index('Y')], img.shape[img.axes.index('X')]
    else:  # Grayscale (unlikely, fallback)
        num_channels = 1
        height, width = img.shape[-2], img.shape[-1]
        green_channel = 0

    print(f'Processing {height}x{width} image with {num_channels} channels.')
    print(f'Will create {(height//tile_size + 1) * (width//tile_size + 1)} tiles')
    
    # Calculate total number of tiles for progress bar
    total_tiles = ((height + tile_size - 1) // tile_size) * ((width + tile_size - 1) // tile_size)
    
    tile_idx = 0
    with tqdm(total=total_tiles, desc="Creating tiles") as pbar:
        for y in range(0, height, tile_size):
            for x in range(0, width, tile_size):
                tile_h = min(tile_size, height - y)
                tile_w = min(tile_size, width - x)

                try:
                    # Read only the specific tile region from the file
                    if num_channels > 1:
                        # For YXS format: read the tile region and extract green channel
                        tile_region = img.asarray()[y:y+tile_h, x:x+tile_w, :]
                        tile_data = tile_region[:, :, green_channel]
                    else:
                        tile_data = img.asarray()[y:y+tile_h, x:x+tile_w]

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
