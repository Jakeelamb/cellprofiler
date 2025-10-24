import tifffile
import numpy as np
import os

file_path = 'Process_374_raw.ome.btf'
test_size = 1000  # Small test region
green_channel = 1
output_dir = 'test_tiles'

os.makedirs(output_dir, exist_ok=True)

print("Testing small region tiling...")

with tifffile.TiffFile(file_path) as tif:
    img = tif.series[0]
    if 'S' in img.axes:
        num_channels = img.shape[img.axes.index('S')]
        height, width = img.shape[img.axes.index('Y')], img.shape[img.axes.index('X')]
    else:
        num_channels = 1
        height, width = img.shape[-2], img.shape[-1]
        green_channel = 0

    print(f'Image dimensions: {height}x{width} with {num_channels} channels')
    
    # Test reading a small region
    test_region = img.asarray(key=(slice(0, test_size), slice(0, test_size), slice(None)))
    print(f'Test region shape: {test_region.shape}')
    
    if num_channels > 1:
        green_test = test_region[:, :, green_channel]
    else:
        green_test = test_region
    
    print(f'Green channel shape: {green_test.shape}')
    print(f'Green channel dtype: {green_test.dtype}')
    print(f'Green channel range: {green_test.min()} to {green_test.max()}')
    
    # Save test tile
    output_file = os.path.join(output_dir, 'test_green_tile.tif')
    tifffile.imwrite(output_file, green_test, photometric='minisblack', compression='zlib')
    print(f'Test tile saved: {output_file}')

print("Test complete!")




