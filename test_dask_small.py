import dask.array as da
import tifffile
import numpy as np
import os

file_path = 'Process_374_raw.ome.btf'
test_size = 1000
green_channel = 1
output_dir = 'test_tiles'

os.makedirs(output_dir, exist_ok=True)

print("Testing Dask approach with small region...")

try:
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
        
        # Create Dask array with small chunks for testing
        dask_img = da.from_array(img.asarray(), chunks=(test_size, test_size, num_channels))
        print(f"Dask array created: {dask_img.shape}")
        
        # Extract green channel
        if num_channels > 1:
            green_dask = dask_img[:, :, green_channel]
        else:
            green_dask = dask_img
        
        print(f"Green channel shape: {green_dask.shape}")
        
        # Test reading a small region
        test_region = green_dask[0:test_size, 0:test_size].compute()
        print(f'Test region shape: {test_region.shape}')
        print(f'Test region dtype: {test_region.dtype}')
        print(f'Test region range: {test_region.min()} to {test_region.max()}')
        
        # Save test tile
        output_file = os.path.join(output_dir, 'test_green_tile_dask.tif')
        tifffile.imwrite(output_file, test_region, photometric='minisblack', compression='zlib')
        print(f'Test tile saved: {output_file}')

except Exception as e:
    print(f"Error: {e}")
    print("Dask approach failed, trying alternative method...")

print("Test complete!")


