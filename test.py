import tifffile

file_path = 'Process_374_raw.ome.btf'
with tifffile.TiffFile(file_path) as tif:
    print('Dimensions:', tif.series[0].shape)  # e.g., (channels, height, width) or similar
    print('Axes:', tif.series[0].axes)  # e.g., 'CYX' for channel, y, x
    print('Bit depth:', tif.series[0].dtype)
    ome_xml = tif.ome_metadata
    if ome_xml:
        print('OME XML snippet:', ome_xml[:1000])  # First 1000 chars; full XML may be long
    else:
        print('No OME metadata found.')
