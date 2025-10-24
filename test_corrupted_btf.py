#!/usr/bin/env python3
"""
Test script for corrupted BTF files
==================================

This script tests the corrupted OME workaround for BTF files with severely
corrupted OME metadata.

Usage:
    python test_corrupted_btf.py
"""

import sys
import os
from pathlib import Path

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))

from hpc_btf_processor import BTFProcessor

def test_corrupted_btf():
    """Test processing with corrupted OME workaround."""
    
    # Configuration for corrupted OME workaround
    config = {
        'input_dir': './input',
        'output_dir': './output',
        'tile_size': 1024,  # Smaller tiles for testing
        'green_channel': 1,
        'log_level': 'INFO',
        'log_file': 'corrupted_btf_test.log',
        'corrupted_ome_workaround': True  # Enable workaround
    }
    
    print("Testing corrupted OME workaround...")
    print(f"Configuration: {config}")
    
    try:
        processor = BTFProcessor(config)
        
        # Test with a single file
        btf_files = processor.get_btf_files()
        if not btf_files:
            print("No BTF files found in input directory")
            return
            
        print(f"Found {len(btf_files)} BTF files")
        
        # Process first file only for testing
        test_file = btf_files[0]
        print(f"Testing with file: {test_file}")
        
        # Get image info
        height, width, num_channels, axes = processor.get_image_info(test_file)
        print(f"Image info: {height}x{width}, {num_channels} channels, axes: {axes}")
        
        # Test extracting a small tile
        print("Testing tile extraction...")
        tile_data = processor.extract_tile(test_file, 0, 0, 100, 100, 1, num_channels, axes)
        print(f"Tile extracted successfully: {tile_data.shape}, dtype: {tile_data.dtype}")
        
        print("✅ Corrupted OME workaround test passed!")
        
    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    test_corrupted_btf()
