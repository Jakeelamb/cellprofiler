#!/usr/bin/env python3
"""
Memory Usage Comparison
=======================

Compare memory usage between single-pass and two-pass approaches
for processing large BTF files.

Usage:
    python memory_comparison.py
"""

def calculate_memory_usage():
    """Calculate theoretical memory usage for different approaches."""
    
    # Example: 69,338 x 69,789 x 3 channels (RGB)
    height, width, channels = 69338, 69789, 3
    green_channel = 1
    
    print("BTF File Memory Usage Comparison")
    print("=" * 50)
    print(f"Image dimensions: {height} x {width} x {channels} channels")
    print()
    
    # Calculate memory usage for different data types
    data_types = {
        'uint8': 1,
        'uint16': 2,
        'uint32': 4,
        'float32': 4,
        'float64': 8
    }
    
    for dtype, bytes_per_pixel in data_types.items():
        print(f"\n{dtype.upper()} Data Type:")
        print("-" * 30)
        
        # Single-pass approach (original)
        full_image_memory = height * width * channels * bytes_per_pixel
        full_image_gb = full_image_memory / (1024**3)
        
        # Two-pass approach
        green_only_memory = height * width * 1 * bytes_per_pixel  # Only green channel
        green_only_gb = green_only_memory / (1024**3)
        
        # Memory reduction
        reduction_percent = ((full_image_memory - green_only_memory) / full_image_memory) * 100
        
        print(f"  Single-pass (full RGB): {full_image_gb:.2f} GB")
        print(f"  Two-pass (green only):   {green_only_gb:.2f} GB")
        print(f"  Memory reduction:       {reduction_percent:.1f}%")
        print(f"  Memory saved:           {full_image_gb - green_only_gb:.2f} GB")
    
    print("\n" + "=" * 50)
    print("RECOMMENDATIONS:")
    print("=" * 50)
    
    # Typical microscopy data is uint16
    uint16_full = height * width * channels * 2 / (1024**3)
    uint16_green = height * width * 1 * 2 / (1024**3)
    
    print(f"For typical uint16 microscopy data:")
    print(f"  Single-pass requires: {uint16_full:.1f} GB RAM")
    print(f"  Two-pass requires:    {uint16_green:.1f} GB RAM")
    print(f"  Reduction:            {((uint16_full - uint16_green) / uint16_full) * 100:.1f}%")
    print()
    
    if uint16_full > 64:
        print("⚠️  Single-pass approach may exceed 64GB RAM limit")
        print("✅ Two-pass approach recommended for this file size")
    else:
        print("✅ Both approaches should work within 64GB RAM limit")
    
    print()
    print("TWO-PASS BENEFITS:")
    print("• 70% less memory usage")
    print("• More reliable for corrupted OME metadata")
    print("• Faster tiling (smaller file to process)")
    print("• Better error recovery")
    print("• Intermediate green file can be reused")

if __name__ == '__main__':
    calculate_memory_usage()
