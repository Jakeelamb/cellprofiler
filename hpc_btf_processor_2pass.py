#!/usr/bin/env python3
"""
Two-Pass BTF File Processor
===========================

A memory-efficient pipeline that first extracts the green channel to a new file,
then tiles the smaller green-only file. This reduces memory usage by ~70%.

Usage:
    python hpc_btf_processor_2pass.py --input-dir /path/to/btf/files --output-dir /path/to/output
    python hpc_btf_processor_2pass.py --config config.yaml
"""

import argparse
import logging
import os
import sys
import time
import psutil
import shutil
from pathlib import Path
from typing import List, Tuple, Optional, Dict
import yaml
import tifffile
import numpy as np
from tqdm import tqdm
import json
from datetime import datetime


class TwoPassBTFProcessor:
    """Two-pass processor for BTF files with corrupted OME metadata."""
    
    def __init__(self, config: dict):
        """Initialize processor with configuration."""
        self.config = config
        self.setup_logging()
        self.validate_config()
        
    def setup_logging(self):
        """Setup logging configuration."""
        log_level = getattr(logging, self.config.get('log_level', 'INFO').upper())
        log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        
        logging.basicConfig(
            level=log_level,
            format=log_format,
            handlers=[
                logging.FileHandler(self.config['log_file']),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def validate_config(self):
        """Validate configuration parameters."""
        required_keys = ['input_dir', 'output_dir', 'tile_size', 'green_channel']
        for key in required_keys:
            if key not in self.config:
                raise ValueError(f"Missing required configuration key: {key}")
                
        # Ensure directories exist
        Path(self.config['input_dir']).mkdir(parents=True, exist_ok=True)
        Path(self.config['output_dir']).mkdir(parents=True, exist_ok=True)
        
        self.logger.info(f"Configuration validated successfully")
        
    def get_file_size_mb(self, file_path: Path) -> float:
        """Get file size in MB."""
        return file_path.stat().st_size / (1024 * 1024)
    
    def get_memory_usage_mb(self) -> Dict[str, float]:
        """Get current memory usage in MB."""
        process = psutil.Process()
        memory_info = process.memory_info()
        return {
            'rss_mb': memory_info.rss / (1024 * 1024),  # Resident Set Size
            'vms_mb': memory_info.vms / (1024 * 1024),  # Virtual Memory Size
            'percent': process.memory_percent()
        }
    
    def log_memory_usage(self, stage: str):
        """Log current memory usage."""
        memory = self.get_memory_usage_mb()
        self.logger.info(f"[{stage}] Memory usage: RSS={memory['rss_mb']:.1f}MB, "
                        f"VMS={memory['vms_mb']:.1f}MB, Percent={memory['percent']:.1f}%")
        
    def log_file_sizes(self, files: Dict[str, Path], stage: str):
        """Log file sizes for tracking."""
        self.logger.info(f"[{stage}] File sizes:")
        for name, path in files.items():
            if path.exists():
                size_mb = self.get_file_size_mb(path)
                self.logger.info(f"  {name}: {size_mb:.2f} MB ({path.name})")
            else:
                self.logger.info(f"  {name}: File not found ({path.name})")
        
    def get_btf_files(self) -> List[Path]:
        """Get list of BTF files to process."""
        input_dir = Path(self.config['input_dir'])
        btf_files = list(input_dir.glob('*.btf')) + list(input_dir.glob('*.ome.btf'))
        
        if not btf_files:
            raise FileNotFoundError(f"No BTF files found in {input_dir}")
            
        self.logger.info(f"Found {len(btf_files)} BTF files to process")
        return sorted(btf_files)
        
    def extract_green_channel(self, file_path: Path) -> Path:
        """Extract green channel to a new file (Pass 1)."""
        self.logger.info(f"Pass 1: Extracting green channel from {file_path.name}")
        
        # Log initial memory and file sizes
        self.log_memory_usage("PASS1_START")
        self.log_file_sizes({"Original BTF": file_path}, "PASS1_START")
        
        # Create output path for green channel file
        green_file_path = Path(self.config['output_dir']) / f"{file_path.stem}_green.tif"
        
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                
                # Read the entire image (memory intensive but necessary for corrupted OME)
                self.logger.warning("Reading entire image to extract green channel...")
                self.log_memory_usage("PASS1_READING_IMAGE")
                
                start_time = time.time()
                full_image = tifffile.imread(file_path, key=0)
                read_time = time.time() - start_time
                
                self.logger.info(f"Image read completed in {read_time:.2f} seconds")
                self.log_memory_usage("PASS1_IMAGE_LOADED")
                
                # Log image properties
                self.logger.info(f"Original image properties:")
                self.logger.info(f"  Shape: {full_image.shape}")
                self.logger.info(f"  Data type: {full_image.dtype}")
                self.logger.info(f"  Memory size: {full_image.nbytes / (1024**3):.2f} GB")
                
                # Extract green channel
                self.logger.info(f"Extracting green channel (index {self.config['green_channel']})...")
                if len(full_image.shape) == 3:
                    green_channel = full_image[:, :, self.config['green_channel']]
                    self.logger.info(f"Green channel extracted from RGB image")
                else:
                    # Single channel, use as is
                    green_channel = full_image
                    self.logger.info(f"Single channel image, using as green channel")
                
                # Log green channel properties
                self.logger.info(f"Green channel properties:")
                self.logger.info(f"  Shape: {green_channel.shape}")
                self.logger.info(f"  Data type: {green_channel.dtype}")
                self.logger.info(f"  Memory size: {green_channel.nbytes / (1024**3):.2f} GB")
                
                # Calculate size reduction
                original_size = full_image.nbytes
                green_size = green_channel.nbytes
                reduction_percent = ((original_size - green_size) / original_size) * 100
                self.logger.info(f"Size reduction: {reduction_percent:.1f}% "
                               f"({original_size / (1024**3):.2f} GB → {green_size / (1024**3):.2f} GB)")
                
                # Save green channel as new TIFF file (without OME metadata)
                self.logger.info(f"Saving green channel to {green_file_path}")
                self.log_memory_usage("PASS1_SAVING_GREEN")
                
                save_start = time.time()
                tifffile.imwrite(
                    green_file_path, 
                    green_channel, 
                    photometric='minisblack', 
                    compression='zlib',
                    metadata=None  # Remove all metadata to avoid OME corruption
                )
                save_time = time.time() - save_start
                
                self.logger.info(f"Green channel saved in {save_time:.2f} seconds")
                self.log_memory_usage("PASS1_GREEN_SAVED")
                
                # Log final file sizes
                self.log_file_sizes({
                    "Original BTF": file_path,
                    "Green Channel": green_file_path
                }, "PASS1_COMPLETE")
                
                return green_file_path
                
        except Exception as e:
            self.logger.error(f"Failed to extract green channel from {file_path}: {e}")
            raise
            
    def tile_green_file(self, green_file_path: Path, original_file_path: Path) -> dict:
        """Tile the green channel file (Pass 2)."""
        self.logger.info(f"Pass 2: Tiling green channel file {green_file_path.name}")
        
        # Log initial memory and file sizes
        self.log_memory_usage("PASS2_START")
        self.log_file_sizes({"Green Channel": green_file_path}, "PASS2_START")
        
        try:
            # Get image dimensions
            with tifffile.TiffFile(green_file_path) as tif:
                img = tif.series[0]
                height, width = img.shape
                
            self.logger.info(f"Green file dimensions: {height}x{width}")
            
            # Create output directory for tiles
            file_stem = original_file_path.stem.replace('.ome', '')
            tiles_dir = Path(self.config['output_dir']) / file_stem / 'green_tiles'
            tiles_dir.mkdir(parents=True, exist_ok=True)
            
            # Process tiles
            tile_size = self.config['tile_size']
            total_tiles = ((height + tile_size - 1) // tile_size) * ((width + tile_size - 1) // tile_size)
            
            self.logger.info(f"Creating {total_tiles} tiles of size {tile_size}x{tile_size}")
            self.logger.info(f"Tile grid: {((height + tile_size - 1) // tile_size)} rows x {((width + tile_size - 1) // tile_size)} cols")
            
            tile_info = []
            tile_idx = 0
            total_tile_size_mb = 0
            
            start_time = time.time()
            
            with tqdm(total=total_tiles, desc=f"Tiling {green_file_path.name}") as pbar:
                for y in range(0, height, tile_size):
                    for x in range(0, width, tile_size):
                        tile_h = min(tile_size, height - y)
                        tile_w = min(tile_size, width - x)
                        
                        try:
                            # Read tile from green file (much more memory efficient)
                            with tifffile.TiffFile(green_file_path) as tif:
                                img = tif.series[0]
                                try:
                                    # Try key-based slicing first
                                    tile_data = img.asarray(
                                        key=(slice(y, y+tile_h), slice(x, x+tile_w))
                                    )
                                except Exception as slicing_error:
                                    # Fallback: read entire image and slice (memory intensive but works)
                                    self.logger.warning(f"Key-based slicing failed, using direct array slicing: {slicing_error}")
                                    full_image = img.asarray()
                                    tile_data = full_image[y:y+tile_h, x:x+tile_w]
                            
                            # Save tile
                            tile_filename = f"{file_stem}_green_tile_{tile_idx:04d}_{y}_{x}_{tile_h}x{tile_w}.tif"
                            tile_path = tiles_dir / tile_filename
                            
                            tifffile.imwrite(
                                tile_path, 
                                tile_data, 
                                photometric='minisblack', 
                                compression='zlib'
                            )
                            
                            # Track tile size
                            tile_size_mb = self.get_file_size_mb(tile_path)
                            total_tile_size_mb += tile_size_mb
                            
                            # Record tile information
                            tile_info.append({
                                'tile_id': tile_idx,
                                'filename': tile_filename,
                                'position': (y, x),
                                'size': (tile_h, tile_w),
                                'file_path': str(tile_path),
                                'file_size_mb': tile_size_mb
                            })
                            
                            tile_idx += 1
                            
                            # Log progress every 100 tiles
                            if tile_idx % 100 == 0:
                                self.log_memory_usage(f"PASS2_TILE_{tile_idx}")
                                avg_tile_size = total_tile_size_mb / tile_idx
                                self.logger.info(f"Processed {tile_idx}/{total_tiles} tiles, "
                                              f"avg tile size: {avg_tile_size:.2f} MB")
                            
                        except Exception as e:
                            self.logger.error(f"Error processing tile at ({y},{x}): {e}")
                            continue
                        
                        pbar.update(1)
            
            tiling_time = time.time() - start_time
            
            # Log final statistics
            self.logger.info(f"Tiling completed in {tiling_time:.2f} seconds")
            self.logger.info(f"Average time per tile: {tiling_time/tile_idx:.3f} seconds")
            self.logger.info(f"Total tiles created: {tile_idx}")
            self.logger.info(f"Total tile size: {total_tile_size_mb:.2f} MB")
            self.logger.info(f"Average tile size: {total_tile_size_mb/tile_idx:.2f} MB")
            
            self.log_memory_usage("PASS2_COMPLETE")
            
            # Log final file sizes
            self.log_file_sizes({
                "Green Channel": green_file_path,
                "Tiles Directory": tiles_dir
            }, "PASS2_COMPLETE")
            
            # Save metadata
            metadata = {
                'source_file': str(original_file_path),
                'green_file': str(green_file_path),
                'original_dimensions': (height, width),
                'green_channel': self.config['green_channel'],
                'tile_size': tile_size,
                'total_tiles': tile_idx,
                'total_tile_size_mb': total_tile_size_mb,
                'average_tile_size_mb': total_tile_size_mb / tile_idx,
                'tiling_time_seconds': tiling_time,
                'processing_time': time.time(),
                'timestamp': datetime.now().isoformat(),
                'tiles': tile_info
            }
            
            metadata_file = tiles_dir.parent / f'{file_stem}_metadata.json'
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)
            
            self.logger.info(f"Successfully tiled {green_file_path.name}: {tile_idx} tiles created")
            return metadata
            
        except Exception as e:
            self.logger.error(f"Failed to tile {green_file_path.name}: {e}")
            raise
            
    def process_single_file(self, file_path: Path) -> dict:
        """Process a single BTF file using two-pass approach."""
        start_time = time.time()
        self.logger.info(f"Processing file: {file_path.name}")
        
        # Log initial file size
        original_size_mb = self.get_file_size_mb(file_path)
        self.logger.info(f"Original BTF file size: {original_size_mb:.2f} MB")
        
        try:
            # Pass 1: Extract green channel
            green_file_path = self.extract_green_channel(file_path)
            green_size_mb = self.get_file_size_mb(green_file_path)
            
            # Pass 2: Tile the green channel file
            metadata = self.tile_green_file(green_file_path, file_path)
            
            # Calculate final statistics
            total_tiles = metadata['total_tiles']
            total_tile_size_mb = metadata['total_tile_size_mb']
            avg_tile_size_mb = metadata['average_tile_size_mb']
            
            # Log comprehensive summary
            self.logger.info("=" * 60)
            self.logger.info("PROCESSING SUMMARY")
            self.logger.info("=" * 60)
            self.logger.info(f"Original BTF file: {original_size_mb:.2f} MB")
            self.logger.info(f"Green channel file: {green_size_mb:.2f} MB")
            self.logger.info(f"Total tiles created: {total_tiles}")
            self.logger.info(f"Total tile size: {total_tile_size_mb:.2f} MB")
            self.logger.info(f"Average tile size: {avg_tile_size_mb:.2f} MB")
            self.logger.info(f"Size reduction (BTF → Green): {((original_size_mb - green_size_mb) / original_size_mb) * 100:.1f}%")
            self.logger.info(f"Size reduction (BTF → Tiles): {((original_size_mb - total_tile_size_mb) / original_size_mb) * 100:.1f}%")
            self.logger.info(f"Processing time: {time.time() - start_time:.2f} seconds")
            self.logger.info("=" * 60)
            
            # Clean up intermediate green file if requested
            if self.config.get('cleanup_intermediate', True):
                self.logger.info(f"Cleaning up intermediate file: {green_file_path}")
                green_file_path.unlink()
                self.logger.info("Intermediate green channel file removed")
            
            metadata['total_processing_time'] = time.time() - start_time
            metadata['original_file_size_mb'] = original_size_mb
            metadata['green_file_size_mb'] = green_size_mb
            
            self.logger.info(f"Successfully processed {file_path.name} in {metadata['total_processing_time']:.2f} seconds")
            return metadata
            
        except Exception as e:
            self.logger.error(f"Failed to process {file_path.name}: {e}")
            raise
            
    def run(self):
        """Run the complete two-pass processing pipeline."""
        self.logger.info("Starting two-pass BTF processing pipeline")
        start_time = time.time()
        
        try:
            btf_files = self.get_btf_files()
            results = []
            
            for file_path in btf_files:
                try:
                    result = self.process_single_file(file_path)
                    results.append(result)
                except Exception as e:
                    self.logger.error(f"Failed to process {file_path}: {e}")
                    continue
            
            # Save summary
            summary = {
                'total_files': len(btf_files),
                'successful_files': len(results),
                'failed_files': len(btf_files) - len(results),
                'total_processing_time': time.time() - start_time,
                'timestamp': datetime.now().isoformat(),
                'results': results
            }
            
            summary_file = Path(self.config['output_dir']) / 'two_pass_processing_summary.json'
            with open(summary_file, 'w') as f:
                json.dump(summary, f, indent=2)
            
            self.logger.info(f"Two-pass pipeline completed: {len(results)}/{len(btf_files)} files processed successfully")
            self.logger.info(f"Summary saved to: {summary_file}")
            
        except Exception as e:
            self.logger.error(f"Two-pass pipeline failed: {e}")
            raise


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def create_default_config() -> dict:
    """Create default configuration."""
    return {
        'input_dir': './input',
        'output_dir': './output',
        'tile_size': 2048,
        'green_channel': 1,
        'log_level': 'INFO',
        'log_file': 'two_pass_btf_processing.log',
        'cleanup_intermediate': True
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description='Two-Pass BTF File Processor')
    parser.add_argument('--config', type=str, help='Path to configuration YAML file')
    parser.add_argument('--input-dir', type=str, help='Input directory containing BTF files')
    parser.add_argument('--output-dir', type=str, help='Output directory for processed tiles')
    parser.add_argument('--tile-size', type=int, default=2048, help='Tile size (default: 2048)')
    parser.add_argument('--green-channel', type=int, default=1, help='Green channel index (default: 1)')
    parser.add_argument('--log-level', type=str, default='INFO', 
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='Logging level')
    parser.add_argument('--cleanup-intermediate', type=str, default='true',
                       help='Clean up intermediate green channel files (true/false)')
    
    args = parser.parse_args()
    
    # Load configuration
    if args.config:
        config = load_config(args.config)
    else:
        config = create_default_config()
        
    # Override with command line arguments
    if args.input_dir:
        config['input_dir'] = args.input_dir
    if args.output_dir:
        config['output_dir'] = args.output_dir
    if args.tile_size:
        config['tile_size'] = args.tile_size
    if args.green_channel:
        config['green_channel'] = args.green_channel
    if args.log_level:
        config['log_level'] = args.log_level
    if args.cleanup_intermediate:
        config['cleanup_intermediate'] = args.cleanup_intermediate.lower() == 'true'
        
    # Ensure log file path is absolute
    if not os.path.isabs(config['log_file']):
        config['log_file'] = os.path.join(config['output_dir'], config['log_file'])
    
    # Run processor
    processor = TwoPassBTFProcessor(config)
    processor.run()


if __name__ == '__main__':
    main()
