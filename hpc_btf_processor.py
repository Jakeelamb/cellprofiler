#!/usr/bin/env python3
"""
HPC BTF File Processor
=====================

A robust pipeline for processing large BTF files on HPC clusters.
Extracts green channels and creates memory-efficient tiles for CellProfiler analysis.

Usage:
    python hpc_btf_processor.py --input-dir /path/to/btf/files --output-dir /path/to/output
    python hpc_btf_processor.py --config config.yaml
"""

import argparse
import logging
import os
import sys
import time
from pathlib import Path
from typing import List, Tuple, Optional
import yaml
import tifffile
import numpy as np
from tqdm import tqdm
import json
from datetime import datetime


class BTFProcessor:
    """Main processor class for BTF files on HPC clusters."""
    
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
        
    def get_btf_files(self) -> List[Path]:
        """Get list of BTF files to process."""
        input_dir = Path(self.config['input_dir'])
        btf_files = list(input_dir.glob('*.btf')) + list(input_dir.glob('*.ome.btf'))
        
        if not btf_files:
            raise FileNotFoundError(f"No BTF files found in {input_dir}")
            
        self.logger.info(f"Found {len(btf_files)} BTF files to process")
        return sorted(btf_files)
        
    def get_image_info(self, file_path: Path) -> Tuple[int, int, int, str]:
        """Extract image dimensions and channel information."""
        try:
            # Suppress OME warnings for corrupted metadata
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                
                with tifffile.TiffFile(file_path) as tif:
                    # Try to get info from series first
                    try:
                        img = tif.series[0]
                        
                        if 'S' in img.axes:  # RGB (YXS axes)
                            num_channels = img.shape[img.axes.index('S')]
                            height, width = img.shape[img.axes.index('Y')], img.shape[img.axes.index('X')]
                            axes = 'YXS'
                        elif 'C' in img.axes:  # CYX format
                            num_channels = img.shape[img.axes.index('C')]
                            height, width = img.shape[img.axes.index('Y')], img.shape[img.axes.index('X')]
                            axes = 'CYX'
                        else:  # Grayscale
                            num_channels = 1
                            height, width = img.shape[-2], img.shape[-1]
                            axes = 'YX'
                            
                        return height, width, num_channels, axes
                        
                    except Exception as series_error:
                        # Fallback: read from first page directly
                        self.logger.warning(f"Series reading failed, using direct page reading: {series_error}")
                        page = tif.pages[0]
                        page_data = page.asarray()
                        
                        if len(page_data.shape) == 3:  # RGB format
                            height, width, num_channels = page_data.shape
                            axes = 'YXS'
                        elif len(page_data.shape) == 2:  # Grayscale
                            height, width = page_data.shape
                            num_channels = 1
                            axes = 'YX'
                        else:  # Unexpected format
                            height, width = page_data.shape[-2], page_data.shape[-1]
                            num_channels = 1
                            axes = 'YX'
                            
                        return height, width, num_channels, axes
                
        except Exception as e:
            self.logger.error(f"Error reading image info from {file_path}: {e}")
            raise
            
    def create_output_structure(self, file_path: Path) -> Path:
        """Create organized output directory structure."""
        file_stem = file_path.stem.replace('.ome', '')  # Remove .ome if present
        output_dir = Path(self.config['output_dir']) / file_stem / 'green_tiles'
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir
        
    def process_single_file(self, file_path: Path) -> dict:
        """Process a single BTF file."""
        start_time = time.time()
        self.logger.info(f"Processing file: {file_path.name}")
        
        try:
            # Get image information
            height, width, num_channels, axes = self.get_image_info(file_path)
            self.logger.info(f"Image dimensions: {height}x{width}, Channels: {num_channels}, Axes: {axes}")
            
            # Create output directory
            output_dir = self.create_output_structure(file_path)
            
            # Determine green channel index
            green_channel = self.config['green_channel']
            if green_channel >= num_channels:
                self.logger.warning(f"Green channel {green_channel} not available, using channel 0")
                green_channel = 0
                
            # Process tiles
            tile_size = self.config['tile_size']
            total_tiles = ((height + tile_size - 1) // tile_size) * ((width + tile_size - 1) // tile_size)
            
            self.logger.info(f"Creating {total_tiles} tiles of size {tile_size}x{tile_size}")
            
            tile_info = []
            tile_idx = 0
            
            with tqdm(total=total_tiles, desc=f"Processing {file_path.name}") as pbar:
                for y in range(0, height, tile_size):
                    for x in range(0, width, tile_size):
                        tile_h = min(tile_size, height - y)
                        tile_w = min(tile_size, width - x)
                        
                        try:
                            # Read tile data
                            tile_data = self.extract_tile(
                                file_path, y, x, tile_h, tile_w, 
                                green_channel, num_channels, axes
                            )
                            
                            # Save tile
                            tile_filename = self.create_tile_filename(
                                file_path, tile_idx, y, x, tile_h, tile_w
                            )
                            tile_path = output_dir / tile_filename
                            
                            tifffile.imwrite(
                                tile_path, 
                                tile_data, 
                                photometric='minisblack', 
                                compression='zlib'
                            )
                            
                            # Record tile information
                            tile_info.append({
                                'tile_id': tile_idx,
                                'filename': tile_filename,
                                'position': (y, x),
                                'size': (tile_h, tile_w),
                                'file_path': str(tile_path)
                            })
                            
                            tile_idx += 1
                            
                        except Exception as e:
                            self.logger.error(f"Error processing tile at ({y},{x}): {e}")
                            continue
                        
                        pbar.update(1)
            
            # Save metadata
            metadata = {
                'source_file': str(file_path),
                'original_dimensions': (height, width),
                'num_channels': num_channels,
                'axes': axes,
                'green_channel': green_channel,
                'tile_size': tile_size,
                'total_tiles': tile_idx,
                'processing_time': time.time() - start_time,
                'timestamp': datetime.now().isoformat(),
                'tiles': tile_info
            }
            
            metadata_file = output_dir.parent / f'{file_path.stem}_metadata.json'
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)
            
            self.logger.info(f"Successfully processed {file_path.name}: {tile_idx} tiles created")
            return metadata
            
        except Exception as e:
            self.logger.error(f"Failed to process {file_path.name}: {e}")
            raise
            
    def extract_tile(self, file_path: Path, y: int, x: int, h: int, w: int, 
                    green_channel: int, num_channels: int, axes: str) -> np.ndarray:
        """Extract a single tile from the BTF file."""
        try:
            # Try the memory-efficient approach first
            with tifffile.TiffFile(file_path) as tif:
                img = tif.series[0]
                
                # Use memory-efficient key-based slicing to avoid loading entire image
                if axes == 'YXS':  # RGB format
                    tile_region = img.asarray(
                        key=(slice(y, y+h), slice(x, x+w), slice(None))
                    )
                    return tile_region[:, :, green_channel]
                elif axes == 'CYX':  # Channel first format
                    tile_region = img.asarray(
                        key=(slice(None), slice(y, y+h), slice(x, x+w))
                    )
                    return tile_region[green_channel, :, :]
                else:  # Grayscale
                    return img.asarray(
                        key=(slice(y, y+h), slice(x, x+w))
                    )
        except Exception as e:
            # Fallback: use direct page reading for corrupted OME files
            self.logger.warning(f"Key-based slicing failed, using direct page reading: {e}")
            return self._extract_tile_direct_pages(file_path, y, x, h, w, green_channel, num_channels, axes)
    
    def _extract_tile_direct_pages(self, file_path: Path, y: int, x: int, h: int, w: int, 
                                  green_channel: int, num_channels: int, axes: str) -> np.ndarray:
        """Extract tile using direct page reading for corrupted OME files."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            with tifffile.TiffFile(file_path) as tif:
                # Read the first page to get basic info
                page = tif.pages[0]
                page_data = page.asarray()
                
                # Calculate which pages we need to read
                if axes == 'YXS':  # RGB format - single page with all channels
                    # Extract the tile region and green channel
                    tile_region = page_data[y:y+h, x:x+w, :]
                    return tile_region[:, :, green_channel]
                elif axes == 'CYX':  # Channel first format - multiple pages
                    # For CYX format, we need to read the specific channel page
                    if green_channel < len(tif.pages):
                        channel_page = tif.pages[green_channel]
                        channel_data = channel_page.asarray()
                        return channel_data[y:y+h, x:x+w]
                    else:
                        # Fallback: read first page and assume it's the green channel
                        return page_data[y:y+h, x:x+w]
                else:  # Grayscale
                    return page_data[y:y+h, x:x+w]
                
    def create_tile_filename(self, file_path: Path, tile_idx: int, y: int, x: int, 
                           h: int, w: int) -> str:
        """Create standardized tile filename."""
        file_stem = file_path.stem.replace('.ome', '')
        return f"{file_stem}_green_tile_{tile_idx:04d}_{y}_{x}_{h}x{w}.tif"
        
    def run(self):
        """Run the complete processing pipeline."""
        self.logger.info("Starting BTF processing pipeline")
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
            
            summary_file = Path(self.config['output_dir']) / 'processing_summary.json'
            with open(summary_file, 'w') as f:
                json.dump(summary, f, indent=2)
            
            self.logger.info(f"Pipeline completed: {len(results)}/{len(btf_files)} files processed successfully")
            self.logger.info(f"Summary saved to: {summary_file}")
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {e}")
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
        'log_file': 'btf_processing.log'
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description='HPC BTF File Processor')
    parser.add_argument('--config', type=str, help='Path to configuration YAML file')
    parser.add_argument('--input-dir', type=str, help='Input directory containing BTF files')
    parser.add_argument('--output-dir', type=str, help='Output directory for processed tiles')
    parser.add_argument('--tile-size', type=int, default=2048, help='Tile size (default: 2048)')
    parser.add_argument('--green-channel', type=int, default=1, help='Green channel index (default: 1)')
    parser.add_argument('--log-level', type=str, default='INFO', 
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='Logging level')
    
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
        
    # Ensure log file path is absolute
    if not os.path.isabs(config['log_file']):
        config['log_file'] = os.path.join(config['output_dir'], config['log_file'])
    
    # Run processor
    processor = BTFProcessor(config)
    processor.run()


if __name__ == '__main__':
    main()
