#!/usr/bin/env python3
"""
Green Channel Extractor
======================

Extracts only the green channel from BTF files without any tiling.
Creates clean TIFF files without OME metadata corruption.

Usage:
    python extract_green_only.py --input-dir /path/to/btf/files --output-dir /path/to/output
    python extract_green_only.py --config config.yaml
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


class GreenChannelExtractor:
    """Extract green channels from BTF files."""
    
    def __init__(self, config: dict):
        """Initialize extractor with configuration."""
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
        required_keys = ['input_dir', 'output_dir', 'green_channel']
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
        """Extract green channel to a new file."""
        self.logger.info(f"Extracting green channel from {file_path.name}")
        
        # Log initial memory and file sizes
        self.log_memory_usage("EXTRACT_START")
        self.log_file_sizes({"Original BTF": file_path}, "EXTRACT_START")
        
        # Create output path for green channel file
        green_file_path = Path(self.config['output_dir']) / f"{file_path.stem}_green.tif"
        
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                
                # Read the entire image (memory intensive but necessary for corrupted OME)
                self.logger.warning("Reading entire image to extract green channel...")
                self.log_memory_usage("EXTRACT_READING_IMAGE")
                
                start_time = time.time()
                full_image = tifffile.imread(file_path, key=0)
                read_time = time.time() - start_time
                
                self.logger.info(f"Image read completed in {read_time:.2f} seconds")
                self.log_memory_usage("EXTRACT_IMAGE_LOADED")
                
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
                self.log_memory_usage("EXTRACT_SAVING_GREEN")
                
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
                self.log_memory_usage("EXTRACT_GREEN_SAVED")
                
                # Log final file sizes
                self.log_file_sizes({
                    "Original BTF": file_path,
                    "Green Channel": green_file_path
                }, "EXTRACT_COMPLETE")
                
                return green_file_path
                
        except Exception as e:
            self.logger.error(f"Failed to extract green channel from {file_path}: {e}")
            raise
            
    def process_single_file(self, file_path: Path) -> dict:
        """Process a single BTF file to extract green channel."""
        start_time = time.time()
        self.logger.info(f"Processing file: {file_path.name}")
        
        # Log initial file size
        original_size_mb = self.get_file_size_mb(file_path)
        self.logger.info(f"Original BTF file size: {original_size_mb:.2f} MB")
        
        try:
            # Extract green channel
            green_file_path = self.extract_green_channel(file_path)
            green_size_mb = self.get_file_size_mb(green_file_path)
            
            # Calculate final statistics
            reduction_percent = ((original_size_mb - green_size_mb) / original_size_mb) * 100
            
            # Log comprehensive summary
            self.logger.info("=" * 60)
            self.logger.info("GREEN CHANNEL EXTRACTION SUMMARY")
            self.logger.info("=" * 60)
            self.logger.info(f"Original BTF file: {original_size_mb:.2f} MB")
            self.logger.info(f"Green channel file: {green_size_mb:.2f} MB")
            self.logger.info(f"Size reduction: {reduction_percent:.1f}%")
            self.logger.info(f"Processing time: {time.time() - start_time:.2f} seconds")
            self.logger.info(f"Output file: {green_file_path}")
            self.logger.info("=" * 60)
            
            metadata = {
                'source_file': str(file_path),
                'green_file': str(green_file_path),
                'original_file_size_mb': original_size_mb,
                'green_file_size_mb': green_size_mb,
                'size_reduction_percent': reduction_percent,
                'processing_time': time.time() - start_time,
                'timestamp': datetime.now().isoformat()
            }
            
            return metadata
            
        except Exception as e:
            self.logger.error(f"Failed to process {file_path.name}: {e}")
            raise
            
    def run(self):
        """Run the green channel extraction pipeline."""
        self.logger.info("Starting green channel extraction pipeline")
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
            
            summary_file = Path(self.config['output_dir']) / 'green_extraction_summary.json'
            with open(summary_file, 'w') as f:
                json.dump(summary, f, indent=2)
            
            self.logger.info(f"Green channel extraction completed: {len(results)}/{len(btf_files)} files processed successfully")
            self.logger.info(f"Summary saved to: {summary_file}")
            
        except Exception as e:
            self.logger.error(f"Green channel extraction pipeline failed: {e}")
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
        'green_channel': 1,
        'log_level': 'INFO',
        'log_file': 'green_extraction.log'
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description='Green Channel Extractor')
    parser.add_argument('--config', type=str, help='Path to configuration YAML file')
    parser.add_argument('--input-dir', type=str, help='Input directory containing BTF files')
    parser.add_argument('--output-dir', type=str, help='Output directory for green channel files')
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
    if args.green_channel:
        config['green_channel'] = args.green_channel
    if args.log_level:
        config['log_level'] = args.log_level
        
    # Ensure log file path is absolute
    if not os.path.isabs(config['log_file']):
        config['log_file'] = os.path.join(config['output_dir'], config['log_file'])
    
    # Run extractor
    extractor = GreenChannelExtractor(config)
    extractor.run()


if __name__ == '__main__':
    main()
