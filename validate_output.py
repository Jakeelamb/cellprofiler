#!/usr/bin/env python3
"""
Output Validation Tool
======================

Validates the output of BTF processing pipeline to ensure:
- All tiles are readable
- File naming convention is correct
- Metadata is consistent
- No data corruption

Usage:
    python validate_output.py --output-dir /path/to/output
    python validate_output.py --config config.yaml
"""

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple
import tifffile
import numpy as np
from tqdm import tqdm


class OutputValidator:
    """Validates BTF processing output."""
    
    def __init__(self, output_dir: str, log_level: str = 'INFO'):
        self.output_dir = Path(output_dir)
        self.setup_logging(log_level)
        
    def setup_logging(self, log_level: str):
        """Setup logging configuration."""
        level = getattr(logging, log_level.upper())
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.output_dir / 'validation.log'),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def find_metadata_files(self) -> List[Path]:
        """Find all metadata JSON files."""
        return list(self.output_dir.rglob('*_metadata.json'))
        
    def validate_metadata_file(self, metadata_path: Path) -> Dict:
        """Validate a single metadata file."""
        try:
            with open(metadata_path, 'r') as f:
                metadata = json.load(f)
                
            # Check required fields
            required_fields = [
                'source_file', 'original_dimensions', 'num_channels',
                'tile_size', 'total_tiles', 'tiles'
            ]
            
            missing_fields = [field for field in required_fields if field not in metadata]
            if missing_fields:
                return {
                    'valid': False,
                    'errors': [f"Missing required fields: {missing_fields}"]
                }
                
            # Validate tile information
            tile_errors = []
            for i, tile in enumerate(metadata['tiles']):
                tile_path = Path(tile['file_path'])
                if not tile_path.exists():
                    tile_errors.append(f"Tile {i}: File not found - {tile_path}")
                    
            return {
                'valid': len(tile_errors) == 0,
                'errors': tile_errors,
                'metadata': metadata
            }
            
        except Exception as e:
            return {
                'valid': False,
                'errors': [f"Error reading metadata: {e}"]
            }
            
    def validate_tile_file(self, tile_path: Path) -> Dict:
        """Validate a single tile file."""
        try:
            # Check if file exists and is readable
            if not tile_path.exists():
                return {'valid': False, 'error': 'File does not exist'}
                
            # Try to read the image
            with tifffile.TiffFile(tile_path) as tif:
                img = tif.series[0]
                shape = img.shape
                dtype = img.dtype
                
            # Basic validation
            if len(shape) != 2:
                return {'valid': False, 'error': f'Expected 2D image, got {len(shape)}D'}
                
            if dtype not in [np.uint8, np.uint16, np.uint32, np.float32, np.float64]:
                return {'valid': False, 'error': f'Unsupported data type: {dtype}'}
                
            # Check for reasonable dimensions
            if shape[0] == 0 or shape[1] == 0:
                return {'valid': False, 'error': 'Empty image'}
                
            return {
                'valid': True,
                'shape': shape,
                'dtype': str(dtype)
            }
            
        except Exception as e:
            return {'valid': False, 'error': f'Error reading file: {e}'}
            
    def validate_naming_convention(self, tile_path: Path) -> Dict:
        """Validate tile naming convention."""
        filename = tile_path.name
        
        # Expected format: {source}_green_tile_{tile_id:04d}_{y}_{x}_{h}x{w}.tif
        parts = filename.replace('.tif', '').split('_')
        
        if len(parts) < 6:
            return {'valid': False, 'error': 'Filename does not match expected format'}
            
        if parts[-2] != 'green':
            return {'valid': False, 'error': 'Filename missing "green" identifier'}
            
        if parts[-3] != 'tile':
            return {'valid': False, 'error': 'Filename missing "tile" identifier'}
            
        try:
            tile_id = int(parts[-4])
            y = int(parts[-3])
            x = int(parts[-2])
            dimensions = parts[-1].split('x')
            h, w = int(dimensions[0]), int(dimensions[1])
            
            return {
                'valid': True,
                'tile_id': tile_id,
                'position': (y, x),
                'dimensions': (h, w)
            }
            
        except (ValueError, IndexError) as e:
            return {'valid': False, 'error': f'Error parsing filename: {e}'}
            
    def run_validation(self) -> Dict:
        """Run complete validation."""
        self.logger.info("Starting output validation")
        
        results = {
            'total_metadata_files': 0,
            'valid_metadata_files': 0,
            'total_tiles': 0,
            'valid_tiles': 0,
            'naming_errors': 0,
            'file_errors': 0,
            'summary': {}
        }
        
        # Find and validate metadata files
        metadata_files = self.find_metadata_files()
        results['total_metadata_files'] = len(metadata_files)
        
        self.logger.info(f"Found {len(metadata_files)} metadata files")
        
        for metadata_path in tqdm(metadata_files, desc="Validating metadata"):
            validation_result = self.validate_metadata_file(metadata_path)
            if validation_result['valid']:
                results['valid_metadata_files'] += 1
            else:
                self.logger.error(f"Metadata validation failed for {metadata_path}: {validation_result['errors']}")
                
        # Find and validate tile files
        tile_files = list(self.output_dir.rglob('*.tif'))
        results['total_tiles'] = len(tile_files)
        
        self.logger.info(f"Found {len(tile_files)} tile files")
        
        for tile_path in tqdm(tile_files, desc="Validating tiles"):
            # Validate file
            file_validation = self.validate_tile_file(tile_path)
            if file_validation['valid']:
                results['valid_tiles'] += 1
            else:
                results['file_errors'] += 1
                self.logger.error(f"Tile validation failed for {tile_path}: {file_validation['error']}")
                
            # Validate naming
            naming_validation = self.validate_naming_convention(tile_path)
            if not naming_validation['valid']:
                results['naming_errors'] += 1
                self.logger.error(f"Naming validation failed for {tile_path}: {naming_validation['error']}")
                
        # Calculate summary
        results['summary'] = {
            'metadata_success_rate': results['valid_metadata_files'] / max(results['total_metadata_files'], 1),
            'tile_success_rate': results['valid_tiles'] / max(results['total_tiles'], 1),
            'overall_success': (results['valid_metadata_files'] == results['total_metadata_files'] and 
                              results['valid_tiles'] == results['total_tiles'] and 
                              results['naming_errors'] == 0)
        }
        
        self.logger.info("Validation completed")
        self.logger.info(f"Metadata files: {results['valid_metadata_files']}/{results['total_metadata_files']} valid")
        self.logger.info(f"Tile files: {results['valid_tiles']}/{results['total_tiles']} valid")
        self.logger.info(f"Naming errors: {results['naming_errors']}")
        self.logger.info(f"File errors: {results['file_errors']}")
        
        return results


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description='Validate BTF processing output')
    parser.add_argument('--output-dir', type=str, required=True, help='Output directory to validate')
    parser.add_argument('--log-level', type=str, default='INFO', 
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='Logging level')
    
    args = parser.parse_args()
    
    validator = OutputValidator(args.output_dir, args.log_level)
    results = validator.run_validation()
    
    # Save results
    results_file = Path(args.output_dir) / 'validation_results.json'
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nValidation results saved to: {results_file}")
    
    if results['summary']['overall_success']:
        print("✅ All validations passed!")
        sys.exit(0)
    else:
        print("❌ Some validations failed. Check logs for details.")
        sys.exit(1)


if __name__ == '__main__':
    main()
