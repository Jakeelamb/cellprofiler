#!/usr/bin/env python3
"""
BTF Pipeline Manager
===================

Utility script for managing the BTF processing pipeline on HPC clusters.
Provides commands for setup, monitoring, and cleanup.

Usage:
    python pipeline_manager.py setup --input-dir /path/to/btf/files
    python pipeline_manager.py submit --config config.yaml
    python pipeline_manager.py status --job-id 12345
    python pipeline_manager.py validate --output-dir /path/to/output
"""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional
import yaml


class PipelineManager:
    """Manages BTF processing pipeline operations."""
    
    def __init__(self):
        self.work_dir = Path.cwd()
        
    def setup_environment(self, input_dir: str, output_dir: str, config_file: Optional[str] = None):
        """Setup the processing environment."""
        print("Setting up BTF processing environment...")
        
        # Create directories
        Path(input_dir).mkdir(parents=True, exist_ok=True)
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Create default config if not provided
        if not config_file:
            config_file = self.work_dir / "config.yaml"
            self.create_default_config(config_file, input_dir, output_dir)
            
        print(f"✅ Environment setup complete")
        print(f"   Input directory: {input_dir}")
        print(f"   Output directory: {output_dir}")
        print(f"   Config file: {config_file}")
        
    def create_default_config(self, config_file: Path, input_dir: str, output_dir: str):
        """Create default configuration file."""
        config = {
            'input_dir': str(input_dir),
            'output_dir': str(output_dir),
            'tile_size': 2048,
            'green_channel': 1,
            'log_level': 'INFO',
            'log_file': 'btf_processing.log'
        }
        
        with open(config_file, 'w') as f:
            yaml.dump(config, f, default_flow_style=False)
            
    def submit_job(self, config_file: str, sbatch_file: str = "run_btf_processing.sbatch"):
        """Submit job to SLURM scheduler."""
        print(f"Submitting job with config: {config_file}")
        
        # Check if config file exists
        if not Path(config_file).exists():
            print(f"❌ Config file not found: {config_file}")
            return None
            
        # Check if sbatch file exists
        if not Path(sbatch_file).exists():
            print(f"❌ SLURM script not found: {sbatch_file}")
            return None
            
        # Submit job
        try:
            result = subprocess.run(
                ['sbatch', sbatch_file],
                capture_output=True,
                text=True,
                check=True
            )
            
            job_id = result.stdout.strip().split()[-1]
            print(f"✅ Job submitted successfully: {job_id}")
            
            # Save job info
            job_info = {
                'job_id': job_id,
                'config_file': config_file,
                'submit_time': time.time(),
                'status': 'PENDING'
            }
            
            with open('job_info.json', 'w') as f:
                json.dump(job_info, f, indent=2)
                
            return job_id
            
        except subprocess.CalledProcessError as e:
            print(f"❌ Failed to submit job: {e}")
            print(f"   Error: {e.stderr}")
            return None
            
    def check_job_status(self, job_id: str) -> Dict:
        """Check status of submitted job."""
        try:
            result = subprocess.run(
                ['squeue', '-j', job_id, '--format=%i,%T,%M,%N', '--noheader'],
                capture_output=True,
                text=True,
                check=True
            )
            
            if result.stdout.strip():
                parts = result.stdout.strip().split(',')
                return {
                    'job_id': parts[0],
                    'status': parts[1],
                    'runtime': parts[2],
                    'nodes': parts[3]
                }
            else:
                return {'job_id': job_id, 'status': 'COMPLETED', 'runtime': 'N/A', 'nodes': 'N/A'}
                
        except subprocess.CalledProcessError:
            return {'job_id': job_id, 'status': 'UNKNOWN', 'runtime': 'N/A', 'nodes': 'N/A'}
            
    def monitor_job(self, job_id: str, interval: int = 60):
        """Monitor job progress."""
        print(f"Monitoring job {job_id} (checking every {interval} seconds)")
        
        while True:
            status = self.check_job_status(job_id)
            print(f"[{time.strftime('%H:%M:%S')}] Job {job_id}: {status['status']} - {status['runtime']}")
            
            if status['status'] in ['COMPLETED', 'FAILED', 'CANCELLED']:
                break
                
            time.sleep(interval)
            
        print(f"Job {job_id} finished with status: {status['status']}")
        return status
        
    def validate_output(self, output_dir: str):
        """Validate processing output."""
        print(f"Validating output in: {output_dir}")
        
        try:
            result = subprocess.run(
                [sys.executable, 'validate_output.py', '--output-dir', output_dir],
                capture_output=True,
                text=True,
                check=True
            )
            
            print("✅ Validation completed successfully")
            print(result.stdout)
            
        except subprocess.CalledProcessError as e:
            print(f"❌ Validation failed: {e}")
            print(f"   Error: {e.stderr}")
            
    def cleanup_failed_jobs(self, days: int = 7):
        """Clean up failed job files older than specified days."""
        print(f"Cleaning up failed job files older than {days} days...")
        
        cutoff_time = time.time() - (days * 24 * 60 * 60)
        cleaned = 0
        
        for log_file in Path('.').glob('btf_processing_*.out'):
            if log_file.stat().st_mtime < cutoff_time:
                log_file.unlink()
                cleaned += 1
                
        for log_file in Path('.').glob('btf_processing_*.err'):
            if log_file.stat().st_mtime < cutoff_time:
                log_file.unlink()
                cleaned += 1
                
        print(f"✅ Cleaned up {cleaned} old log files")
        
    def show_usage_stats(self, output_dir: str):
        """Show disk usage statistics."""
        print(f"Disk usage for {output_dir}:")
        
        try:
            result = subprocess.run(
                ['du', '-sh', output_dir],
                capture_output=True,
                text=True,
                check=True
            )
            print(f"   Total size: {result.stdout.strip().split()[0]}")
            
            # Count files
            tif_files = list(Path(output_dir).rglob('*.tif'))
            json_files = list(Path(output_dir).rglob('*.json'))
            
            print(f"   TIFF files: {len(tif_files)}")
            print(f"   JSON files: {len(json_files)}")
            
        except subprocess.CalledProcessError as e:
            print(f"❌ Failed to get usage stats: {e}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description='BTF Pipeline Manager')
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Setup command
    setup_parser = subparsers.add_parser('setup', help='Setup processing environment')
    setup_parser.add_argument('--input-dir', required=True, help='Input directory for BTF files')
    setup_parser.add_argument('--output-dir', required=True, help='Output directory for processed tiles')
    setup_parser.add_argument('--config', help='Custom config file path')
    
    # Submit command
    submit_parser = subparsers.add_parser('submit', help='Submit job to SLURM')
    submit_parser.add_argument('--config', required=True, help='Configuration file')
    submit_parser.add_argument('--sbatch', default='run_btf_processing.sbatch', help='SLURM script file')
    
    # Status command
    status_parser = subparsers.add_parser('status', help='Check job status')
    status_parser.add_argument('--job-id', required=True, help='Job ID to check')
    
    # Monitor command
    monitor_parser = subparsers.add_parser('monitor', help='Monitor job progress')
    monitor_parser.add_argument('--job-id', required=True, help='Job ID to monitor')
    monitor_parser.add_argument('--interval', type=int, default=60, help='Check interval in seconds')
    
    # Validate command
    validate_parser = subparsers.add_parser('validate', help='Validate processing output')
    validate_parser.add_argument('--output-dir', required=True, help='Output directory to validate')
    
    # Cleanup command
    cleanup_parser = subparsers.add_parser('cleanup', help='Clean up old files')
    cleanup_parser.add_argument('--days', type=int, default=7, help='Days to keep files')
    
    # Stats command
    stats_parser = subparsers.add_parser('stats', help='Show usage statistics')
    stats_parser.add_argument('--output-dir', required=True, help='Output directory to analyze')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return
        
    manager = PipelineManager()
    
    if args.command == 'setup':
        manager.setup_environment(args.input_dir, args.output_dir, args.config)
    elif args.command == 'submit':
        manager.submit_job(args.config, args.sbatch)
    elif args.command == 'status':
        status = manager.check_job_status(args.job_id)
        print(f"Job {args.job_id}: {status['status']} - {status['runtime']}")
    elif args.command == 'monitor':
        manager.monitor_job(args.job_id, args.interval)
    elif args.command == 'validate':
        manager.validate_output(args.output_dir)
    elif args.command == 'cleanup':
        manager.cleanup_failed_jobs(args.days)
    elif args.command == 'stats':
        manager.show_usage_stats(args.output_dir)


if __name__ == '__main__':
    main()
