#!/usr/bin/env python3
"""
Memory Monitoring Script
========================

Monitors memory usage during BTF processing and logs to a separate file.
Can be run alongside the main processing script.

Usage:
    python monitor_memory.py --pid <process_id> --log-file memory_monitor.log
    python monitor_memory.py --process-name hpc_btf_processor_2pass.py
"""

import argparse
import psutil
import time
import logging
import sys
from datetime import datetime
from pathlib import Path


class MemoryMonitor:
    """Monitor memory usage of a specific process."""
    
    def __init__(self, log_file: str = "memory_monitor.log"):
        self.setup_logging(log_file)
        
    def setup_logging(self, log_file: str):
        """Setup logging for memory monitoring."""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def get_process_memory(self, process):
        """Get memory information for a process."""
        try:
            memory_info = process.memory_info()
            return {
                'rss_mb': memory_info.rss / (1024 * 1024),
                'vms_mb': memory_info.vms / (1024 * 1024),
                'percent': process.memory_percent(),
                'cpu_percent': process.cpu_percent(),
                'num_threads': process.num_threads(),
                'status': process.status()
            }
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            return None
            
    def get_system_memory(self):
        """Get system-wide memory information."""
        memory = psutil.virtual_memory()
        return {
            'total_gb': memory.total / (1024**3),
            'available_gb': memory.available / (1024**3),
            'used_gb': memory.used / (1024**3),
            'percent': memory.percent
        }
        
    def monitor_process_by_pid(self, pid: int, interval: float = 5.0):
        """Monitor a process by PID."""
        try:
            process = psutil.Process(pid)
            self.logger.info(f"Monitoring process PID {pid}: {process.name()}")
        except psutil.NoSuchProcess:
            self.logger.error(f"Process with PID {pid} not found")
            return
            
        self.monitor_process(process, interval)
        
    def monitor_process_by_name(self, process_name: str, interval: float = 5.0):
        """Monitor a process by name."""
        processes = []
        for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
            try:
                if process_name in proc.info['name'] or any(process_name in arg for arg in proc.info['cmdline']):
                    processes.append(proc)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
                
        if not processes:
            self.logger.error(f"No processes found matching '{process_name}'")
            return
            
        if len(processes) > 1:
            self.logger.warning(f"Found {len(processes)} processes matching '{process_name}', monitoring all")
            
        for process in processes:
            self.logger.info(f"Monitoring process PID {process.pid}: {process.name()}")
            self.monitor_process(process, interval)
            
    def monitor_process(self, process, interval: float):
        """Monitor a single process."""
        start_time = time.time()
        
        self.logger.info("=" * 80)
        self.logger.info("MEMORY MONITORING STARTED")
        self.logger.info("=" * 80)
        
        try:
            while True:
                # Get process memory
                proc_memory = self.get_process_memory(process)
                if proc_memory is None:
                    self.logger.warning("Process no longer exists or access denied")
                    break
                    
                # Get system memory
                sys_memory = self.get_system_memory()
                
                # Calculate elapsed time
                elapsed = time.time() - start_time
                
                # Log memory information
                self.logger.info(f"[{elapsed:.1f}s] Process Memory: "
                               f"RSS={proc_memory['rss_mb']:.1f}MB, "
                               f"VMS={proc_memory['vms_mb']:.1f}MB, "
                               f"CPU={proc_memory['cpu_percent']:.1f}%, "
                               f"Threads={proc_memory['num_threads']}")
                
                self.logger.info(f"[{elapsed:.1f}s] System Memory: "
                               f"Used={sys_memory['used_gb']:.1f}GB/{sys_memory['total_gb']:.1f}GB "
                               f"({sys_memory['percent']:.1f}%), "
                               f"Available={sys_memory['available_gb']:.1f}GB")
                
                # Check for memory warnings
                if proc_memory['rss_mb'] > 10000:  # > 10GB
                    self.logger.warning(f"HIGH MEMORY USAGE: {proc_memory['rss_mb']:.1f}MB")
                    
                if sys_memory['percent'] > 90:
                    self.logger.warning(f"HIGH SYSTEM MEMORY USAGE: {sys_memory['percent']:.1f}%")
                
                time.sleep(interval)
                
        except KeyboardInterrupt:
            self.logger.info("Memory monitoring stopped by user")
        except Exception as e:
            self.logger.error(f"Error during monitoring: {e}")
        finally:
            self.logger.info("=" * 80)
            self.logger.info("MEMORY MONITORING ENDED")
            self.logger.info("=" * 80)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description='Monitor memory usage during BTF processing')
    parser.add_argument('--pid', type=int, help='Process ID to monitor')
    parser.add_argument('--process-name', type=str, help='Process name to monitor')
    parser.add_argument('--log-file', type=str, default='memory_monitor.log', help='Log file path')
    parser.add_argument('--interval', type=float, default=5.0, help='Monitoring interval in seconds')
    
    args = parser.parse_args()
    
    if not args.pid and not args.process_name:
        parser.error("Must specify either --pid or --process-name")
        
    monitor = MemoryMonitor(args.log_file)
    
    if args.pid:
        monitor.monitor_process_by_pid(args.pid, args.interval)
    else:
        monitor.monitor_process_by_name(args.process_name, args.interval)


if __name__ == '__main__':
    main()
