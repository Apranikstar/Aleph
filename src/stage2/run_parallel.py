#!/usr/bin/env python3
"""
Parallel processing pipeline for stage2 conversion.
Processes multiple ROOT files in parallel, splitting each into chunks.
"""

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple
import yaml
import uproot
from math import ceil
from dataclasses import dataclass
from datetime import datetime


@dataclass
class ProcessingConfig:
    """Configuration for parallel processing."""
    input_dir: str
    output_dir: str
    stage2_script: str = "stage2.py"
    config_file: str = "config.yaml"
    n_workers: int = 64
    delete_parts_after_merge: bool = True
    n_final_files: int = 1
    chunk_size: int = None  # If None, auto-calculate based on n_workers


def get_num_entries(root_file: str, tree_name: str = "events") -> int:
    """Get number of entries in ROOT file."""
    try:
        with uproot.open(root_file) as f:
            tree = f[tree_name]
            return tree.num_entries
    except Exception as e:
        print(f"Warning: Could not read {root_file}: {e}")
        return 0


def build_jobs(config: ProcessingConfig) -> List[Tuple[List[str], str]]:
    """
    Build list of processing jobs.
    
    Returns:
        List of (command, label) tuples
    """
    jobs = []
    input_path = Path(config.input_dir)
    output_path = Path(config.output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find all ROOT files
    root_files = sorted(input_path.glob("*.root"))
    
    if not root_files:
        print(f"WARNING: No ROOT files found in {config.input_dir}")
        return jobs
    
    print(f"Found {len(root_files)} ROOT file(s) to process")
    
    for root_file in root_files:
        filename = root_file.name
        base = root_file.stem
        
        # Get number of entries
        n_entries = get_num_entries(str(root_file))
        if n_entries == 0:
            print(f"Skipping {filename} (no entries)")
            continue
        
        # Calculate chunk size
        if config.chunk_size is not None:
            chunk_size = config.chunk_size
            n_chunks = ceil(n_entries / chunk_size)
        else:
            n_chunks = config.n_workers
            chunk_size = ceil(n_entries / n_chunks)
        
        print(f"  {filename}: {n_entries} entries → {n_chunks} chunks of ~{chunk_size} entries")
        
        # Create job for each chunk
        for i in range(n_chunks):
            start = i * chunk_size
            end = min((i + 1) * chunk_size, n_entries)
            
            if start >= end:
                break
            
            output_file = output_path / f"{base}_part{i:03d}.root"
            
            cmd = [
                "python",
                config.stage2_script,
                str(root_file),
                str(output_file),
                str(start),
                str(end),
                "--config",
                config.config_file
            ]
            
            jobs.append((cmd, f"{base}_part{i:03d}"))
    
    return jobs


def run_job(job: Tuple[List[str], str]) -> Tuple[str, int, str, str]:
    """
    Run a single processing job.
    
    Returns:
        Tuple of (label, returncode, stdout, stderr)
    """
    cmd, label = job
    
    start_time = datetime.now()
    result = subprocess.run(cmd, capture_output=True, text=True)
    duration = (datetime.now() - start_time).total_seconds()
    
    status = "✓" if result.returncode == 0 else "✗"
    print(f"[{status}] {label} ({duration:.1f}s)")
    
    return label, result.returncode, result.stdout, result.stderr


def merge_outputs(config: ProcessingConfig) -> None:
    """Merge part files into final output files."""
    print("\n" + "="*70)
    print("Merging output files")
    print("="*70)
    
    input_path = Path(config.input_dir)
    output_path = Path(config.output_dir)
    
    # Group part files by base name
    root_files = sorted(input_path.glob("*.root"))
    
    for root_file in root_files:
        base = root_file.stem
        part_pattern = f"{base}_part*.root"
        part_files = sorted(output_path.glob(part_pattern))
        
        if not part_files:
            print(f"  No parts found for {base}")
            continue
        
        print(f"\n  Processing {base}: {len(part_files)} parts")
        
        # Split into n_final_files groups
        chunk_size = ceil(len(part_files) / config.n_final_files)
        
        for i in range(config.n_final_files):
            chunk_files = part_files[i*chunk_size:(i+1)*chunk_size]
            
            if not chunk_files:
                continue
            
            # Determine output filename
            if config.n_final_files == 1:
                merged_file = output_path / f"{base}_s2.root"
            else:
                merged_file = output_path / f"{base}_s2_{i}.root"
            
            print(f"    Merging {len(chunk_files)} parts → {merged_file.name}")
            
            # Run hadd
            cmd = ["hadd", "-f", str(merged_file)] + [str(f) for f in chunk_files]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"    ERROR: hadd failed")
                print(f"    {result.stderr}")
            else:
                print(f"    ✓ Merged successfully")
                
                # Delete part files if requested
                if config.delete_parts_after_merge:
                    for part_file in chunk_files:
                        try:
                            part_file.unlink()
                        except Exception as e:
                            print(f"    Warning: Could not delete {part_file.name}: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Parallel processing pipeline for stage2 conversion"
    )
    parser.add_argument(
        "input_dir",
        type=str,
        help="Directory containing input ROOT files"
    )
    parser.add_argument(
        "output_dir",
        type=str,
        help="Directory for output ROOT files"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=64,
        help="Number of parallel workers (default: 64)"
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config.yaml",
        help="Path to configuration file (default: config.yaml)"
    )
    parser.add_argument(
        "--stage2-script",
        type=str,
        default="stage2.py",
        help="Path to stage2 processing script (default: stage2.py)"
    )
    parser.add_argument(
        "--keep-parts",
        action="store_true",
        help="Keep part files after merging (default: delete them)"
    )
    parser.add_argument(
        "--n-final-files",
        type=int,
        default=1,
        help="Number of final merged files per sample (default: 1)"
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=None,
        help="Number of entries per chunk (default: auto-calculate from workers)"
    )
    
    args = parser.parse_args()
    
    # Create processing configuration
    proc_config = ProcessingConfig(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        stage2_script=args.stage2_script,
        config_file=args.config,
        n_workers=args.workers,
        delete_parts_after_merge=not args.keep_parts,
        n_final_files=args.n_final_files,
        chunk_size=args.chunk_size
    )
    
    # Print configuration
    print("="*70)
    print("Stage2 Parallel Processing Pipeline")
    print("="*70)
    print(f"Input directory:  {proc_config.input_dir}")
    print(f"Output directory: {proc_config.output_dir}")
    print(f"Workers:          {proc_config.n_workers}")
    print(f"Config file:      {proc_config.config_file}")
    print(f"Stage2 script:    {proc_config.stage2_script}")
    print(f"Final files:      {proc_config.n_final_files} per sample")
    print(f"Keep parts:       {not proc_config.delete_parts_after_merge}")
    print("="*70)
    print()
    
    # Build job list
    jobs = build_jobs(proc_config)
    
    if not jobs:
        print("No jobs to process!")
        sys.exit(1)
    
    print(f"\nPrepared {len(jobs)} processing jobs")
    print("="*70)
    print()
    
    # Run jobs in parallel
    failed_jobs = []
    
    with ThreadPoolExecutor(max_workers=proc_config.n_workers) as executor:
        futures = {executor.submit(run_job, job): job for job in jobs}
        
        for future in as_completed(futures):
            job = futures[future]
            try:
                label, returncode, stdout, stderr = future.result()
                
                if returncode != 0:
                    failed_jobs.append((label, stderr))
                    
            except Exception as e:
                cmd, label = job
                print(f"[✗] {label} - Exception: {e}")
                failed_jobs.append((label, str(e)))
    
    # Print summary
    print("\n" + "="*70)
    print("Processing Summary")
    print("="*70)
    print(f"Total jobs:    {len(jobs)}")
    print(f"Successful:    {len(jobs) - len(failed_jobs)}")
    print(f"Failed:        {len(failed_jobs)}")
    
    if failed_jobs:
        print("\nFailed jobs:")
        for label, error in failed_jobs:
            print(f"  - {label}")
            if error:
                print(f"    {error[:200]}")
    
    # Merge outputs
    if len(failed_jobs) < len(jobs):  # Only merge if we have some successful jobs
        merge_outputs(proc_config)
    
    print("\n" + "="*70)
    print("All done!")
    print("="*70)
    
    # Exit with error code if any jobs failed
    sys.exit(1 if failed_jobs else 0)


if __name__ == "__main__":
    main()
