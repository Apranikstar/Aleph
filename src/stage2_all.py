#!/usr/bin/env python3
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
import uproot
import glob
from math import ceil
import re

# --- user configuration ---
input_dir = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/1994/zqq/stage1/1.0.0"
output_dir = "./output"
stage2_script = "stage2.py"
ncpus = 64  # Use all available threads
delete_parts_after_merge = True
n_final_files = 1  # Number of final merged output files per sample
# ---------------------------

os.makedirs(output_dir, exist_ok=True)

def get_nentries(root_file):
    """Return number of entries in the first TTree of the ROOT file."""
    try:
        with uproot.open(root_file) as f:
            tree_name = next((k for k in f.keys() if hasattr(f[k], "num_entries")), list(f.keys())[0])
            tree = f[tree_name]
            return tree.num_entries
    except Exception as e:
        print(f"Error reading {root_file}: {e}")
        return 0

def natural_sort_key(filename):
    """Extract the part number for proper numerical sorting."""
    match = re.search(r'_part(\d+)\.root$', filename)
    if match:
        return int(match.group(1))
    return 0

def build_jobs_for_sample(filename):
    """Prepare jobs for a single input ROOT file."""
    jobs = []
    if not filename.endswith(".root"):
        return jobs
    
    input_path = os.path.join(input_dir, filename)
    base = os.path.splitext(filename)[0]
    n_entries = get_nentries(input_path)
    
    if n_entries == 0:
        print(f"Skipping {filename} - no entries found")
        return jobs
        
    chunk_size = (n_entries + ncpus - 1) // ncpus
    
    for i in range(ncpus):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, n_entries)
        if start >= end:
            break
        output_file = os.path.join(output_dir, f"{base}_part{i}.root")
        cmd = ["python", stage2_script, input_path, output_file, str(start), str(end)]
        jobs.append((cmd, f"{base}_part{i}", base))
    return jobs

def run_job(job):
    """Run one job command."""
    cmd, label, sample_name = job
    print(f"[START] {label}: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            print(f"[ERROR] {label}:\n{result.stderr.strip()}")
        else:
            print(f"[DONE] {label}")
        return result.returncode, sample_name
    except subprocess.TimeoutExpired:
        print(f"[TIMEOUT] {label} - job took longer than 10 minutes")
        return -1, sample_name
    except Exception as e:
        print(f"[EXCEPTION] {label}: {e}")
        return -1, sample_name

def merge_sample(base):
    """Merge parts for a single sample."""
    part_files = glob.glob(os.path.join(output_dir, f"{base}_part*.root"))
    part_files = sorted(part_files, key=natural_sort_key)
    
    if not part_files:
        return
    
    print(f"\n{'='*60}")
    print(f"Merging sample: {base}")
    print(f"Found {len(part_files)} part files")
    print(f"Sequential order: {', '.join([f'part{natural_sort_key(f)}' for f in part_files])}")
    
    # Split part_files into n_final_files chunks
    chunk_size = ceil(len(part_files) / n_final_files)
    
    for i in range(n_final_files):
        chunk_files = part_files[i*chunk_size:(i+1)*chunk_size]
        if not chunk_files:
            continue
        
        merged_file = os.path.join(output_dir, f"{base}_s2_{i}.root")
        part_range = f"part{natural_sort_key(chunk_files[0])}-part{natural_sort_key(chunk_files[-1])}"
        
        print(f"\n  Chunk {i}: Merging {len(chunk_files)} files ({part_range})")
        print(f"  Output: {os.path.basename(merged_file)}")
        print(f"  Files: {' → '.join([os.path.basename(f) for f in chunk_files])}")
        
        cmd = ["hadd", "-f", merged_file] + chunk_files
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"  [ERROR] Failed to merge {base} chunk {i}:")
            print(f"  {result.stderr.strip()}")
        else:
            print(f"  [OK] Chunk {i} merged successfully!")
            if delete_parts_after_merge:
                for f in chunk_files:
                    try:
                        os.remove(f)
                    except Exception as e:
                        print(f"  Warning: could not remove {os.path.basename(f)}: {e}")

def main():
    all_files = [f for f in os.listdir(input_dir) if f.endswith(".root")]
    total_samples = len(all_files)
    
    print(f"{'='*60}")
    print(f"Configuration:")
    print(f"  Input directory: {input_dir}")
    print(f"  Output directory: {output_dir}")
    print(f"  Parallel workers: {ncpus}")
    print(f"  Final files per sample: {n_final_files}")
    print(f"  Total samples: {total_samples}")
    print(f"{'='*60}\n")
    
    for idx, filename in enumerate(all_files, 1):
        base = os.path.splitext(filename)[0]
        print(f"\n{'#'*60}")
        print(f"# Processing sample {idx}/{total_samples}: {base}")
        print(f"{'#'*60}")
        
        jobs = build_jobs_for_sample(filename)
        print(f"Created {len(jobs)} jobs for {base}")
        
        with ThreadPoolExecutor(max_workers=ncpus) as executor:
            results = list(executor.map(run_job, jobs))
        
        failed = sum(1 for r, _ in results if r != 0)
        if failed:
            print(f"\n⚠ {failed}/{len(jobs)} jobs failed for {base}")
        else:
            print(f"\n✓ All {len(jobs)} jobs succeeded for {base}")
        
        # Merge immediately after processing
        merge_sample(base)
    
    print(f"\n{'='*60}")
    print("✓ All samples processed!")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
