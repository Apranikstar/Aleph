#!/usr/bin/env python3
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
import uproot
import glob

# --- user configuration ---
input_dir = "/afs/cern.ch/work/h/hfatehi/Aleph/prod/database/s1"
output_dir = "/afs/cern.ch/work/h/hfatehi/Aleph/prod/database/s2"
stage2_script = "stage2.py"
ncpus = 4  # number of parallel threads
delete_parts_after_merge = False  # set to True to remove part files after merging
# ---------------------------

os.makedirs(output_dir, exist_ok=True)

def get_nentries(root_file):
    """Return number of entries in the first TTree of the ROOT file."""
    with uproot.open(root_file) as f:
        # try to find a TTree key
        tree_name = next((k for k in f.keys() if hasattr(f[k], "num_entries")), list(f.keys())[0])
        tree = f[tree_name]
        return tree.num_entries

def build_jobs():
    """Prepare all (command, label) jobs across input ROOT files."""
    jobs = []
    for filename in os.listdir(input_dir):
        if not filename.endswith(".root"):
            continue
        input_path = os.path.join(input_dir, filename)
        base = os.path.splitext(filename)[0]

        n_entries = get_nentries(input_path)
        chunk_size = (n_entries + ncpus - 1) // ncpus

        for i in range(ncpus):
            start = i * chunk_size
            end = min((i + 1) * chunk_size, n_entries)
            if start >= end:
                break
            output_file = os.path.join(output_dir, f"{base}_part{i}.root")
            cmd = ["python", stage2_script, input_path, output_file, str(start), str(end)]
            jobs.append((cmd, f"{base}_part{i}"))
    return jobs

def run_job(job):
    """Run one job command."""
    cmd, label = job
    print(f"[START] {label}: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] {label}:\n{result.stderr.strip()}")
    else:
        print(f"[DONE] {label}")
    return result.returncode

def merge_outputs():
    """Merge all _part*.root files into _s2.root per sample."""
    print("\n=== Merging output parts per sample ===")
    for filename in os.listdir(input_dir):
        if not filename.endswith(".root"):
            continue
        base = os.path.splitext(filename)[0]
        part_files = sorted(glob.glob(os.path.join(output_dir, f"{base}_part*.root")))
        if not part_files:
            continue

        merged_file = os.path.join(output_dir, f"{base}_s2.root")
        print(f"Merging {len(part_files)} parts → {merged_file}")
        cmd = ["hadd", "-f", merged_file] + part_files
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"[ERROR] Failed to merge {base}: {result.stderr.strip()}")
        else:
            print(f"[OK] {base} merged successfully.")
            if delete_parts_after_merge:
                for f in part_files:
                    try:
                        os.remove(f)
                    except Exception as e:
                        print(f"Warning: could not remove {f}: {e}")

def main():
    jobs = build_jobs()
    print(f"Prepared {len(jobs)} total jobs across all ROOT files.")
    with ThreadPoolExecutor(max_workers=ncpus) as executor:
        results = list(executor.map(run_job, jobs))

    print("\n=== All jobs completed ===")
    failed = sum(1 for r in results if r != 0)
    if failed:
        print(f"{failed} jobs failed.")
    else:
        print("All jobs succeeded.")

    merge_outputs()
    print("\nAll done!")

if __name__ == "__main__":
    main()
