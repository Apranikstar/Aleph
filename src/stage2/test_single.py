#!/usr/bin/env python3
"""
Quick diagnostic script to test stage2.py on a single file chunk.
This will show the full error message.
"""

import subprocess
import sys
import os

# Your actual directories
INPUT_DIR = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/1994/zqq/stage1/temp/"
OUTPUT_DIR = "/eos/user/h/hfatehi/aleph-temp"

# Make sure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Test with one of your files
input_file = os.path.join(INPUT_DIR, "Zbb.root")
output_file = os.path.join(OUTPUT_DIR, "test_output.root")
start = 0
end = 100  # Just test 100 entries

cmd = [
    "python", "stage2.py",
    input_file,
    output_file,
    str(start),
    str(end),
    "--config", "config.yaml"
]

print("Testing stage2.py with a small sample...")
print("="*80)
print(f"Input:  {input_file}")
print(f"Output: {output_file}")
print(f"Entries: {start} to {end}")
print("="*80)
print()
print("Command:")
print(" ".join(cmd))
print()
print("="*80)
print()

result = subprocess.run(cmd, capture_output=True, text=True)

if result.stdout:
    print("STDOUT:")
    print(result.stdout)
    print()

if result.stderr:
    print("STDERR:")
    print(result.stderr)
    print()

print("="*80)
print(f"Return code: {result.returncode}")
print("="*80)

if result.returncode != 0:
    print("\n❌ FAILED - See error above")
    sys.exit(1)
else:
    print("\n✅ SUCCESS!")
    print(f"Output file created: {output_file}")
    
    # Check the output file
    try:
        import uproot
        with uproot.open(output_file) as f:
            tree = f['tree']
            print(f"Output tree has {tree.num_entries} entries")
            print(f"Branches: {list(tree.keys())[:10]}...")
    except Exception as e:
        print(f"Could not inspect output file: {e}")
