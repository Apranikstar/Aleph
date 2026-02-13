#!/usr/bin/env python3
"""
Steering file for parallel jet flavor tagging processing.
Configure your processing parameters here and run.
"""

from pathlib import Path
import sys
import subprocess


# ============================================================================
# USER CONFIGURATION - EDIT THESE PARAMETERS
# ============================================================================

# Input/Output directories
INPUT_DIR = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/1994/zqq/stage1/temp"
OUTPUT_DIR = "/eos/user/h/hfatehi/aleph-temp2/"

# Processing parameters
N_WORKERS = 64                      # Number of parallel workers
CHUNK_SIZE = None                   # Entries per chunk (None = auto-calculate)
N_FINAL_FILES = 1                   # Number of final merged files per sample

# File management
DELETE_PARTS_AFTER_MERGE = True     # Delete intermediate part files after merging
KEEP_LOGS = True                    # Keep stdout/stderr logs

# Script paths (usually don't need to change these)
STAGE2_SCRIPT = "stage2.py"
CONFIG_FILE = "config.yaml"

# Optional: Filter input files (None = process all .root files)
# Examples:
#   FILE_PATTERN = "Zbb*.root"      # Only Zbb files
#   FILE_PATTERN = "Zcc*.root"      # Only Zcc files
#   FILE_PATTERN = None             # All .root files
FILE_PATTERN = None

# Optional: Dry run (print commands without executing)
DRY_RUN = False

# ============================================================================
# END OF USER CONFIGURATION
# ============================================================================


def run_processing():
    """Execute the parallel processing pipeline."""
    
    # Build command
    cmd = [
        "python",
        "run_parallel.py",
        INPUT_DIR,
        OUTPUT_DIR,
        "--workers", str(N_WORKERS),
        "--config", CONFIG_FILE,
        "--stage2-script", STAGE2_SCRIPT,
        "--n-final-files", str(N_FINAL_FILES),
    ]
    
    # Add optional parameters
    if CHUNK_SIZE is not None:
        cmd.extend(["--chunk-size", str(CHUNK_SIZE)])
    
    if not DELETE_PARTS_AFTER_MERGE:
        cmd.append("--keep-parts")
    
    # Print configuration
    print("=" * 80)
    print("JET FLAVOR TAGGING - PARALLEL PROCESSING STEERING")
    print("=" * 80)
    print()
    print("Configuration:")
    print(f"  Input directory:       {INPUT_DIR}")
    print(f"  Output directory:      {OUTPUT_DIR}")
    print(f"  Number of workers:     {N_WORKERS}")
    print(f"  Chunk size:            {CHUNK_SIZE if CHUNK_SIZE else 'auto'}")
    print(f"  Final files per sample: {N_FINAL_FILES}")
    print(f"  Delete intermediate:   {DELETE_PARTS_AFTER_MERGE}")
    print(f"  File pattern:          {FILE_PATTERN if FILE_PATTERN else 'all .root files'}")
    print(f"  Config file:           {CONFIG_FILE}")
    print(f"  Stage2 script:         {STAGE2_SCRIPT}")
    print()
    print("Command:")
    print(f"  {' '.join(cmd)}")
    print()
    print("=" * 80)
    print()
    
    # Validate paths
    input_path = Path(INPUT_DIR)
    if not input_path.exists():
        print(f"ERROR: Input directory does not exist: {INPUT_DIR}")
        sys.exit(1)
    
    config_path = Path(CONFIG_FILE)
    if not config_path.exists():
        print(f"ERROR: Config file does not exist: {CONFIG_FILE}")
        sys.exit(1)
    
    stage2_path = Path(STAGE2_SCRIPT)
    if not stage2_path.exists():
        print(f"ERROR: Stage2 script does not exist: {STAGE2_SCRIPT}")
        sys.exit(1)
    
    # Count input files
    if FILE_PATTERN:
        input_files = list(input_path.glob(FILE_PATTERN))
    else:
        input_files = list(input_path.glob("*.root"))
    
    if not input_files:
        print(f"WARNING: No input files found matching pattern: {FILE_PATTERN or '*.root'}")
        response = input("Continue anyway? [y/N]: ")
        if response.lower() != 'y':
            sys.exit(0)
    else:
        print(f"Found {len(input_files)} input file(s):")
        for f in sorted(input_files)[:10]:  # Show first 10
            print(f"  - {f.name}")
        if len(input_files) > 10:
            print(f"  ... and {len(input_files) - 10} more")
        print()
    
    # Dry run mode
    if DRY_RUN:
        print("DRY RUN MODE - Command would be executed:")
        print(f"  {' '.join(cmd)}")
        print()
        print("Set DRY_RUN = False in this steering file to actually run the processing.")
        return
    
    # Confirm before running
    print("Ready to start processing.")
    response = input("Proceed? [Y/n]: ")
    if response.lower() == 'n':
        print("Aborted by user.")
        sys.exit(0)
    
    print()
    print("=" * 80)
    print("Starting processing...")
    print("=" * 80)
    print()
    
    # Run the command
    try:
        if KEEP_LOGS:
            # Run with real-time output
            result = subprocess.run(cmd)
            exit_code = result.returncode
        else:
            # Run with captured output
            result = subprocess.run(cmd, capture_output=True, text=True)
            exit_code = result.returncode
            if exit_code != 0:
                print("STDERR:")
                print(result.stderr)
    
    except KeyboardInterrupt:
        print("\n\nProcessing interrupted by user (Ctrl+C)")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nERROR: Processing failed with exception: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Report final status
    print()
    print("=" * 80)
    if exit_code == 0:
        print("✓ PROCESSING COMPLETED SUCCESSFULLY")
        print("=" * 80)
        print()
        print(f"Output files are in: {OUTPUT_DIR}")
    else:
        print("✗ PROCESSING FAILED")
        print("=" * 80)
        print()
        print(f"Exit code: {exit_code}")
        print("Check the error messages above for details.")
    
    sys.exit(exit_code)


def main():
    """Main entry point."""
    
    # Simple command-line interface
    if len(sys.argv) > 1:
        if sys.argv[1] in ["-h", "--help"]:
            print(__doc__)
            print("\nUsage:")
            print("  python steering.py              # Run with parameters set in this file")
            print("  python steering.py --help       # Show this help")
            print()
            print("Edit the configuration parameters at the top of this file.")
            sys.exit(0)
    
    run_processing()


if __name__ == "__main__":
    main()
