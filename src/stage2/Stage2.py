#!/usr/bin/env python3
"""
Modern stage2 processing script using uproot.
Converts event-level jet data to jet-level output with flavor labels.
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import uproot
import awkward as ak
import yaml


def load_config(config_path: str = "config.yaml") -> dict:
    """Load YAML configuration file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def infer_jet_flavor(filename: str, flavor_mapping: dict) -> tuple[str, str]:
    """
    Infer jet flavor from filename based on flavor mapping.
    
    Args:
        filename: Input ROOT file name
        flavor_mapping: Dictionary mapping file patterns to flavor labels
    
    Returns:
        Tuple of (matched_pattern, flavor_label)
    
    Raises:
        ValueError: If no flavor pattern matches the filename
    """
    for pattern, label in flavor_mapping.items():
        if pattern in filename:
            return pattern, label
    
    patterns = ", ".join(flavor_mapping.keys())
    raise ValueError(
        f"Could not infer jet flavor from filename: {filename}\n"
        f"Expected one of: {patterns}"
    )


def sanitize_array(arr):
    """
    Convert NaN and None values to 0.0 in awkward arrays.
    
    Args:
        arr: Awkward array that may contain NaN or None values
    
    Returns:
        Sanitized array with NaN/None replaced by 0.0
    """
    # Handle None/null values
    arr = ak.fill_none(arr, 0.0)
    
    # Handle NaN values
    arr = ak.where(np.isnan(arr), 0.0, arr)
    
    return arr


def process_chunk(
    input_file: str,
    output_file: str,
    config: dict,
    entry_start: int,
    entry_stop: int
) -> None:
    """
    Process a chunk of events from input file and write jet-level output.
    
    Args:
        input_file: Path to input ROOT file
        output_file: Path to output ROOT file
        config: Configuration dictionary
        entry_start: First entry to process
        entry_stop: Last entry to process (exclusive)
    """
    # Infer jet flavor
    flavor_pattern, flavor_label = infer_jet_flavor(
        Path(input_file).name, 
        config['flavor_mapping']
    )
    
    # Open input file and get tree
    with uproot.open(input_file) as infile:
        tree = infile[config['input_tree_name']]
        
        # Get number of entries
        n_entries = tree.num_entries
        
        # Validate entry range
        if entry_stop > n_entries:
            raise ValueError(
                f"Requested entry range [{entry_start}, {entry_stop}) exceeds "
                f"available entries ({n_entries})"
            )
        
        # Get all variable names we need to read
        jet_vars = config['jet_variables']
        pfcand_vars = config['pfcand_variables']
        all_vars = jet_vars + pfcand_vars
        
        # Read the data for this chunk
        arrays = tree.arrays(
            all_vars,
            entry_start=entry_start,
            entry_stop=entry_stop,
            library="ak"
        )
    
    # Prepare output data - we'll flatten jets across all events
    output_data = {}
    
    # Get number of jets per event
    n_jets_per_event = ak.num(arrays[jet_vars[0]], axis=0)
    total_jets = ak.sum(n_jets_per_event)
    
    # Create flavor label branches (one for each flavor, only matched one is 1)
    for pattern, label in config['flavor_mapping'].items():
        is_matched = int(pattern == flavor_pattern)
        output_data[label] = np.full(total_jets, is_matched, dtype=np.int32)
    
    # Flatten jet-level variables
    for var in jet_vars:
        output_data[var] = ak.to_numpy(ak.flatten(arrays[var]))
    
    # Process pfcand-level variables (vector of vectors -> jagged array per jet)
    # First, we need to flatten the outer dimension (events) but keep inner (constituents)
    for var in pfcand_vars:
        # arrays[var] has shape: [[event1_jet1_constituents, event1_jet2_constituents], [event2_jet1_constituents], ...]
        # We want: [event1_jet1_constituents, event1_jet2_constituents, event2_jet1_constituents, ...]
        flattened = ak.flatten(arrays[var], axis=0)
        
        # Sanitize NaN and None values
        flattened = sanitize_array(flattened)
        
        output_data[var] = flattened
    
    # Add jet_npfcand - number of constituents per jet
    output_data['jet_npfcand'] = ak.to_numpy(ak.num(output_data[pfcand_vars[0]], axis=1))
    
    # Write output file
    with uproot.recreate(output_file) as outfile:
        outfile[config['output_tree_name']] = output_data


def main():
    parser = argparse.ArgumentParser(
        description="Process jet data from event-level to jet-level format"
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Input ROOT file path"
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="Output ROOT file path"
    )
    parser.add_argument(
        "entry_start",
        type=int,
        help="First entry to process"
    )
    parser.add_argument(
        "entry_stop",
        type=int,
        help="Last entry to process (exclusive)"
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config.yaml",
        help="Path to configuration YAML file (default: config.yaml)"
    )
    
    args = parser.parse_args()
    
    # Load configuration
    try:
        config = load_config(args.config)
    except FileNotFoundError:
        print(f"ERROR: Configuration file not found: {args.config}")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"ERROR: Failed to parse configuration file: {e}")
        sys.exit(1)
    
    # Process the chunk
    try:
        process_chunk(
            args.input_file,
            args.output_file,
            config,
            args.entry_start,
            args.entry_stop
        )
    except Exception as e:
        print(f"ERROR: Processing failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
