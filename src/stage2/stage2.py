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
    """Infer jet flavor from filename based on flavor mapping."""
    for pattern, label in flavor_mapping.items():
        if pattern in filename:
            return pattern, label
    
    patterns = ", ".join(flavor_mapping.keys())
    raise ValueError(
        f"Could not infer jet flavor from filename: {filename}\n"
        f"Expected one of: {patterns}"
    )


def sanitize_array(arr):
    """Convert NaN and None values to 0.0 in awkward arrays."""
    arr = ak.fill_none(arr, 0.0)
    arr = ak.where(np.isnan(arr), 0.0, arr)
    return arr


def process_chunk(
    input_file: str,
    output_file: str,
    config: dict,
    entry_start: int,
    entry_stop: int,
    flavor_label: str = None
) -> None:
    """Process a chunk of events from input file and write jet-level output."""
    
    if flavor_label is None:
        flavor_pattern, flavor_label = infer_jet_flavor(
            Path(input_file).name, 
            config['flavor_mapping']
        )
    else:
        valid_labels = set(config['flavor_mapping'].values())
        if flavor_label not in valid_labels:
            raise ValueError(f"Invalid flavor label: {flavor_label}")
        flavor_pattern = flavor_label 

    with uproot.open(input_file) as infile:
        tree = infile[config['input_tree_name']]
        n_entries = tree.num_entries
        
        if entry_stop > n_entries:
            entry_stop = n_entries
        
        jet_vars = config['jet_variables']
        pfcand_vars = config['pfcand_variables']
        all_vars = jet_vars + pfcand_vars
        
        arrays = tree.arrays(
            all_vars,
            entry_start=entry_start,
            entry_stop=entry_stop,
            library="ak"
        )

    output_data = {}
    
    # Flatten the event dimension to get a flat list of jets
    # This turns (events, jets, ...) into (total_jets, ...)
    for var in jet_vars:
        output_data[var] = ak.to_numpy(ak.flatten(arrays[var], axis=1))
    
    total_jets = len(output_data[jet_vars[0]])

    # Create flavor labels
    for pattern, label in config['flavor_mapping'].items():
        is_matched = int(pattern == flavor_pattern)
        output_data[label] = np.full(total_jets, is_matched, dtype=np.int32)

    # Process pfcand-level variables
    for var in pfcand_vars:
        # Step 1: Flatten axis 1 (the jet dimension)
        # Input shape: [event][jet][particle]
        # Output shape: [total_jets][particle] -> This is "var * float32"
        flattened = ak.flatten(arrays[var], axis=1)
        
        # Step 2: Sanitize
        flattened = sanitize_array(flattened)
        
        # Step 3: Ensure float32 for ROOT compatibility
        output_data[var] = ak.values_astype(flattened, np.float32)

    # Add jet_npfcand (count particles in each jet)
    # We use axis=1 on the already event-flattened array
    output_data['jet_npfcand'] = ak.to_numpy(ak.num(output_data[pfcand_vars[0]], axis=1)).astype(np.int32)

    # Write output file
    with uproot.recreate(output_file) as outfile:
        outfile[config['output_tree_name']] = output_data


def main():
    parser = argparse.ArgumentParser(description="Process jet data")
    parser.add_argument("input_file", type=str)
    parser.add_argument("output_file", type=str)
    parser.add_argument("entry_start", type=int)
    parser.add_argument("entry_stop", type=int)
    parser.add_argument("--config", type=str, default="config.yaml")
    parser.add_argument("--flavor-label", type=str, default=None)
    
    args = parser.parse_args()
    
    try:
        config = load_config(args.config)
        process_chunk(
            args.input_file,
            args.output_file,
            config,
            args.entry_start,
            args.entry_stop,
            args.flavor_label
        )
        print(f"Successfully wrote {args.output_file}")
    except Exception as e:
        print(f"ERROR: Processing failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
