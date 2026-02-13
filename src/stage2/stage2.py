#!/usr/bin/env python3
"""
Stage2 processing using PyROOT (simpler for jagged arrays).
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import uproot
import awkward as ak
import yaml
from array import array
from ROOT import TFile, TTree


def load_config(config_path: str = "config.yaml") -> dict:
    """Load YAML configuration file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def infer_jet_flavor(filename: str, flavor_mapping: dict) -> tuple[str, str]:
    """Infer jet flavor from filename."""
    for pattern, label in flavor_mapping.items():
        if pattern in filename:
            return pattern, label
    
    patterns = ", ".join(flavor_mapping.keys())
    raise ValueError(
        f"Could not infer jet flavor from filename: {filename}\n"
        f"Expected one of: {patterns}"
    )


def sanitize_value(val):
    """Convert NaN/None to 0.0."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return 0.0
    return float(val)


def process_chunk(
    input_file: str,
    output_file: str,
    config: dict,
    entry_start: int,
    entry_stop: int,
    flavor_label: str = None
) -> None:
    """Process a chunk using PyROOT."""
    
    # Infer jet flavor
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
    
    # Get variable lists
    jet_vars = config['jet_variables']
    pfcand_vars = config['pfcand_variables']
    all_vars = jet_vars + pfcand_vars
    
    # Read data with uproot
    with uproot.open(input_file) as infile:
        tree = infile[config['input_tree_name']]
        
        if entry_stop > tree.num_entries:
            raise ValueError(
                f"Requested entry range [{entry_start}, {entry_stop}) exceeds "
                f"available entries ({tree.num_entries})"
            )
        
        arrays = tree.arrays(
            all_vars,
            entry_start=entry_start,
            entry_stop=entry_stop,
            library="ak"
        )
    
    # Create output file with PyROOT
    outfile = TFile(output_file, "RECREATE")
    outtree = TTree(config['output_tree_name'], config['output_tree_name'])
    
    # Create branches for flavor labels
    flavor_branches = {}
    for pattern, label in config['flavor_mapping'].items():
        is_matched = int(pattern == flavor_pattern)
        flavor_branches[label] = array('i', [is_matched])
        outtree.Branch(label, flavor_branches[label], f"{label}/I")
    
    # Create branches for jet variables
    jet_branches = {}
    for var in jet_vars:
        jet_branches[var] = array('f', [0.0])
        outtree.Branch(var, jet_branches[var], f"{var}/F")
    
    # Create branches for pfcand variables (jagged)
    max_pfcand = 500
    jet_npfcand = array('i', [0])
    outtree.Branch("jet_npfcand", jet_npfcand, "jet_npfcand/I")
    
    pfcand_branches = {}
    for var in pfcand_vars:
        pfcand_branches[var] = array('f', max_pfcand * [0.0])
        outtree.Branch(var, pfcand_branches[var], f"{var}[jet_npfcand]/F")
    
    # Process events
    n_events = entry_stop - entry_start
    for i in range(n_events):
        # Get number of jets in this event
        n_jets = len(arrays[jet_vars[0]][i])
        
        # Loop over jets
        for j in range(n_jets):
            # Fill flavor labels
            for label, branch in flavor_branches.items():
                # Already set, no need to change
                pass
            
            # Fill jet variables
            for var in jet_vars:
                jet_branches[var][0] = float(arrays[var][i][j])
            
            # Fill pfcand variables
            n_const = len(arrays[pfcand_vars[0]][i][j])
            jet_npfcand[0] = min(n_const, max_pfcand)
            
            for var in pfcand_vars:
                constituents = arrays[var][i][j]
                for k in range(min(n_const, max_pfcand)):
                    pfcand_branches[var][k] = sanitize_value(constituents[k])
                # Fill remaining with zeros
                for k in range(n_const, max_pfcand):
                    pfcand_branches[var][k] = 0.0
            
            # Fill tree for this jet
            outtree.Fill()
    
    # Write and close
    outtree.Write()
    outfile.Close()


def main():
    parser = argparse.ArgumentParser(
        description="Process jet data from event-level to jet-level format"
    )
    parser.add_argument("input_file", type=str, help="Input ROOT file path")
    parser.add_argument("output_file", type=str, help="Output ROOT file path")
    parser.add_argument("entry_start", type=int, help="First entry to process")
    parser.add_argument("entry_stop", type=int, help="Last entry to process (exclusive)")
    parser.add_argument("--config", type=str, default="config.yaml", help="Config file")
    parser.add_argument("--flavor-label", type=str, default=None, help="Override flavor detection")
    
    args = parser.parse_args()
    
    try:
        config = load_config(args.config)
    except FileNotFoundError:
        print(f"ERROR: Configuration file not found: {args.config}")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"ERROR: Failed to parse configuration file: {e}")
        sys.exit(1)
    
    try:
        process_chunk(
            args.input_file,
            args.output_file,
            config,
            args.entry_start,
            args.entry_stop,
            args.flavor_label
        )
    except Exception as e:
        print(f"ERROR: Processing failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
