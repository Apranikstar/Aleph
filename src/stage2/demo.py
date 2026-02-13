#!/usr/bin/env python3
"""
Example script showing how to use the modernized pipeline.
This demonstrates the key features and typical usage patterns.
"""

import yaml
import numpy as np
import uproot
import awkward as ak


def create_example_config():
    """Create an example config file."""
    config = {
        'flavor_mapping': {
            'Zbb': 'recojet_isB',
            'Zcc': 'recojet_isC',
            'Zuu': 'recojet_isU',
        },
        'jet_variables': [
            'jet_p',
            'jet_e',
            'jet_mass',
        ],
        'pfcand_variables': [
            'pfcand_erel_log',
            'pfcand_btagSip3dSig',
            'pfcand_dxy',
        ],
        'input_tree_name': 'events',
        'output_tree_name': 'tree',
    }
    
    with open('example_config.yaml', 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    
    print("Created example_config.yaml")


def create_dummy_input_file():
    """Create a dummy input file for testing."""
    # Create fake event-level data
    n_events = 10
    
    data = {}
    
    # Jet variables (vector per event)
    data['jet_p'] = ak.Array([
        np.random.uniform(20, 80, size=np.random.randint(1, 4))
        for _ in range(n_events)
    ])
    data['jet_e'] = ak.Array([
        np.random.uniform(20, 80, size=len(data['jet_p'][i]))
        for i in range(n_events)
    ])
    data['jet_mass'] = ak.Array([
        np.random.uniform(0, 10, size=len(data['jet_p'][i]))
        for i in range(n_events)
    ])
    
    # PFCand variables (vector of vectors per event)
    data['pfcand_erel_log'] = ak.Array([
        [
            np.random.uniform(-3, 0, size=np.random.randint(5, 20))
            for _ in range(len(data['jet_p'][i]))
        ]
        for i in range(n_events)
    ])
    
    data['pfcand_btagSip3dSig'] = ak.Array([
        [
            np.random.uniform(-5, 5, size=len(data['pfcand_erel_log'][i][j]))
            for j in range(len(data['pfcand_erel_log'][i]))
        ]
        for i in range(n_events)
    ])
    
    # Add some NaN values to test sanitization
    data['pfcand_dxy'] = ak.Array([
        [
            np.random.choice([np.random.uniform(-0.5, 0.5), np.nan], size=len(data['pfcand_erel_log'][i][j]), p=[0.9, 0.1])
            for j in range(len(data['pfcand_erel_log'][i]))
        ]
        for i in range(n_events)
    ])
    
    # Write to ROOT file
    with uproot.recreate('example_Zbb_input.root') as f:
        f['events'] = data
    
    print("Created example_Zbb_input.root")
    print(f"  - {n_events} events")
    print(f"  - {ak.sum(ak.num(data['jet_p']))} total jets")
    print(f"  - Some NaN values in pfcand_dxy for testing")


def inspect_output_file(filename):
    """Inspect the output file to verify structure."""
    print(f"\nInspecting {filename}:")
    print("-" * 60)
    
    with uproot.open(filename) as f:
        tree = f['tree']
        
        print(f"Branches: {tree.keys()}")
        print(f"\nNumber of entries (jets): {tree.num_entries}")
        
        # Read a few branches
        arrays = tree.arrays(['jet_p', 'jet_e', 'recojet_isB', 'recojet_isC', 'jet_npfcand'], library='ak')
        
        print(f"\nFirst 5 jets:")
        print(f"  jet_p: {arrays['jet_p'][:5]}")
        print(f"  jet_e: {arrays['jet_e'][:5]}")
        print(f"  recojet_isB: {arrays['recojet_isB'][:5]}")
        print(f"  recojet_isC: {arrays['recojet_isC'][:5]}")
        print(f"  jet_npfcand: {arrays['jet_npfcand'][:5]}")
        
        # Check pfcand variables
        pfcand = tree.arrays(['pfcand_erel_log', 'pfcand_dxy'], library='ak')
        print(f"\nFirst jet constituents:")
        print(f"  pfcand_erel_log[0]: {pfcand['pfcand_erel_log'][0]}")
        print(f"  pfcand_dxy[0]: {pfcand['pfcand_dxy'][0]}")
        
        # Verify no NaN values
        has_nan = ak.any(ak.is_none(pfcand['pfcand_dxy'], axis=1))
        print(f"\nNaN check:")
        print(f"  pfcand_dxy has None values: {has_nan}")
        # Flatten and check for actual NaN
        flat_dxy = ak.flatten(pfcand['pfcand_dxy'])
        print(f"  pfcand_dxy has NaN values: {np.any(np.isnan(ak.to_numpy(flat_dxy)))}")


def main():
    print("=" * 60)
    print("Example Usage of Modern Jet Flavor Tagging Pipeline")
    print("=" * 60)
    print()
    
    # Step 1: Create example config
    print("Step 1: Creating example configuration...")
    create_example_config()
    print()
    
    # Step 2: Create dummy input file
    print("Step 2: Creating example input file...")
    create_dummy_input_file()
    print()
    
    # Step 3: Show how to run processing
    print("Step 3: Running processing...")
    print("\nCommand to run:")
    print("  python stage2.py example_Zbb_input.root example_output.root 0 10 --config example_config.yaml")
    print()
    
    import subprocess
    result = subprocess.run([
        'python', 'stage2.py',
        'example_Zbb_input.root',
        'example_output.root',
        '0', '10',
        '--config', 'example_config.yaml'
    ], capture_output=True, text=True)
    
    if result.returncode == 0:
        print("✓ Processing completed successfully")
        print()
        
        # Step 4: Inspect output
        print("Step 4: Inspecting output file...")
        inspect_output_file('example_output.root')
        
    else:
        print("✗ Processing failed:")
        print(result.stderr)
    
    print()
    print("=" * 60)
    print("Example complete!")
    print("=" * 60)
    print()
    print("Next steps:")
    print("  1. Edit config.yaml for your actual variables")
    print("  2. Run run_parallel.py for large-scale processing")
    print("  3. Check README.md for full documentation")


if __name__ == "__main__":
    main()
