#!/usr/bin/env python3
#python check_labels.py --input '/eos/user/h/hfatehi/aleph/trainingFiles/aleph-temp2/Z*.root' --sample
import uproot
import awkward as ak
import numpy as np
from glob import glob
from tqdm import tqdm
import argparse

def check_labels_in_file(filename):
    """Check label consistency in a single ROOT file"""
    try:
        with uproot.open(filename) as f:
            # Get the tree
            tree = f.get("tree") if "tree" in f else f[list(f.keys())[0]]
            
            # Read the label branches - now checking for individual U,D,S labels
            branches = tree.arrays(["recojet_isB", "recojet_isC", "recojet_isU", "recojet_isD", "recojet_isS", "jet_pT"], 
                                  library="ak")
            
            reco_jet_isB = branches["recojet_isB"]
            reco_jet_isC = branches["recojet_isC"]
            reco_jet_isU = branches["recojet_isU"]
            reco_jet_isD = branches["recojet_isD"]
            reco_jet_isS = branches["recojet_isS"]
            jet_pt = branches["jet_pT"]
            
            # Create combined UDS label
            reco_jet_isUDS = (reco_jet_isU + reco_jet_isD + reco_jet_isS) > 0
            
            # Sum of labels (should be 1 for each jet)
            label_sum = reco_jet_isB + reco_jet_isC + reco_jet_isUDS
            
            # Apply pT cut
            pt_cut = jet_pt > 10
            
            # Jets with no label (sum == 0)
            no_label = label_sum == 0
            no_label_after_pt = (label_sum == 0) & pt_cut
            
            # Jets with multiple labels (sum > 1)
            multi_label = label_sum > 1
            multi_label_after_pt = (label_sum > 1) & pt_cut
            
            # Label distribution
            b_jets = reco_jet_isB == 1
            c_jets = reco_jet_isC == 1
            uds_jets = reco_jet_isUDS == 1
            
            results = {
                'file': filename,
                'total_jets': len(jet_pt),
                'total_jets_pt10': ak.sum(pt_cut),
                'no_label': ak.sum(no_label),
                'no_label_pt10': ak.sum(no_label_after_pt),
                'multi_label': ak.sum(multi_label),
                'multi_label_pt10': ak.sum(multi_label_after_pt),
                'b_jets': ak.sum(b_jets),
                'b_jets_pt10': ak.sum(b_jets & pt_cut),
                'c_jets': ak.sum(c_jets),
                'c_jets_pt10': ak.sum(c_jets & pt_cut),
                'uds_jets': ak.sum(uds_jets),
                'uds_jets_pt10': ak.sum(uds_jets & pt_cut),
                # Also track individual U,D,S for debugging
                'u_jets': ak.sum(reco_jet_isU),
                'd_jets': ak.sum(reco_jet_isD),
                's_jets': ak.sum(reco_jet_isS),
            }
            
            return results
            
    except Exception as e:
        print(f"Error processing {filename}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description='Check jet labels in ROOT files')
    parser.add_argument('--input', '-i', type=str, required=True,
                       help='Input file pattern (e.g., "/eos/user/.../Z*.root")')
    parser.add_argument('--sample', '-s', action='store_true',
                       help='Group results by sample type (Zbb, Zcc, etc.)')
    args = parser.parse_args()
    
    # Get list of files
    files = glob(args.input)
    print(f"Found {len(files)} files")
    
    if args.sample:
        # Group by sample type
        samples = {}
        for f in files:
            if 'Zbb' in f:
                samples.setdefault('Zbb', []).append(f)
            elif 'Zcc' in f:
                samples.setdefault('Zcc', []).append(f)
            elif 'Zdd' in f:
                samples.setdefault('Zdd', []).append(f)
            elif 'Zss' in f:
                samples.setdefault('Zss', []).append(f)
            elif 'Zuu' in f:
                samples.setdefault('Zuu', []).append(f)
            else:
                samples.setdefault('Other', []).append(f)
        
        # Analyze each sample type
        all_results = {}
        for sample_name, sample_files in samples.items():
            print(f"\n{'='*60}")
            print(f"Analyzing {sample_name} ({len(sample_files)} files)")
            print(f"{'='*60}")
            
            sample_results = []
            for f in tqdm(sample_files, desc=f"Processing {sample_name}"):
                result = check_labels_in_file(f)
                if result:
                    sample_results.append(result)
            
            if sample_results:
                # Aggregate results
                total = {
                    'total_jets': sum(r['total_jets'] for r in sample_results),
                    'total_jets_pt10': sum(r['total_jets_pt10'] for r in sample_results),
                    'no_label': sum(r['no_label'] for r in sample_results),
                    'no_label_pt10': sum(r['no_label_pt10'] for r in sample_results),
                    'multi_label': sum(r['multi_label'] for r in sample_results),
                    'multi_label_pt10': sum(r['multi_label_pt10'] for r in sample_results),
                    'b_jets': sum(r['b_jets'] for r in sample_results),
                    'b_jets_pt10': sum(r['b_jets_pt10'] for r in sample_results),
                    'c_jets': sum(r['c_jets'] for r in sample_results),
                    'c_jets_pt10': sum(r['c_jets_pt10'] for r in sample_results),
                    'uds_jets': sum(r['uds_jets'] for r in sample_results),
                    'uds_jets_pt10': sum(r['uds_jets_pt10'] for r in sample_results),
                    'u_jets': sum(r['u_jets'] for r in sample_results),
                    'd_jets': sum(r['d_jets'] for r in sample_results),
                    's_jets': sum(r['s_jets'] for r in sample_results),
                }
                all_results[sample_name] = total
                
                # Print summary
                print(f"\nResults for {sample_name}:")
                print(f"  Total jets: {total['total_jets']}")
                print(f"  Jets with pT > 10: {total['total_jets_pt10']}")
                print(f"  Jets with NO label: {total['no_label']} ({total['no_label']/total['total_jets']*100:.2f}%)")
                print(f"  Jets with NO label (pT>10): {total['no_label_pt10']} ({total['no_label_pt10']/total['total_jets_pt10']*100:.2f}%)")
                print(f"  Jets with MULTIPLE labels: {total['multi_label']} ({total['multi_label']/total['total_jets']*100:.2f}%)")
                print(f"  B jets: {total['b_jets']} (pT>10: {total['b_jets_pt10']})")
                print(f"  C jets: {total['c_jets']} (pT>10: {total['c_jets_pt10']})")
                print(f"  UDS jets: {total['uds_jets']} (pT>10: {total['uds_jets_pt10']})")
                print(f"    - U jets: {total['u_jets']}")
                print(f"    - D jets: {total['d_jets']}")
                print(f"    - S jets: {total['s_jets']}")
        
        # Overall summary
        if all_results:
            print(f"\n{'='*60}")
            print("OVERALL SUMMARY")
            print(f"{'='*60}")
            
            overall = {
                'total_jets': sum(r['total_jets'] for r in all_results.values()),
                'total_jets_pt10': sum(r['total_jets_pt10'] for r in all_results.values()),
                'no_label': sum(r['no_label'] for r in all_results.values()),
                'no_label_pt10': sum(r['no_label_pt10'] for r in all_results.values()),
                'multi_label': sum(r['multi_label'] for r in all_results.values()),
                'multi_label_pt10': sum(r['multi_label_pt10'] for r in all_results.values()),
            }
            
            print(f"Total jets across all samples: {overall['total_jets']}")
            print(f"Jets with pT > 10: {overall['total_jets_pt10']}")
            print(f"Total jets with NO label: {overall['no_label']} ({overall['no_label']/overall['total_jets']*100:.2f}%)")
            print(f"Jets with NO label (pT>10): {overall['no_label_pt10']} ({overall['no_label_pt10']/overall['total_jets_pt10']*100:.2f}%)")
            print(f"Total jets with MULTIPLE labels: {overall['multi_label']} ({overall['multi_label']/overall['total_jets']*100:.2f}%)")

if __name__ == "__main__":
    main()
