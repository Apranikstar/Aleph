import os
import uproot
import awkward as ak
import numpy as np
import pickle
from pathlib import Path

# ======================================================================
#  CONFIGURATION
# ======================================================================

# Directory containing ROOT files
ROOT_DIR = "/eos/experiment/aleph/EDM4HEP/MC/1994/QQB"
OUTPUT_DIR = "/eos/user/h/hfatehi/dEdxPickles/"

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define particle species
particles = {
    'protons': 2212,
    'pions': 211,
    'kaons': 321,
    'electrons': 11,
    'muons': 13
}

# Initialize storage for all particles
particle_data = {name: {
    'pads': {'type': [], 'value': [], 'error': [], 'momentum': []},
    'wires': {'type': [], 'value': [], 'error': [], 'momentum': []}
} for name in particles.keys()}

# ======================================================================
#  PROCESS ROOT FILES
# ======================================================================

# Get all ROOT files in directory
root_files = list(Path(ROOT_DIR).glob("*.root"))
print(f"Found {len(root_files)} ROOT files in {ROOT_DIR}")

for file_idx, root_file in enumerate(root_files):
    print(f"\nProcessing file {file_idx + 1}/{len(root_files)}: {root_file.name}")
    
    try:
        file = uproot.open(root_file)
        tree = file["events"]
        
        # Read all necessary data
        arrays = tree.arrays([
            "MCParticles/MCParticles.PDG",
            "MCParticles/MCParticles.momentum.x",
            "MCParticles/MCParticles.momentum.y",
            "MCParticles/MCParticles.momentum.z",
            "dEdxPads/dEdxPads.dQdx.type",
            "dEdxPads/dEdxPads.dQdx.value",
            "dEdxPads/dEdxPads.dQdx.error",
            "dEdxWires/dEdxWires.dQdx.type",
            "dEdxWires/dEdxWires.dQdx.value",
            "dEdxWires/dEdxWires.dQdx.error",
            "_dEdxPads_track/_dEdxPads_track.index",
            "_dEdxWires_track/_dEdxWires_track.index",
            "trackMCLink/trackMCLink.weight",
            "_trackMCLink_from/_trackMCLink_from.index",
            "_trackMCLink_to/_trackMCLink_to.index",
        ], library="ak")
        
        total_events = len(arrays)
        print(f"  Total events in file: {total_events}")
        print(f"  Processing all events...")
        
        # Process ALL events
        for evt in range(total_events):
            if (evt + 1) % 1000 == 0:
                print(f"    Processed {evt + 1}/{total_events} events...")
            
            pdg = arrays["MCParticles/MCParticles.PDG"][evt]
            px = arrays["MCParticles/MCParticles.momentum.x"][evt]
            py = arrays["MCParticles/MCParticles.momentum.y"][evt]
            pz = arrays["MCParticles/MCParticles.momentum.z"][evt]
            
            # Calculate momentum magnitude
            p_mag = np.sqrt(px**2 + py**2 + pz**2)
            
            # dEdx data
            dedx_pads_type = arrays["dEdxPads/dEdxPads.dQdx.type"][evt]
            dedx_pads_val = arrays["dEdxPads/dEdxPads.dQdx.value"][evt]
            dedx_pads_err = arrays["dEdxPads/dEdxPads.dQdx.error"][evt]
            dedx_pads_track = arrays["_dEdxPads_track/_dEdxPads_track.index"][evt]
            
            dedx_wires_type = arrays["dEdxWires/dEdxWires.dQdx.type"][evt]
            dedx_wires_val = arrays["dEdxWires/dEdxWires.dQdx.value"][evt]
            dedx_wires_err = arrays["dEdxWires/dEdxWires.dQdx.error"][evt]
            dedx_wires_track = arrays["_dEdxWires_track/_dEdxWires_track.index"][evt]
            
            # Track-MC links
            track_idx = arrays["_trackMCLink_from/_trackMCLink_from.index"][evt]
            mc_idx = arrays["_trackMCLink_to/_trackMCLink_to.index"][evt]
            weight = arrays["trackMCLink/trackMCLink.weight"][evt]
            
            # Create track to MC mapping
            track_to_mc = {}
            for t_idx, m_idx, w in zip(track_idx, mc_idx, weight):
                if t_idx not in track_to_mc or w > track_to_mc[t_idx][1]:
                    track_to_mc[t_idx] = (m_idx, w)
            
            # Process Pads
            for track, d_type, d_val, d_err in zip(dedx_pads_track, dedx_pads_type, dedx_pads_val, dedx_pads_err):
                if track in track_to_mc:
                    mc_particle_idx, _ = track_to_mc[track]
                    particle_pdg = abs(pdg[mc_particle_idx])
                    
                    # Check which particle type this is
                    for name, pdg_code in particles.items():
                        if particle_pdg == pdg_code:
                            particle_data[name]['pads']['type'].append(d_type)
                            particle_data[name]['pads']['value'].append(d_val)
                            particle_data[name]['pads']['error'].append(d_err)
                            particle_data[name]['pads']['momentum'].append(p_mag[mc_particle_idx])
                            break
            
            # Process Wires
            for track, d_type, d_val, d_err in zip(dedx_wires_track, dedx_wires_type, dedx_wires_val, dedx_wires_err):
                if track in track_to_mc:
                    mc_particle_idx, _ = track_to_mc[track]
                    particle_pdg = abs(pdg[mc_particle_idx])
                    
                    # Check which particle type this is
                    for name, pdg_code in particles.items():
                        if particle_pdg == pdg_code:
                            particle_data[name]['wires']['type'].append(d_type)
                            particle_data[name]['wires']['value'].append(d_val)
                            particle_data[name]['wires']['error'].append(d_err)
                            particle_data[name]['wires']['momentum'].append(p_mag[mc_particle_idx])
                            break
        
        print(f"  ✓ Completed all {total_events} events")
        file.close()
        
    except Exception as e:
        print(f"  ✗ Error processing {root_file.name}: {e}")
        continue

# ======================================================================
#  SAVE PICKLE FILES
# ======================================================================

print("\n" + "="*70)
print("SAVING PICKLE FILES")
print("="*70)

for name in particles.keys():
    # Convert lists to numpy arrays
    for detector in ['pads', 'wires']:
        for key in ['type', 'value', 'error', 'momentum']:
            particle_data[name][detector][key] = np.array(particle_data[name][detector][key])
    
    # Save to pickle file
    output_file = os.path.join(OUTPUT_DIR, f"{name}.pkl")
    with open(output_file, 'wb') as f:
        pickle.dump(particle_data[name], f)
    
    print(f"\n{name.capitalize()}:")
    print(f"  Pads measurements: {len(particle_data[name]['pads']['value'])}")
    print(f"  Wires measurements: {len(particle_data[name]['wires']['value'])}")
    print(f"  Saved to: {output_file}")

print("\n" + "="*70)
print("DONE!")
print("="*70)
