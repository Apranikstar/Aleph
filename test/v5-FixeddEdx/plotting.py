### 1-10k sample events picked


import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import os

# ======================================================================
#  MATPLOTLIB CONFIGURATION
# ======================================================================

def configure_matplotlib():
    os.environ['PATH'] = '/cvmfs/sft.cern.ch/lcg/external/texlive/2020/bin/x86_64-linux:' + os.environ['PATH']

    rcParams.update({
        "pgf.texsystem": "pdflatex",
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8,
        "figure.titlesize": 12,
        "lines.linewidth": 1.2,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6,
        "ytick.minor.width": 0.6,
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "xtick.minor.size": 2.5,
        "ytick.minor.size": 2.5,
        "legend.framealpha": 0.0,
        "legend.edgecolor": "none",
        "legend.fancybox": False,
    })


# ======================================================================
#  LOAD AND FLATTEN DATA
# ======================================================================

def load_and_flatten_data(filename, n_events=None, seed=42):
    """
    Load and flatten dE/dx data from ROOT file.
    
    Parameters:
    -----------
    filename : str
        Path to ROOT file
    n_events : int or None
        Number of events to randomly sample. If None, use all events. Default: None
    seed : int
        Random seed for reproducible subsampling. Default: 42
    """
    file = uproot.open(filename)
    tree = file["events"]
    
    # Set random seed for reproducibility
    np.random.seed(seed)
    
    # Determine which events to load
    total_events = tree.num_entries
    if n_events is not None and n_events < total_events:
        event_indices = np.sort(np.random.choice(total_events, size=n_events, replace=False))
        print(f"  Randomly selecting {n_events} out of {total_events} events")
    else:
        event_indices = np.arange(total_events)
        if n_events is not None:
            print(f"  Requested {n_events} events but file only has {total_events}, using all")

    # Load branches for selected events
    p_all = tree["pfcand_p"].array(library="np", entry_start=0, entry_stop=total_events)[event_indices]

    pads_values_all = {
        "ch": tree["pfcand_dEdx_ch_value_pads"].array(library="np", entry_start=0, entry_stop=total_events)[event_indices],
        "el": tree["pfcand_dEdx_el_value_pads"].array(library="np", entry_start=0, entry_stop=total_events)[event_indices],
        "mu": tree["pfcand_dEdx_mu_value_pads"].array(library="np", entry_start=0, entry_stop=total_events)[event_indices],
    }

    wires_values_all = {
        "ch": tree["pfcand_dEdx_ch_wires_value"].array(library="np", entry_start=0, entry_stop=total_events)[event_indices],
        "el": tree["pfcand_dEdx_el_wires_value"].array(library="np", entry_start=0, entry_stop=total_events)[event_indices],
        "mu": tree["pfcand_dEdx_mu_wires_value"].array(library="np", entry_start=0, entry_stop=total_events)[event_indices],
    }

    # Flatten nested RVecs and filter valid data
    data = {}
    for label in ["ch", "el", "mu"]:
        p_flat_pads, pads_values_flat = [], []
        p_flat_wires, wires_values_flat = [], []

        for ev_idx in range(len(p_all)):
            p_event = p_all[ev_idx]
            pads_values_event = pads_values_all[label][ev_idx]
            wires_values_event = wires_values_all[label][ev_idx]

            for jet_idx in range(len(pads_values_event)):
                pads_values_jet = pads_values_event[jet_idx]
                wires_values_jet = wires_values_event[jet_idx]
                p_jet = p_event[jet_idx]

                for const_idx in range(len(pads_values_jet)):
                    p_val = float(p_jet[const_idx])
                    
                    # Pads
                    pads_val = float(pads_values_jet[const_idx])
                    if pads_val > 0:
                        p_flat_pads.append(p_val)
                        pads_values_flat.append(pads_val)
                    
                    # Wires
                    wires_val = float(wires_values_jet[const_idx])
                    if wires_val > 0:
                        p_flat_wires.append(p_val)
                        wires_values_flat.append(wires_val)

        # Convert to numpy arrays
        p_pads_arr = np.array(p_flat_pads)
        pads_values_arr = np.array(pads_values_flat)
        p_wires_arr = np.array(p_flat_wires)
        wires_values_arr = np.array(wires_values_flat)

        data[label] = {
            "p_pads": p_pads_arr,
            "pads_values": pads_values_arr,
            "p_wires": p_wires_arr,
            "wires_values": wires_values_arr,
        }

    return data


# ======================================================================
#  PARTICLE COLORS AND LABELS
# ======================================================================

particle_info = {
    "ch": {
        "pads_color": "#912A0B",
        "wires_color": "#0B918D",
        "label": r"Charged Hadrons ($\pi^\pm, K^\pm, p$)"
    },
    "el": {
        "pads_color": "#F0175C",
        "wires_color": "#2F0B91",
        "label": r"Electrons ($e^\pm$)"
    },
    "mu": {
        "pads_color": "#FF8C00",
        "wires_color": "#8B008B",
        "label": r"Muons ($\mu^\pm$)"
    }
}

# Flavor information
flavor_info = {
    "Zuu": {"label": r"$u\bar{u}$", "order": 0},
    "Zdd": {"label": r"$d\bar{d}$", "order": 1},
    "Zss": {"label": r"$s\bar{s}$", "order": 2},
    "Zcc": {"label": r"$c\bar{c}$", "order": 3},
    "Zbb": {"label": r"$b\bar{b}$", "order": 4},
    "Data": {"label": r"1994 Data", "order": 5},
}


# ======================================================================
#  GRID PLOTTING FUNCTIONS
# ======================================================================

def plot_grid_by_particle(all_data, particle_type, detector_type):
    """
    Create a 2x3 grid plot showing all flavors for a specific particle type and detector.
    
    Parameters:
    -----------
    all_data : dict
        Dictionary containing data for all flavors
    particle_type : str
        "ch", "el", or "mu"
    detector_type : str
        "pads" or "wires"
    """
    
    # Sort flavors by order
    sorted_flavors = sorted(flavor_info.keys(), key=lambda x: flavor_info[x]["order"])
    
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    axes = axes.flatten()
    
    info = particle_info[particle_type]
    color = info[f"{detector_type}_color"]
    
    for idx, flavor in enumerate(sorted_flavors):
        ax = axes[idx]
        data = all_data[flavor][particle_type]
        
        p_key = f"p_{detector_type}"
        values_key = f"{detector_type}_values"
        
        ax.scatter(
            data[p_key], 
            data[values_key], 
            s=0.3, 
            color=color, 
            alpha=0.4,
            rasterized=True
        )
        
        ax.set_xlabel(r"$p$ [GeV]")
        ax.set_ylabel(r"$\mathrm{d}E/\mathrm{d}x$ " + f"({detector_type.capitalize()})")
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 5)
        
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which="both", direction="in", top=True, right=True)
        
        # Title for each subplot
        ax.text(
            0.05, 0.95,
            flavor_info[flavor]["label"],
            transform=ax.transAxes,
            fontsize=11,
            va="top",
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='none')
        )
    
    # Overall title
    detector_name = detector_type.capitalize()
    particle_label = info["label"]
    fig.suptitle(
        r"\textbf{ALEPH} \textit{Preliminary} -- " + particle_label + f" ({detector_name})",
        fontsize=14,
        y=0.995
    )
    
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    
    output_name = f"dEdx{detector_name}_vs_p_{particle_type}_grid"
    plt.savefig(f"{output_name}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_name}.pdf", dpi=600, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved {output_name}.pdf/.png")


def plot_grid_combined_by_flavor(all_data, flavor, detector_type):
    """
    Create a combined plot showing all particle types for a specific flavor and detector.
    
    Parameters:
    -----------
    all_data : dict
        Dictionary containing data for all flavors
    flavor : str
        Flavor name (e.g., "Zuu", "Data")
    detector_type : str
        "pads" or "wires"
    """
    
    fig, ax = plt.subplots(figsize=(7, 6))
    
    for particle_type in ["ch", "el", "mu"]:
        info = particle_info[particle_type]
        color = info[f"{detector_type}_color"]
        data = all_data[flavor][particle_type]
        
        p_key = f"p_{detector_type}"
        values_key = f"{detector_type}_values"
        
        ax.scatter(
            data[p_key], 
            data[values_key], 
            s=0.5, 
            color=color, 
            alpha=0.3,
            label=info["label"],
            rasterized=True
        )
    
    ax.set_xlabel(r"$p$ [GeV]")
    ax.set_ylabel(r"$\mathrm{d}E/\mathrm{d}x$ " + f"({detector_type.capitalize()})")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="both", direction="in", top=True, right=True)
    
    # Header text
    ax.text(
        0.05, 0.95,
        r"\textbf{ALEPH} \textit{Preliminary}",
        transform=ax.transAxes,
        fontsize=13,
        va="top",
    )
    ax.text(
        0.05, 0.89,
        flavor_info[flavor]["label"] + f" ({detector_type.capitalize()})",
        transform=ax.transAxes,
        fontsize=11,
        va="top",
    )
    
    ax.legend(loc="upper right", frameon=False, markerscale=5)
    
    plt.tight_layout()
    
    output_name = f"dEdx{detector_type.capitalize()}_allParticles_{flavor}"
    plt.savefig(f"{output_name}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_name}.pdf", dpi=600, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved {output_name}.pdf/.png")


def plot_full_grid_combined(all_data, detector_type):
    """
    Create a 2x3 grid showing all particles combined for each flavor.
    
    Parameters:
    -----------
    all_data : dict
        Dictionary containing data for all flavors
    detector_type : str
        "pads" or "wires"
    """
    
    sorted_flavors = sorted(flavor_info.keys(), key=lambda x: flavor_info[x]["order"])
    
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes = axes.flatten()
    
    for idx, flavor in enumerate(sorted_flavors):
        ax = axes[idx]
        
        for particle_type in ["ch", "el", "mu"]:
            info = particle_info[particle_type]
            color = info[f"{detector_type}_color"]
            data = all_data[flavor][particle_type]
            
            p_key = f"p_{detector_type}"
            values_key = f"{detector_type}_values"
            
            # Only add label for first subplot
            label = info["label"] if idx == 0 else None
            
            ax.scatter(
                data[p_key], 
                data[values_key], 
                s=0.3, 
                color=color, 
                alpha=0.3,
                label=label,
                rasterized=True
            )
        
        ax.set_xlabel(r"$p$ [GeV]")
        ax.set_ylabel(r"$\mathrm{d}E/\mathrm{d}x$ " + f"({detector_type.capitalize()})")
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 5)
        
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which="both", direction="in", top=True, right=True)
        
        # Title for each subplot
        ax.text(
            0.05, 0.95,
            flavor_info[flavor]["label"],
            transform=ax.transAxes,
            fontsize=11,
            va="top",
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='none')
        )
    
    # Add legend to first subplot
    axes[0].legend(loc="upper right", frameon=False, markerscale=5, fontsize=8)
    
    # Overall title
    detector_name = detector_type.capitalize()
    fig.suptitle(
        r"\textbf{ALEPH} \textit{Preliminary} -- All Particles " + f"({detector_name})",
        fontsize=14,
        y=0.995
    )
    
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    
    output_name = f"dEdx{detector_name}_allParticles_allFlavors_grid"
    plt.savefig(f"{output_name}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_name}.pdf", dpi=600, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved {output_name}.pdf/.png")


# ======================================================================
#  MAIN
# ======================================================================

if __name__ == "__main__":
    configure_matplotlib()
    
    # Define file paths
    data_file = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedData/1994/stage1/2dEdx/1994.root"
    mc_dir = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/1994/zqq/stage1/2dEdx/"
    
    mc_files = {
        "Zuu": os.path.join(mc_dir, "Zuu.root"),
        "Zdd": os.path.join(mc_dir, "Zdd.root"),
        "Zss": os.path.join(mc_dir, "Zss.root"),
        "Zcc": os.path.join(mc_dir, "Zcc.root"),
        "Zbb": os.path.join(mc_dir, "Zbb.root"),
    }
    
    # Load all data
    print("Loading data files...")
    all_data = {}
    
    # Number of events to use (randomly sampled)
    n_events = 1000
    
    # Load MC files
    for flavor, filepath in mc_files.items():
        print(f"\nLoading {flavor}...")
        all_data[flavor] = load_and_flatten_data(filepath, n_events=n_events, seed=42)
    
    # Load Data
    print(f"\nLoading 1994 Data...")
    all_data["Data"] = load_and_flatten_data(data_file, n_events=n_events, seed=42)
    
    print("\n" + "="*70)
    print("Creating grid plots...")
    print("="*70)
    
    # Create grid plots for each particle type (Pads and Wires)
    for particle_type in ["ch", "el", "mu"]:
        for detector_type in ["pads", "wires"]:
            print(f"\nCreating {particle_type} {detector_type} grid...")
            plot_grid_by_particle(all_data, particle_type, detector_type)
    
    # Create combined plots for each flavor
    for flavor in flavor_info.keys():
        for detector_type in ["pads", "wires"]:
            print(f"\nCreating combined {flavor} {detector_type} plot...")
            plot_grid_combined_by_flavor(all_data, flavor, detector_type)
    
    # Create full grid plots (all particles, all flavors)
    for detector_type in ["pads", "wires"]:
        print(f"\nCreating full grid for {detector_type}...")
        plot_full_grid_combined(all_data, detector_type)
    
    print("\n" + "="*70)
    print("All plots created successfully!")
    print("="*70)
