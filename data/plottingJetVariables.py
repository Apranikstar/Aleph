
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os

# --- Matplotlib config ---
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "axes.labelsize": 15,
    "axes.titlesize": 16,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 12,
    "figure.dpi": 300,
    "lines.linewidth": 1.5,
})

# --- Files and MC info ---
files = {
    r"$Z \to d\bar{d}$": {"path": "Zdd.root", "type": "mc"},
    r"$Z \to u\bar{u}$": {"path": "Zuu.root", "type": "mc"},
        r"$Z \to s\bar{s}$": {"path": "Zss.root", "type": "mc"},
        r"$Z \to c\bar{c}$": {"path": "Zcc.root", "type": "mc"},
        r"$Z \to b\bar{b}$": {"path": "Zbb.root", "type": "mc"},

    "LEP1-1994 (Class 16–17)": {"path": "1994.root", "type": "data"},

}

colors = ['blue', 'green', 'orange', 'purple', 'brown']
linestyles = ['-', '--', '-.', ':', (0, (3, 3, 1, 3))]

L_data = 42.64  # pb^-1
mc_xsec = {
    r"$Z \to b\bar{b}$": 2.81e5,
    r"$Z \to c\bar{c}$": 2.81e5 * 0.8,
    r"$Z \to s\bar{s}$": 2.81e5,
    r"$Z \to d\bar{d}$": 2.81e5,
    r"$Z \to u\bar{u}$": 2.81e5 * 0.8,
}



def plot_variable_stacked_custom(var, bins=50, logy=False, save=False, xmin=None, xmax=None, density=False):
    print(f"📊 Plotting (custom stacked): {var}")

    mc_data = {}
    data_arr = None

    # --- Load data ---
    for label, info in files.items():
        with uproot.open(info["path"]) as f:
            arr = f["events"][var].array(library="ak")
        if isinstance(arr, ak.Array) and arr.ndim > 1:
            arr = ak.flatten(arr)
        arr = np.array(arr, dtype=float)
        arr = arr[np.isfinite(arr)]
        if info["type"] == "data":
            data_arr = arr
        else:
            mc_data[label] = arr

    if data_arr is None:
        raise RuntimeError("❌ No data sample found.")

    # --- Histogram range ---
    if xmin is None:
        xmin = np.min(data_arr)
    if xmax is None:
        xmax = np.max(data_arr)
    hist_range = (xmin, xmax)
    edges = np.linspace(xmin, xmax, bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # --- Define custom stacking order ---
    stack_order = [r"$Z \to d\bar{d}$", r"$Z \to u\bar{u}$", r"$Z \to s\bar{s}$",
                   r"$Z \to c\bar{c}$", r"$Z \to b\bar{b}$"]
    mc_counts = []
    mc_labels = []
    mc_colors = []
    mc_weights_list = []

    for i, label in enumerate(stack_order):
        arr = mc_data[label]
        sigma = mc_xsec[label]
        N_gen = len(arr)
        w = (L_data * sigma) / N_gen
        mc_weights_list.append(np.ones_like(arr) * w)
        mc_counts.append(arr)
        mc_labels.append(label)
        mc_colors.append(colors[i % len(colors)])

    # --- Global scale factor ---
    total_mc_weighted = sum(np.sum(w) for w in mc_weights_list)
    scale_factor = len(data_arr) / total_mc_weighted
    mc_weights_list = [w * scale_factor for w in mc_weights_list]

    # --- Figure setup ---
    fig = plt.figure(figsize=(14, 14), dpi=300)
    gs = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)
    ax = fig.add_subplot(gs[0])
    ax_ratio = fig.add_subplot(gs[1], sharex=ax)

    # --- Plot stacked MC ---
    ax.hist(mc_counts, bins=edges, weights=mc_weights_list, stacked=True,
            color=mc_colors, label=mc_labels, edgecolor='black')

    # --- Plot data ---
    y_data_counts, _ = np.histogram(data_arr, bins=edges)
    y_err = np.sqrt(y_data_counts)
    ax.errorbar(centers, y_data_counts, yerr=y_err, fmt='o', color='black',
                label="LEP1-1994 (Class 16–17)")

    # --- Stats box ---
    text_stats = []
    text_stats.append(f"Data: N={len(data_arr):,}  μ={np.mean(data_arr):.3f}  σ={np.std(data_arr):.3f}")
    for i, arr in enumerate(mc_counts):
        counts_sum_scaled = int(np.sum(mc_weights_list[i]))
        mean = np.average(arr, weights=np.ones_like(arr) * mc_xsec[mc_labels[i]] / len(arr))
        std = np.sqrt(np.average((arr - mean) ** 2, weights=np.ones_like(arr) * mc_xsec[mc_labels[i]] / len(arr)))
        text_stats.append(f"{mc_labels[i]}: Σw={counts_sum_scaled:,}  μ={mean:.3f}  σ={std:.3f}")
    stats_line = "     |     ".join(text_stats)
    ax.text(0.5, 1.02, stats_line,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment="bottom",
            horizontalalignment="center",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.85))

    # --- Decorations ---
    ylabel_text = "Density" if density else "Entries"
    ax.set_ylabel(ylabel_text, fontweight="bold")
    ax.legend(frameon=False)
    if logy:
        ax.set_yscale("log")
    ax.grid(alpha=0.3)

    # --- Ratio plot ---
    mc_sum_counts, _ = np.histogram(np.concatenate(mc_counts), bins=edges, weights=np.concatenate(mc_weights_list))
    mask = mc_sum_counts > 0
    ratio = np.zeros_like(mc_sum_counts)
    ratio_err = np.zeros_like(mc_sum_counts)
    ratio[mask] = y_data_counts[mask] / mc_sum_counts[mask]
    ratio_err[mask] = y_err[mask] / mc_sum_counts[mask]

    ax_ratio.errorbar(centers, ratio, yerr=ratio_err, fmt='o', color='black')
    ax_ratio.axhline(1.0, color='red', linestyle='--', linewidth=1)
    xlabel = " ".join(word.capitalize() for word in var.replace("_", " ").split())
    ax_ratio.set_ylabel("Data/MC", fontweight="bold")
    ax_ratio.set_xlabel(xlabel, fontweight="bold")
    ax_ratio.grid(alpha=0.3)
    plt.setp(ax.get_xticklabels(), visible=True)
    ax_ratio.set_ylim(0, 2)


    # --- Save ---
    if save:
        os.makedirs("jets", exist_ok=True)
        fname = f"jets/{var}_PRD_stacked_custom.png"
        plt.savefig(fname, bbox_inches="tight")
        print(f"✅ Saved: {fname}")

    total_mc_scaled = int(np.sum(np.concatenate(mc_weights_list)))
    print(f"Total scaled MC = {total_mc_scaled}, Total data = {len(data_arr)}")

    plt.show()


plot_variable_stacked_custom("event_invariant_mass", bins=50, logy=False, save=True, xmin=50, xmax=120, density= True)
plot_variable_stacked_custom("jet_e", bins=90, logy=False, save=True, xmin=0, xmax=90, density= True)
plot_variable_stacked_custom("jet_p", bins=90, logy=False, save=True, xmin=0, xmax=90, density= True)
plot_variable_stacked_custom("jet_pT", bins=90, logy=False, save=True, xmin=0, xmax=90, density= True)
plot_variable_stacked_custom("jet_ngamma", bins=25, logy=False, save=True, xmin=0, xmax=25, density= True)
plot_variable_stacked_custom("jet_nnhad", bins=15, logy=False, save=True, xmin=0, xmax=15, density= True)
plot_variable_stacked_custom("jet_nchad", bins=25, logy=False, save=True, xmin=0, xmax=25, density= True)
plot_variable_stacked_custom("jet_nel", bins=5, logy=False, save=True, xmin=0, xmax=5, density= True)
plot_variable_stacked_custom("jet_nmu", bins=5, logy=False, save=True, xmin=0, xmax=5, density= True)
