import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os

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

files = {
    "LEP1-1994 (Class 16–17)": {"path": "1994.root", "type": "data"},
    r"$Z \to d\bar{d}$": {"path": "Zdd.root", "type": "mc"},
    r"$Z \to u\bar{u}$": {"path": "Zuu.root", "type": "mc"},
    r"$Z \to s\bar{s}$": {"path": "Zss.root", "type": "mc"},
    r"$Z \to c\bar{c}$": {"path": "Zcc.root", "type": "mc"},
    r"$Z \to b\bar{b}$": {"path": "Zbb.root", "type": "mc"},
}

colors = ['blue', 'green', 'orange', 'purple', 'brown']

L_data = 42.64
mc_xsec = {
    r"$Z \to b\bar{b}$": 2.81e5,
    r"$Z \to c\bar{c}$": 2.81e5 * 0.8,
    r"$Z \to s\bar{s}$": 2.81e5,
    r"$Z \to d\bar{d}$": 2.81e5,
    r"$Z \to u\bar{u}$": 2.81e5 * 0.8,
}

stack_order = [
    r"$Z \to d\bar{d}$",
    r"$Z \to u\bar{u}$",
    r"$Z \to s\bar{s}$",
    r"$Z \to c\bar{c}$",
    r"$Z \to b\bar{b}$"
]


def plot_variable_stacked(var, bins=50, logy=False, save=False, xmin=None, xmax=None):
    print(f"📊 Plotting stacked: {var}")

    mc_data = {}
    data_arr = None

    # Load input
    for label, info in files.items():
        with uproot.open(info["path"]) as f:
            arr = f["events"][var].array(library="ak")

        # ✅ FIX: flatten ALL jagged structures
        arr = ak.flatten(arr, axis=None)
        arr = ak.to_numpy(arr)

        # ✅ Remove invalid values
        arr = arr[np.isfinite(arr)]

        if info["type"] == "data":
            data_arr = arr
        else:
            mc_data[label] = arr

    if data_arr is None:
        raise RuntimeError("❌ No data found")

    if xmin is None: xmin = np.min(data_arr)
    if xmax is None: xmax = np.max(data_arr)

    edges = np.linspace(xmin, xmax, bins + 1)
    centers = 0.5*(edges[1:]+edges[:-1])

    # ✅ MC scaling
    mc_counts = []
    mc_weights = []
    mc_colors = []

    for i, label in enumerate(stack_order):
        arr = mc_data[label]
        sigma = mc_xsec[label]
        N = len(arr)
        weight = (L_data * sigma) / N
        mc_counts.append(arr)
        mc_weights.append(np.ones_like(arr) * weight)
        mc_colors.append(colors[i % len(colors)])

    total_mc_raw = sum(np.sum(w) for w in mc_weights)
    scale = len(data_arr) / total_mc_raw
    mc_weights = [w * scale for w in mc_weights]

    fig = plt.figure(figsize=(14, 14), dpi=300)
    gs = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)
    ax = fig.add_subplot(gs[0])
    ax_ratio = fig.add_subplot(gs[1], sharex=ax)

    # ✅ Stacked MC hist
    ax.hist(mc_counts, bins=edges, stacked=True,
            weights=mc_weights, color=mc_colors,
            edgecolor='black', label=stack_order)

    # ✅ Data
    y_data, _ = np.histogram(data_arr, bins=edges)
    y_err = np.sqrt(y_data)
    ax.errorbar(centers, y_data, yerr=y_err, fmt='o', color='black',
                label="LEP1-1994 (Class 16–17)")

    # ✅ Ratio
    mc_total, _ = np.histogram(np.concatenate(mc_counts), bins=edges,
                               weights=np.concatenate(mc_weights))
    ratio = y_data / mc_total
    ratio_err = y_err / mc_total

    ax_ratio.errorbar(centers, ratio, yerr=ratio_err,
                      fmt='o', color='black')
    ax_ratio.axhline(1.0, linestyle="--", color="red")

    ax.set_ylabel("Entries", fontweight="bold")
    ax_ratio.set_ylabel("Data/MC", fontweight="bold")
    ax_ratio.set_xlabel(var.replace("_"," ").capitalize(), fontweight="bold")

    ax.legend(frameon=False)
    if logy: ax.set_yscale("log")
    ax.grid(alpha=0.3)
    ax_ratio.grid(alpha=0.3)
    plt.setp(ax.get_xticklabels(), visible=False)
    ax_ratio.set_ylim(0, 4)

    if save:
        os.makedirs("pfcand", exist_ok=True)
        fname = f"pfcand/{var}_stacked.png"
        plt.savefig(fname, bbox_inches="tight")
        print(f"✅ Saved: {fname}")

    plt.show()
    print(f"✅ Total Data = {len(data_arr)}, Total scaled MC = {int(np.sum(np.concatenate(mc_weights)))}")


plot_variable_stacked("pfcand_dxy", bins=100,logy=True, save=True, xmin=-1, xmax=1, )#density=True)
plot_variable_stacked("pfcand_dz", bins=100,logy=True, save=True, xmin=-3, xmax=3,)# density=True)
