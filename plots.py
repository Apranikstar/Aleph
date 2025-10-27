import ROOT
import os
import argparse
import math

# ________________________________________________________________________________
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--indir",
        help="path input directory",
        default="/eos/experiment/fcc/ee/analyses/case-studies/aleph/mc/zqq/stage2/v3/",
    )
    parser.add_argument(
        "--outdir",
        help="path output directory",
        default="/eos/user/s/selvaggi/www/test_alephtag_v3",
    )

    args = parser.parse_args()

    # Enable multi-threading
    ROOT.ROOT.EnableImplicitMT()
    ROOT.gROOT.SetBatch(True)

    from config_aleph import variables_pfcand, variables_jet, variables_event, flavors

    input_dir = args.indir
    output_dir = args.outdir

    os.system("mkdir -p {}".format(output_dir))

    # for f in flavors:

    #     sample_a = {
    #         "file": "{}/Z{}{}.root".format(input_dir, f, f),
    #         "flavor": f,
    #         "label": "Aleph MC",
    #     }
    #     sample_b = {
    #         "file": "{}/Z{}{}.root".format(input_dir, f, f),
    #         "flavor": f,
    #         "label": "ALEPH",
    #     }
    #     # We read the tree from the file and create a RDataFrame.
    #     df_a = ROOT.RDataFrame("tree", sample_a["file"])
    #     df_b = ROOT.RDataFrame("tree", sample_b["file"])

    #     print(sample_a["file"])
    #     print(sample_b["file"])

    #     sample_a["histos_pfcand"] = dfhs_pfcand(df_a, variables_pfcand)
    #     sample_b["histos_pfcand"] = dfhs_pfcand(df_b, variables_pfcand)

    #     sample_a["histos_jet"] = dfhs_jet(df_a, variables_jet)
    #     sample_b["histos_jet"] = dfhs_jet(df_b, variables_jet)

    #     #sample_a["histos_event"] = dfhs_event(df_a, variables_event)
    #     #sample_b["histos_event"] = dfhs_event(df_b, variables_event)

    #     # RunGraphs allows to run the event loops of the separate RDataFrame graphs
    #     # concurrently. This results in an improved usage of the available resources
    #     # if each separate RDataFrame can not utilize all available resources, e.g.,
    #     ROOT.RDF.RunGraphs(
    #         # list(sample_a["histos_pfcand"].values())
    #         # + list(sample_b["histos_pfcand"].values())
    #         list(sample_a["histos_jet"].values())
    #         + list(sample_b["histos_jet"].values())
    #         #+ list(sample_a["histos_event"].values())
    #         #+ list(sample_b["histos_event"].values())
    #     )


    #     for var, params in variables_pfcand.items():
    #         plot(sample_a, sample_b, "histos_pfcand", var, params, output_dir)
    #     for var, params in variables_jet.items():
    #         plot(sample_a, sample_b, "histos_jet", var, params, output_dir)

    #     """
    #     for var, params in variables_event.items():
    #         plot(sample_a, sample_b, "histos_event", var, params, output_dir)
    #     """


    version_label = "ALEPH"  # change if you want another label

    samples = []
    for f in flavors:
        file_path = f"{input_dir}/Z{f}{f}.root"
        print(file_path)
        df = ROOT.RDataFrame("tree", file_path)

        sample = {
            "file": file_path,
            "flavor": f,
            "label": version_label,
        }
        sample["histos_pfcand"] = dfhs_pfcand(df, variables_pfcand)
        sample["histos_jet"]    = dfhs_jet(df,    variables_jet)
        # sample["histos_event"]  = dfhs_event(df,  variables_event)

        samples.append(sample)

    # Run all graphs together (both groups so histos are materialized)
    ROOT.RDF.RunGraphs(
        sum(
            [
                list(s["histos_pfcand"].values())
                + list(s["histos_jet"].values())
                # + list(s["histos_event"].values())
                for s in samples
            ],
            []
        )
    )

    # Plot: one figure per variable, all flavors overlaid, ratios to first flavor
    for var, params in variables_pfcand.items():
        plot_multi(samples, "histos_pfcand", var, params, output_dir)

    for var, params in variables_jet.items():
        plot_multi(samples, "histos_jet", var, params, output_dir)

    """
    for var, params in variables_event.items():
        plot_multi(samples, "histos_event", var, params, output_dir)
    """

# _______________________________________________________________________________
def dfhs_pfcand(df, vars):

    ## extract charged particles
    # df_charged = df.Filter("All(abs(pfcand_charge)>0)", "select charged constituents")
    df_charged = df

    ## order constituents in energy
    df_sorted_e = df_charged.Define("e_sorted_id", "Reverse(Argsort(pfcand_e))")

    df_dict = dict()

    for pfcand_var, params in vars.items():
        df_var = df_sorted_e.Redefine(pfcand_var, "Take({}, e_sorted_id)".format(pfcand_var))
        var = pfcand_var.replace("pfcand_", "")
        df_var = df_var.Define(var, "{}[0]".format(pfcand_var))
        df_dict[pfcand_var] = df_var.Histo1D(
            (
                "h_{}".format(var),
                ";{};N_{{Events}}".format(params["title"]),
                params["bin"],
                params["xmin"],
                params["xmax"],
            ),
            var,
        )
    return df_dict


# _______________________________________________________________________________
def dfhs_jet(df, vars):

    ## extract charged particles
    # df_charged = df.Filter("All(abs(pfcand_charge)>0)", "select charged constituents")
    df_dict = dict()
    for jet_var, params in vars.items():
        df_dict[jet_var] = df.Histo1D(
            (
                "h_{}".format(jet_var),
                ";{};N_{{Events}}".format(params["title"]),
                params["bin"],
                params["xmin"],
                params["xmax"],
            ),
            jet_var,
        )
    return df_dict


# _______________________________________________________________________________
def dfhs_event(df, vars):

    ## extract charged particles
    # df_charged = df.Filter("All(abs(pfcand_charge)>0)", "select charged constituents")
    df_dict = dict()
    for event_var, params in vars.items():
        print(event_var)
        df_dict[event_var] = df.Histo1D(
            (
                "h_{}".format(event_var),
                ";{};N_{{Events}}".format(params["title"]),
                params["bin"],
                params["xmin"],
                params["xmax"],
            ),
            event_var,
        )
    return df_dict


# _______________________________________________________________________________
def plot(sample_a, sample_b, histo_coll, var, params, outdir):

    dfh_a = sample_a[histo_coll][var].GetValue()
    dfh_b = sample_b[histo_coll][var].GetValue()

    # Create canvas with pads for main plot and data/MC ratio
    c = ROOT.TCanvas("c", "", 700, 750)

    ROOT.gStyle.SetOptStat(0)
    upper_pad = ROOT.TPad("upper_pad", "", 0, 0.35, 1, 1)
    lower_pad = ROOT.TPad("lower_pad", "", 0, 0, 1, 0.35)
    for p in [upper_pad, lower_pad]:
        p.SetLeftMargin(0.14)
        p.SetRightMargin(0.05)
        p.SetTickx(False)
        p.SetTicky(False)
    upper_pad.SetBottomMargin(0)
    lower_pad.SetTopMargin(0)
    lower_pad.SetBottomMargin(0.3)
    upper_pad.Draw()
    lower_pad.Draw()

    # Draw dfh_a
    upper_pad.cd()
    if params["scale"] == "log":
        upper_pad.SetLogy()
    dfh_a.SetMarkerStyle(20)
    dfh_a.SetMarkerSize(0)
    dfh_a.SetLineWidth(4)
    dfh_a.SetLineColor(ROOT.kGreen + 2)
    dfh_a.GetYaxis().SetLabelSize(0.045)
    dfh_a.GetYaxis().SetTitleSize(0.05)
    dfh_a.SetStats(0)
    dfh_a.SetTitle("")
    dfh_a.Draw("hist")

    # Draw dfh_b
    dfh_b.SetLineColor(ROOT.kRed + 1)
    dfh_b.SetLineStyle(2)
    dfh_b.SetLineWidth(4)
    dfh_b.Draw("hist SAME")

    # Draw ratio
    lower_pad.cd()

    ratio = ROOT.TH1I(
        "zero",
        "",
        params["bin"],
        params["xmin"],
        params["xmax"],
    )
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetLineStyle(2)
    ratio.SetLineWidth(4)
    ratio.SetMinimum(0.0)
    ratio.SetMaximum(2.0)
    ratio.GetXaxis().SetLabelSize(0.08)
    ratio.GetXaxis().SetTitleSize(0.12)
    ratio.GetXaxis().SetTitleOffset(1.0)
    ratio.GetYaxis().SetLabelSize(0.08)
    ratio.GetYaxis().SetTitleSize(0.09)
    ratio.GetYaxis().SetTitle("ratio")
    ratio.GetYaxis().CenterTitle()
    ratio.GetYaxis().SetTitleOffset(0.7)
    # ratio.GetYaxis().SetNdivisions(503, False)
    ratio.GetYaxis().ChangeLabel(-1, -1, 0)
    ratio.GetXaxis().SetTitle(params["title"])
    ratio.Draw("AXIS")

    ratiodata = dfh_a.Clone()
    ratiodata.Sumw2()
    ratiodata.Divide(dfh_b)
    ratiodata.SetLineColor(ROOT.kBlack)
    ratiodata.SetMarkerColor(ROOT.kBlack)
    ratiodata.Draw("same e")

    # Add legend
    upper_pad.cd()
    legend = ROOT.TLegend(0.55, 0.68, 0.926, 0.85)
    legend.SetTextFont(42)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.045)
    legend.SetTextAlign(12)
    legend.AddEntry(dfh_a, "{} ({}-jets)".format(sample_a["label"], sample_a["flavor"]), "l")
    legend.AddEntry(dfh_b, "{} ({}-jets)".format(sample_b["label"], sample_b["flavor"]), "l")
    legend.Draw()

    # Add ALEPH label
    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(72)
    text.SetTextSize(0.05)
    text.DrawLatex(0.14, 0.91, "ALEPH")
    text.SetTextFont(42)
    text.DrawLatex(0.27, 0.91, "(ALEPH Simulation)")
    text.SetTextSize(0.05)
    text.DrawLatex(0.25, 0.78, "e^{+}e^{-} #rightarrow j j")
    text.SetTextSize(0.04)
    text.DrawLatex(0.28, 0.71, "j = u, d, s, c, b")
    text.SetTextSize(0.05)
    text.SetTextAlign(31)
    #text.DrawLatex(0.95, 0.91, "#sqrt{s} = 240 GeV, 5 ab^{-1}")

    # Save the plot
    figpath = "{}/{}_{}.png".format(outdir, sample_a["flavor"], var)
    c.SaveAs(figpath)



# _______________________________________________________________________________
def plot_multi(
    samples,          # list of sample dicts: [{"label": "...", "flavor": "...", histo_coll: {var: RResultPtr<TH1>}}, ...]
    histo_coll,       # e.g. "histos_pfcand" | "histos_jet" | "histos_event"
    var,              # variable name inside the histo_coll dict
    params,           # dict with at least: "bin", "xmin", "xmax", "title", "scale" ("lin"|"log")
    outdir,           # output directory
    colors=None,      # optional list of ROOT colors
    linestyles=None,  # optional list of ROOT line styles
    ratio_range=(0.5, 1.5)  # default ratio panel y-range; will auto-expand if needed
):
    """
    Multi-curve version of your plot():
    - Upper pad: overlays all histograms.
    - Lower pad: draws ratio for every curve i>0 as (hist_i / hist_ref), where hist_ref is the FIRST sample.
    - Legend shows "<label> (<flavor>-jets)" per curve.

    NOTE: If you prefer ref/others instead, swap numerator/denominator where indicated below.
    """

    if not samples or len(samples) < 2:
        raise ValueError("plot_multi requires at least 2 samples (first is the reference).")

    # --- materialize histos (keep original style: GetValue()) ---
    def get_h(sample):
        h = sample[histo_coll][var].GetValue()
        # Ensure it has sumw2 for correct ratio errors
        if not h.GetSumw2N():
            h.Sumw2()
        return h

    H = [get_h(s) for s in samples]
    href = H[0]

    # --- canvas & pads (keeps your layout/style) ---
    c = ROOT.TCanvas("c", "", 700, 750)
    ROOT.gStyle.SetOptStat(0)

    upper_pad = ROOT.TPad("upper_pad", "", 0, 0.35, 1, 1)
    lower_pad = ROOT.TPad("lower_pad", "", 0, 0, 1, 0.35)
    for p in (upper_pad, lower_pad):
        p.SetLeftMargin(0.14)
        p.SetRightMargin(0.05)
        p.SetTickx(False)
        p.SetTicky(False)
    upper_pad.SetBottomMargin(0)
    lower_pad.SetTopMargin(0)
    lower_pad.SetBottomMargin(0.3)

    upper_pad.Draw()
    lower_pad.Draw()

    # --- style helpers ---
    if colors is None:
        colors = [
            ROOT.kGreen+2, ROOT.kRed+1, ROOT.kAzure+1, ROOT.kMagenta+2,
            ROOT.kOrange+7, ROOT.kCyan+2, ROOT.kViolet+1, ROOT.kTeal+4
        ]
    if linestyles is None:
        linestyles = [1, 1,1,1,1,1,1,1,1]

    # --- upper: draw reference first ---
    upper_pad.cd()
    if params.get("scale", "lin") == "log":
        upper_pad.SetLogy()

    # ensure positive minimum for log scale
    if params.get("scale", "lin") == "log":
        # find a minimal positive bin content across all histos
        min_pos = None
        for h in H:
            for b in range(1, h.GetNbinsX()+1):
                val = h.GetBinContent(b)
                if val > 0 and (min_pos is None or val < min_pos):
                    min_pos = val
        if min_pos is None:
            min_pos = 1e-6
        href.SetMinimum(min_pos * 0.5)

    href.SetMarkerStyle(20)
    href.SetMarkerSize(0)
    href.SetLineWidth(4)
    href.SetLineColor(colors[0])
    href.SetLineStyle(linestyles[0] if len(linestyles) > 0 else 1)
    href.GetYaxis().SetLabelSize(0.045)
    href.GetYaxis().SetTitleSize(0.05)
    href.SetStats(0)
    href.SetTitle("")
    href.Draw("hist")

    # --- draw the others on top ---
    for i, h in enumerate(H[1:], start=1):
        h.SetLineColor(colors[i % len(colors)])
        h.SetLineStyle(linestyles[i % len(linestyles)])
        h.SetLineWidth(4)
        h.Draw("hist SAME")

    # --- legend ---
    legend = ROOT.TLegend(0.55, 0.68, 0.926, 0.85)
    legend.SetTextFont(42)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.045)
    legend.SetTextAlign(12)

    for i, s in enumerate(samples):
        h = H[i]
        entry = legend.AddEntry(h, f'{s["label"]} ({s["flavor"]}-jets)', "l")
        # (ROOT uses the histogram's current style/colors for legend swatch)
    legend.Draw()

    # --- FCC-ee / process label (unchanged from your code) ---
    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(72)
    text.SetTextSize(0.05)
    text.DrawLatex(0.14, 0.91, "ALEPH")
    text.SetTextFont(42)
    text.DrawLatex(0.27, 0.91, "(ALEPH Simulation)")
    text.SetTextSize(0.05)
    text.DrawLatex(0.25, 0.78, "e^{+}e^{-} #rightarrow j j")
    text.SetTextSize(0.04)
    text.DrawLatex(0.28, 0.71, "j = u, d, s, c, b")
    text.SetTextSize(0.05)
    text.SetTextAlign(31)
    #text.DrawLatex(0.95, 0.91, "#sqrt{s} = 240 GeV, 5 ab^{-1}")

    # --- lower: ratio axes frame ---
    lower_pad.cd()
    frame = ROOT.TH1I("ratio_frame", "", params["bin"], params["xmin"], params["xmax"])
    frame.SetLineColor(ROOT.kBlack)
    frame.SetLineStyle(2)
    frame.SetLineWidth(4)
    frame.GetXaxis().SetLabelSize(0.08)
    frame.GetXaxis().SetTitleSize(0.12)
    frame.GetXaxis().SetTitleOffset(1.0)
    frame.GetXaxis().SetTitle(params["title"])
    frame.GetYaxis().SetLabelSize(0.08)
    frame.GetYaxis().SetTitleSize(0.09)
    frame.GetYaxis().SetTitle("ratio to first")
    frame.GetYaxis().CenterTitle()
    frame.GetYaxis().SetTitleOffset(0.7)
    frame.GetYaxis().ChangeLabel(-1, -1, 0)

    # Compute ratios to set a good y-range
    ratios = []
    rmin, rmax = ratio_range
    for i, h in enumerate(H[1:], start=1):
        r = h.Clone(f"ratio_{i}")
        r.Sumw2()
        # others / first (swap order if you want ref/others)
        r.Divide(href)
        ratios.append((i, r))
        # scan finite bin contents for dynamic range
        nb = r.GetNbinsX()
        for b in range(1, nb+1):
            val = r.GetBinContent(b)
            err = r.GetBinError(b)
            if math.isfinite(val) and val > 0:
                rmin = min(rmin, val - err)
                rmax = max(rmax, val + err)

    # pad a bit
    if rmax <= rmin:
        rmin, rmax = 0.8, 1.2
    margin = 0.05 * (rmax - rmin)
    frame.SetMinimum(rmin - margin)
    frame.SetMaximum(rmax + margin)
    frame.Draw("AXIS")

    # Draw a unity reference line
    unity = ROOT.TLine(params["xmin"], 1.0, params["xmax"], 1.0)
    unity.SetLineStyle(7)
    unity.SetLineColor(ROOT.kGray+2)
    unity.Draw()

    # Draw ratios
    for i, r in ratios:
        r.SetLineColor(colors[i % len(colors)])
        r.SetLineStyle(linestyles[i % len(linestyles)])
        r.SetLineWidth(3)
        r.SetMarkerStyle(20)
        r.SetMarkerSize(0)
        r.Draw("same hist e")

    # --- save ---
    # name based on var and first sample flavor as reference
    ref_tag = f'{samples[0]["flavor"]}'
    figpath = f"{outdir}/{ref_tag}_{var}_multi.png"
    c.SaveAs(figpath)


if __name__ == "__main__":
    main()