from collections import namedtuple
import ROOT

ConstituentSelection = namedtuple('ConstituentSelection', ['p_min', 'dedx_min', 'dedx_max', 'label'])

selection_dictionary={
    "all_pfcands":ConstituentSelection(p_min=0., dedx_min=0., dedx_max=0., label="All jet constituents"),
    "dedx_data_mc_quality_no_p_cut":ConstituentSelection(p_min=0., dedx_min=0.5, dedx_max=2.5, label="0.5 < dE/dx_{i} < 2.5 for i in jet constituents"),
    "dedx_data_mc_quality_hi_p_cut":ConstituentSelection(p_min=5., dedx_min=0.5, dedx_max=2.5, label="p_{i} > 5 GeV, 0.5 < dE/dx_{i} < 2.5 for i in jet constituents"),
    "dedx_data_mc_quality":ConstituentSelection(p_min=1., dedx_min=0.5, dedx_max=2.5, label="p_{i} > 1 GeV, 0.5 < dE/dx_{i} < 2.5 for i in jet constituents"),
}