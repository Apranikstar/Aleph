import os
import uproot
import math
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

def process_root_file_to_parquet(
    filepath: str,
    file_name: str,
    tree: str,
    output_dir: str,
    chunk_size: int = 100_000
):
    # --- Ensure output directory exists ---
    os.makedirs(output_dir, exist_ok=True)

    pfcand_keys = [
        "pfcand_dEdx_pads_type", "pfcand_dEdx_pads_value", "pfcand_dEdx_pads_error",
        "pfcand_dEdx_wires_type", "pfcand_dEdx_wires_value", "pfcand_dEdx_wires_error",
               
                
        "pfcand_isMu", "pfcand_isEl", "pfcand_isChargedHad", "pfcand_isGamma", "pfcand_isNeutralHad",
        "pfcand_e", "pfcand_p", "pfcand_theta", "pfcand_phi", "pfcand_charge", "pfcand_type",
        "pfcand_erel", "pfcand_erel_log", "pfcand_thetarel", "pfcand_phirel", 

        "pfcand_dxy", "pfcand_dz", "pfcand_phi0", "pfcand_C", "pfcand_ct",
        "pfcand_dptdpt", "pfcand_dxydxy", "pfcand_dzdz", "pfcand_dphidphi", "pfcand_detadeta",
        "pfcand_dxydz", "pfcand_dphidxy", "pfcand_phidz", "pfcand_phictgtheta", "pfcand_dxyctgtheta",
        "pfcand_dlambdadz", "pfcand_cctgtheta", "pfcand_phic", "pfcand_dxyc", "pfcand_cdz",
 
        "pfcand_btagSip2dVal", "pfcand_btagSip2dSig", "pfcand_btagSip3dVal", "pfcand_btagSip3dSig", 
        "pfcand_btagJetDistVal", "pfcand_btagJetDistSig",
    ]

    jet_branches = [
        "jet_mass","jet_e",
        "jet_phi", "jet_theta", "jet_eta", 
        "jet_p", "jet_px", "jet_py", "jet_pz", "jet_pT",
        "jet_nconst",
        "jet_nnhad","jet_nchad",
        "jet_nel", "jet_nmu", 
        "jet_ngamma",

    ]

    # --- Open ROOT file ---
    file_path = os.path.join(filepath, file_name)
    file = uproot.open(file_path)[tree]
    num_events = file.num_entries

    leading_jet_index = 0
    subleading_jet_index = 1

    num_chunks = math.ceil(num_events / chunk_size)

    for chunk_idx in range(num_chunks):
        start = chunk_idx * chunk_size
        stop = min((chunk_idx + 1) * chunk_size, num_events)
        print(f"Processing {file_name}: events {start} → {stop-1}")

        data_dict = {}

        # --- Per-constituent data (list per jet) ---
        for key in pfcand_keys:
            arr = file[key].array(library="np", entry_start=start, entry_stop=stop)
            jets_list = []
            for i in range(len(arr)):
                jets_list.append(arr[i][leading_jet_index].tolist())
                jets_list.append(arr[i][subleading_jet_index].tolist())
            data_dict[key] = jets_list

        # --- Jet-level scalar data ---
        for branch in jet_branches:
            arr = file[branch].array(library="np", entry_start=start, entry_stop=stop)
            branch_list = []
            for i in range(len(arr)):
                branch_list.append(arr[i][leading_jet_index])
                branch_list.append(arr[i][subleading_jet_index])
            data_dict[branch] = branch_list

        # --- Labels ---
        num_jets_total = (stop - start) * 2
        labels = {
    "recojet_isU": 0, "recojet_isD": 0, "recojet_isS": 0,
    "recojet_isC": 0, "recojet_isB": 0
        }


        # Determine label based on file name
        for label in labels:
            labels[label] = [0] * num_jets_total

        if "Zuu" in file_name:
            labels["recojet_isU"] = [1] * num_jets_total
        elif "Zdd" in file_name:
            labels["recojet_isD"] = [1] * num_jets_total
        elif "Zss" in file_name:
            labels["recojet_isS"] = [1] * num_jets_total
        elif "Zcc" in file_name:
            labels["recojet_isC"] = [1] * num_jets_total
        elif "Zbb" in file_name:
            labels["recojet_isB"] = [1] * num_jets_total


        data_dict.update(labels)

        # --- Convert to Arrow Table and save ---
        table = pa.Table.from_pydict(data_dict)
        output_path = os.path.join(
            output_dir,
            f"{os.path.splitext(file_name)[0]}_chunk{chunk_idx}.parquet"
        )
        pq.write_table(table, output_path)
        print(f"Saved {output_path}")

# --- Example usage ---
if __name__ == "__main__":
    tree = "events"
    input_dir = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/mc/zqq/stage1/v2/database"
    output_dir = "/afs/cern.ch/work/h/hfatehi/Aleph/prod/stage2"

    for fname in os.listdir(input_dir):
        if fname.endswith(".root"):
            process_root_file_to_parquet(
                filepath=input_dir,
                file_name=fname,
                tree=tree,
                output_dir=output_dir,
                chunk_size=100_000
            )
