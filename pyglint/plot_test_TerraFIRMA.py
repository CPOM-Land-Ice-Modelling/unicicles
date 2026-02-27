import os
import json
import yaml
import pickle
import xarray as xr
import matplotlib.pyplot as plt

OUTPUT_DIR = "test_outputs"
ICESHEET = "GrIS"   # AIS or GrIS
VARIABLE_TO_PLOT = "total_smb"

def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def main():

    ## Load the data from Glint/BISICLES for comparison
    if ICESHEET == "AIS":
        with open('/gws/ssde/j25b/terrafirma/tm17544/TerraFIRMA_overshoots/processed_data/AIS_data_overview.pkl', 'rb') as file:
            icesheet_d = pickle.load(file) 

        # only keep the first 50 years of data for cs568 to match the pyglint processed data (which only has 50 years of data)
        icesheet_d["cs568"][0] = icesheet_d["cs568"][0].iloc[:9]

    elif ICESHEET == "GrIS":
        with open('/gws/ssde/j25b/terrafirma/tm17544/TerraFIRMA_overshoots/processed_data/GrIS_data_overview.pkl', 'rb') as file:
            icesheet_d = pickle.load(file)

        # Drop the first year of data as its from an errant file
        icesheet_d["cs568"][0] = icesheet_d["cs568"][0].iloc[1:10]

    # Get the data in the right units/variable for plotting
    bike_total_smb = (icesheet_d["cs568"][0]["grounded_SMB"] + icesheet_d["cs568"][0]["floating_SMB"])*918*1e-12

    plt.figure(figsize=(10, 6))

    for run_folder in sorted(os.listdir(OUTPUT_DIR)):

        run_path = os.path.join(OUTPUT_DIR, run_folder)

        if not os.path.isdir(run_path):
            continue

        config_path = os.path.join(run_path, "config_used.yaml")
        metadata_path = os.path.join(run_path, "metadata.json")
        results_path = os.path.join(run_path, "pyglint_smb_data.nc")

        if not os.path.exists(results_path):
            continue

        cfg = load_yaml(config_path)
        meta = load_json(metadata_path)

        # Filter by icesheet
        if cfg.get("icesheet") != ICESHEET:
            continue

        # Don't plot masked data for now
        if cfg.get("gr_fl_mask"):
            continue

        # Load timeseries
        ds = xr.load_dataset(results_path)

        if VARIABLE_TO_PLOT not in ds:
            continue

        time = ds["time"].values
        values = ds[VARIABLE_TO_PLOT].values

        # Build label from metadata + config
        label = (
            f"{meta['git_branch']} | "
            f"{cfg.get('nc_file')[:7]}"
        )

        plt.plot(time[1:], values[1:], label=label)

    # Add Glint/BISICLES data for comparison
    plt.plot(icesheet_d["cs568"][0]["time"]-1850, bike_total_smb, label="Glint/BISICLES", color='black')

    plt.xlabel("Time")
    plt.ylabel(f"{VARIABLE_TO_PLOT} (Gt/a)")
    plt.title(f"{ICESHEET} {VARIABLE_TO_PLOT} Comparison")
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.savefig(f"{OUTPUT_DIR}/pyglint_test_TerraFIRMA_{cfg.get('suite_id')}_{VARIABLE_TO_PLOT}_{ICESHEET}.png", dpi=600)

if __name__ == "__main__":
    
    main()