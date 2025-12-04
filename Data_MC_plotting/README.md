# Data & MC plotting for ALEPH 

Flexible script to make plots with data and MC for the revived ALEPH datasets, all plotting details are set via config files.

## Setup and input preparation

Follow the instructions of the main repo to setup the key4hep stack. 

Common input files can be found on eos in `/eos/user/h/hfatehi/D0fliped-good/` currently. Otherwise, use the scripts in this repo to produce ` stage 1` or `inference` type of input files. 

## Run 

With a given `plotting_config.py`, run the following command:

```
python plot_data_mc.py -c plotting_config
```

## Existing config files

### Stage 1 

The config file `plotting_config_stage1.py` is used to plot variables from the `stage1` ntuples that are used in the tagger. 

### Inference 

The config file `plotting_config_inference.py` is used to plot variables from the output ntuples after running the tagger inference. 

## Config variables

The following variables are set in the config files:

- `inputs_path`: Path to the input ntuples.  
- `outputs_path` : Path where to store the plots. 
- `plots_dict` : This is a dictionary of named tuples specifying the exact plot settings. It is loaded from the helper file `Zqq_plots.py`. Add new dictionaries there as needed. See below for explanation. 
- `data` and `mc_processes`: These define the filenames/samples for data and MC processes plotted. Again they use named dictionaries and are defined in the helper file `Zqq_processes.py`. See below for explanation. 
- `year`: The data taking year, will be added as label on the plots.
- `sel_tag`: A string that specifies the analysis or selection level of events, will be added as label on the plots.
- `lumi` : The luminosity to normalise MC to,  will be added as label on the plots.
- `ecm` : The center of mass energy,  will be added as label on the plots.


### PlotSpecs config

### ProcessSpecs config

