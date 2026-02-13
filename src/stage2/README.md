# Modern Jet Flavor Tagging Pipeline

A modernized pipeline for processing jet flavor tagging data, using `uproot` for efficient ROOT file I/O and YAML for configuration.

## Features

- **Modern Python**: Uses uproot/awkward for ROOT file handling (no PyROOT dependency)
- **YAML Configuration**: Clean, readable configuration instead of Python dictionaries
- **Automatic NaN handling**: Converts NaN and None values to 0.0 in constituent arrays
- **Flexible flavor mapping**: Configure jet flavor labels based on filename patterns
- **Parallel processing**: Process multiple files and chunks in parallel
- **Progress tracking**: Clear progress indicators and error reporting

## Files

- `config.yaml` - Configuration file defining variables and flavor mapping
- `steering.py` - **Main entry point** - Configure and run the pipeline here
- `stage2.py` - Core processing script (converts event-level to jet-level data)
- `run_parallel.py` - Parallel processing engine
- `example.py` - Example/test script demonstrating usage

## Installation

```bash
pip install uproot awkward pyyaml numpy
```

## Configuration

Edit `config.yaml` to customize:

### Flavor Mapping

Maps filename patterns to jet flavor labels:

```yaml
flavor_mapping:
  Zbb: recojet_isB  # Files containing "Zbb" → labeled as B-jets
  Zcc: recojet_isC  # Files containing "Zcc" → labeled as C-jets
  # ... etc
```

### Variables

Define which variables to include:

```yaml
jet_variables:
  - jet_p
  - jet_e
  # ...

pfcand_variables:
  - pfcand_erel_log
  - pfcand_btagSip3dSig
  # ...
```

## Usage

### Quick Start - Using Steering File (Recommended)

The easiest way to run the pipeline is using the steering file:

1. **Edit `steering.py`** - Set your input/output directories and parameters at the top
2. **Run it**: `python steering.py`

```python
# Example configuration in steering.py
INPUT_DIR = "/eos/experiment/fcc/ee/.../stage1/1.1.0/"
OUTPUT_DIR = "/eos/user/h/hfatehi/aleph-v1.1.0-dEdx/"
N_WORKERS = 64
N_FINAL_FILES = 1
```

Then simply run:
```bash
python steering.py
```

The steering file will:
- Validate your paths and configuration
- Show you what will be processed
- Ask for confirmation before starting
- Run the full parallel pipeline
- Report final status

### Single File Processing

Process a single ROOT file (or a chunk of it):

```bash
python stage2.py \
    input.root \
    output.root \
    0 \
    1000 \
    --config config.yaml
```

Arguments:
- `input.root` - Input ROOT file
- `output.root` - Output ROOT file
- `0` - Start entry
- `1000` - Stop entry (exclusive)
- `--config` - Path to config file (optional, default: config.yaml)

### Advanced - Direct Parallel Processing

For more control, you can call `run_parallel.py` directly:

```bash
python run_parallel.py \
    /path/to/input/dir \
    /path/to/output/dir \
    --workers 64 \
    --config config.yaml
```

Options:
- `--workers N` - Number of parallel workers (default: 64)
- `--config FILE` - Configuration file (default: config.yaml)
- `--stage2-script FILE` - Path to stage2.py (default: stage2.py)
- `--keep-parts` - Keep intermediate part files after merging
- `--n-final-files N` - Number of final merged files per sample (default: 1)
- `--chunk-size N` - Entries per chunk (default: auto-calculate)

## Data Structure

### Input (Event-level)

Each event contains multiple jets:
- Jet variables: `[jet1_value, jet2_value, jet3_value]`
- PFCand variables: `[[jet1_const1, jet1_const2, ...], [jet2_const1, ...], ...]`

### Output (Jet-level)

Each row is one jet:
- Jet variables: Single float per row
- PFCand variables: Jagged array (variable-length list) per row
- Flavor labels: One branch per flavor (`recojet_isB`, `recojet_isC`, etc.)

Example output structure:
```
jet_p: [45.2, 38.1, 52.3, ...]
jet_e: [46.5, 39.2, 53.1, ...]
recojet_isB: [1, 1, 1, ...]  # All 1s for Zbb file
recojet_isC: [0, 0, 0, ...]  # All 0s for Zbb file
pfcand_erel_log: [[−1.2, −0.8, ...], [−1.5, −0.9, ...], ...]
jet_npfcand: [15, 23, 18, ...]  # Number of constituents per jet
```

## NaN Handling

The pipeline automatically:
1. Converts `None`/`null` values to `0.0`
2. Converts `NaN` values to `0.0`
3. Applies this to all pfcand-level variables before writing

This ensures clean, numeric data for machine learning applications.

## Performance Notes

- **Chunk size**: Automatically calculated to balance parallelism and overhead
- **Memory**: Processing is done in chunks to avoid loading entire files
- **I/O**: uproot provides faster I/O than PyROOT for array operations
- **Merging**: Uses ROOT's `hadd` command for optimal merging

## Troubleshooting

### "Could not infer jet flavor from filename"

Make sure your input file name contains one of the patterns in `flavor_mapping` (e.g., "Zbb", "Zcc", etc.).

### "No ROOT files found"

Check that your input directory path is correct and contains `.root` files.

### Memory issues

Reduce the number of workers or set a smaller `--chunk-size`.

## Differences from Original Code

1. **No PyROOT**: Uses uproot/awkward instead
2. **Cleaner config**: YAML instead of Python dictionaries
3. **Better error handling**: Meaningful error messages and validation
4. **NaN safety**: Automatic sanitization of problematic values
5. **Type hints**: Modern Python typing for better IDE support
6. **Progress tracking**: Clear indication of processing status
7. **Flexible chunking**: Can specify chunk size or let it auto-calculate

## License

Same as original codebase.
