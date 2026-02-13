# Jet Flavor Tagging - Stage 2 Processing

Modern pipeline for converting event-level jet data to jet-level format for machine learning applications.

## Overview

This pipeline processes ROOT files from stage 1 (event-level) and converts them to stage 2 (jet-level) format:

**Input (Event-level):**
- Each event contains multiple jets
- Jet variables: `[[jet1_val, jet2_val], [jet1_val, jet2_val], ...]`
- PFCand variables: `[[[jet1_const1, jet1_const2, ...], [jet2_const1, ...]], ...]`

**Output (Jet-level):**
- Each row is one jet
- Jet variables: Single value per jet
- PFCand variables: Jagged array (variable-length list) per jet
- Flavor labels: Binary flags (`recojet_isB`, `recojet_isC`, etc.)

## Features

✅ **Modern Python** - Uses `uproot`/`awkward` instead of PyROOT  
✅ **YAML Configuration** - Clean, readable config instead of Python dicts  
✅ **Automatic NaN handling** - Converts NaN/None to 0.0 in constituent arrays  
✅ **Flexible flavor mapping** - Configure jet labels based on filename patterns  
✅ **Parallel processing** - Process multiple files and chunks in parallel  
✅ **Jagged arrays** - Maintains variable-length constituent arrays (no padding)

## Required Files

Download all these files to your working directory:

### Core Files (Required)
1. **`config.yaml`** - Configuration (variables, flavor mapping)
2. **`stage2.py`** - Main processing script
3. **`run_parallel.py`** - Parallel processing engine
4. **`steering.py`** - Easy-to-use steering file

## Quick Start

### 1. Install Dependencies

```bash
pip install uproot numpy akward
```

### 2. Configure

Edit `steering.py`:

```python
INPUT_DIR = "/path/to/your/stage1/files/"
OUTPUT_DIR = "/path/to/output/"  # note that if you use eos it may result in bad writes with high number of threads, so use working dir.
N_WORKERS = 64
```

### 3. Run

```bash
python steering.py
```

## Output

For input `Zbb.root`, get `Zbb_s2.root` with:
- 200 jets from 100 events (example)
- Flavor labels: `recojet_isB=1`, others=0
- All jet and pfcand variables

See full README for details on configuration, troubleshooting, and advanced usage.
