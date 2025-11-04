# Simple Two-Stage UMAP Processor

This is a simplified UMAP processor that follows exactly the two-stage process you described.

## Overview

The processor implements your exact workflow:

1. **Stage 1: Data Preparation (addmetadata)** - Merge CSV metadata into h5ad file
2. **Stage 2: Plot Generation (umapby)** - Generate UMAP plots using existing coordinates

## Files

- `simple_umap_processor.py` - Main processor with two-stage workflow
- `test_simple_umap.py` - Test script to run the processor

## Usage

### Command Line

```bash
python simple_umap_processor.py <h5ad_file> <metadata_file> <output_dir> [base_name]
```

### Example

```bash
python simple_umap_processor.py \
  EVAL_snRNA_no_enriched.h5ad \
  predictions.csv \
  output_directory \
  my_analysis
```

### Python API

```python
from simple_umap_processor import run_two_stage_workflow

results = run_two_stage_workflow(
    h5ad_file="input.h5ad",
    metadata_file="predictions.csv", 
    output_dir="output",
    base_name="analysis"
)
```

## Two-Stage Process

### Stage 1: Data Preparation (addmetadata)

**Purpose**: Merge new cell metadata from CSV into the h5ad file

**What it does**:
1. Loads the main .h5ad file (expected to already contain UMAP coordinates)
2. Loads the predictions.csv file with new metadata indexed by cell barcodes
3. Finds cells common to both files
4. Copies all columns from CSV into the .obs attribute of the h5ad file
5. Saves enriched h5ad file

**Result**: h5ad file with merged metadata (original + new predictions)

### Stage 2: Plot Generation (umapby)

**Purpose**: Generate UMAP plots using existing coordinates

**What it does**:
1. Loads the enriched h5ad file from Stage 1
2. **Critical check**: Verifies 'X_umap' exists in x.obsm
3. Iterates through every column in x.obs
4. For each metadata column, generates appropriate plots:

**Categorical variables** (< 100 unique values):
- `_wolabel.png` - Standard UMAP with legend on side
- `_ondata.png` - UMAP with legend labels on data points  
- `_fontline.png` - Same as _ondata.png but with text outline

**Continuous variables** (≥ 100 unique values):
- `_wolabel.png` - UMAP with continuous color gradient
- `_ondata.png` - Same as _wolabel.png (duplicate for consistency)

**Result**: Folder full of PNG images for every metadata column

## Requirements

- The input h5ad file must already contain UMAP coordinates (`X_umap` in `.obsm`)
- CSV file must have cell barcodes as index matching the h5ad file
- scanpy, pandas, seaborn, matplotlib installed

## Output

### Files Generated

1. **Enriched h5ad file**: `{base_name}_enriched.h5ad`
   - Original data + merged metadata from CSV

2. **UMAP plots**: Multiple PNG files per metadata column
   - Categorical: 3 plots per column (wolabel, ondata, fontline)
   - Continuous: 2 plots per column (wolabel, ondata)

### Example Output Structure

```
output_directory/
├── analysis_enriched.h5ad          # Enriched data file
├── analysis_umap_celltype_wolabel.png
├── analysis_umap_celltype_ondata.png  
├── analysis_umap_celltype_fontline.png
├── analysis_umap_predictions_wolabel.png
├── analysis_umap_predictions_ondata.png
├── analysis_umap_predictions_fontline.png
└── ... (more plots for each metadata column)
```

## Testing

Run the test script to verify everything works:

```bash
python test_simple_umap.py
```

This will automatically find your latest results and run the two-stage process.

## Key Features

- **Simple two-stage process** exactly as you described
- **Automatic detection** of categorical vs continuous variables
- **Proper error handling** for missing UMAP coordinates
- **Clear logging** showing progress through each stage
- **Flexible output** with customizable base names and directories

This processor implements your exact workflow and should work seamlessly with your existing data pipeline!
