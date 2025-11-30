#!/usr/bin/env python3
"""
Simple UMAP Processor - Two Stage Process
Stage 1: Data Preparation (addmetadata) - Merge CSV metadata into h5ad
Stage 2: Plot Generation (umapby) - Generate UMAP plots using existing coordinates
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import seaborn as sns
from pathlib import Path
import logging
import shutil

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def h5ad2addmetadata(indata, metadata):
    """
    Stage 1: Data Preparation - Merge metadata into AnnData object
    
    Args:
        indata: AnnData object from h5ad file
        metadata: pandas DataFrame with new metadata
    
    Returns:
        AnnData object with merged metadata
    """
    # Find common cells between h5ad and metadata
    common = indata.obs.index.intersection(metadata.index)
    
    if len(common) < len(indata.obs.index):
        logger.warning(f"Some cells are discarded due to not in metadata {len(common)} < {len(indata.obs.index)}.")
    
    if len(common) == 0:
        logger.warning("Empty cell barcode in metadata.")
        return indata
    
    # Filter to common cells and add metadata
    indata = indata[common].copy()
    metadata = metadata.loc[common]
    
    # Copy all columns from metadata to obs
    header = metadata.columns.tolist()
    for h in header:
        indata.obs[h] = metadata[h].tolist()
    
    return indata

def stage1_addmetadata(h5ad_file, metadata_file, output_file=None):
    """
    Stage 1: Data Preparation - Add metadata to h5ad file
    
    Args:
        h5ad_file: Path to input .h5ad file
        metadata_file: Path to metadata CSV file
        output_file: Path to output .h5ad file (optional)
    
    Returns:
        Path to the output file
    """
    logger.info("=== Stage 1: Data Preparation (addmetadata) ===")
    logger.info(f"Loading h5ad file: {h5ad_file}")
    indata = sc.read_h5ad(h5ad_file)
    
    logger.info(f"Loading metadata file: {metadata_file}")
    metadata = pd.read_csv(metadata_file, header=0, index_col=0)
    
    logger.info("Merging metadata into h5ad data")
    indata = h5ad2addmetadata(indata, metadata)
    
    if output_file is None:
        output_file = h5ad_file.replace('.h5ad', '_with_metadata.h5ad')
    
    logger.info(f"Saving enriched h5ad file: {output_file}")
    sc.write(filename=output_file, adata=indata)
    
    return output_file

def stage2_umapby(h5ad_file, output_dir, base_name=None, width=5, height=5):
    """
    Stage 2: Plot Generation - Generate UMAP plots using existing coordinates
    
    Args:
        h5ad_file: Path to enriched .h5ad file from Stage 1
        output_dir: Directory to save UMAP plots
        base_name: Base name for output files (optional)
        width: Figure width
        height: Figure height
    
    Returns:
        List of generated plot files
    """
    logger.info("=== Stage 2: Plot Generation (umapby) ===")
    logger.info(f"Loading enriched h5ad file: {h5ad_file}")
    x = sc.read_h5ad(h5ad_file)
    
    # Critical check: UMAP coordinates must already exist
    if 'X_umap' not in x.obsm:
        logger.error('X_umap is missing. UMAP coordinates not found in the data.')
        raise ValueError('X_umap is missing. See scanpy.tl.umap()')
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if base_name is None:
        base_name = Path(h5ad_file).stem
    
    # Set scanpy figure parameters
    sc.set_figure_params(dpi_save=500, figsize=(width, height))
    
    generated_files = []
    
    # Iterate through every column in x.obs
    logger.info(f"Generating plots for {len(x.obs.columns)} metadata columns")
    
    for splitby in x.obs.columns:
        ncolor = len(x.obs[splitby].value_counts())
        
        if ncolor < 100:
            # Categorical variable - Generate 3 plots
            logger.info(f"Generating categorical plots for {splitby} ({ncolor} categories)")
            
            # Plot 1: Standard UMAP with legend on side
            plot1 = f"{base_name}_umap_{splitby}_wolabel.png"
            sc.pl.umap(x, color=splitby, frameon=False, show=False, title=None, save=plot1)
            
            # Plot 2: UMAP with legend labels on data points
            plot2 = f"{base_name}_umap_{splitby}_ondata.png"
            sc.pl.umap(x, color=splitby, frameon=False, show=False, title=None, save=plot2,
                      legend_loc='on data', legend_fontsize='xx-small', legend_fontweight='normal')
            
            # Plot 3: Same as Plot 2 but with text outline for better visibility
            plot3 = f"{base_name}_umap_{splitby}_fontline.png"
            sc.pl.umap(x, color=splitby, frameon=False, show=False, title=None, save=plot3,
                      legend_loc='on data', legend_fontsize='xx-small', legend_fontweight='normal',
                      legend_fontoutline=1)
            
            generated_files.extend([output_dir / f for f in [plot1, plot2, plot3]])
            
        else:
            # Continuous variable - Generate 2 plots
            logger.info(f"Generating continuous plots for {splitby} ({ncolor} values)")
            
            # Use custom palette for many categories
            palette = sns.husl_palette(ncolor)
            
            # Plot 1: UMAP with continuous color gradient
            plot1 = f"{base_name}_umap_{splitby}_wolabel.png"
            sc.pl.umap(x, color=splitby, palette=palette, frameon=False, show=False, title=None, save=plot1)
            
            # Plot 2: Same as Plot 1 (continuous variables don't need on-data labels)
            plot2 = f"{base_name}_umap_{splitby}_ondata.png"
            sc.pl.umap(x, color=splitby, palette=palette, frameon=False, show=False, title=None, save=plot2)
            
            generated_files.extend([output_dir / f for f in [plot1, plot2]])
        
        logger.info(f"Generated plots for metadata column: {splitby}")
    
    # Move plots from figures/ directory to the correct output directory
    figures_dir = Path("figures")
    if figures_dir.exists():
        logger.info("Moving plots from figures/ directory to output directory...")
        for plot_file in figures_dir.glob(f"umap{base_name}_umap_*.png"):
            target_file = output_dir / plot_file.name
            # Use shutil.move() which handles cross-filesystem moves correctly
            # It automatically copies then removes if rename fails (cross-device link error)
            shutil.move(str(plot_file), str(target_file))
            logger.info(f"Moved {plot_file.name} to {target_file}")
        
        # After moving, get the actual moved files
        generated_files = list(output_dir.glob(f"umap{base_name}_umap_*.png"))
    
    logger.info(f"Plot generation completed. Total plots generated: {len(generated_files)}")
    return generated_files

def run_two_stage_workflow(h5ad_file, metadata_file, output_dir, base_name=None):
    """
    Complete two-stage UMAP workflow
    
    Args:
        h5ad_file: Path to input .h5ad file
        metadata_file: Path to metadata CSV file
        output_dir: Directory to save results
        base_name: Base name for output files
    
    Returns:
        Dictionary with paths to generated files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if base_name is None:
        base_name = Path(h5ad_file).stem
    
    # Stage 1: Data Preparation
    enriched_h5ad_file = output_dir / f"{base_name}_enriched.h5ad"
    stage1_addmetadata(h5ad_file, metadata_file, str(enriched_h5ad_file))
    
    # Stage 2: Plot Generation
    plot_files = stage2_umapby(str(enriched_h5ad_file), output_dir, base_name)
    
    return {
        "enriched_h5ad_file": str(enriched_h5ad_file),
        "plot_files": [str(f) for f in plot_files],
        "num_plots": len(plot_files)
    }

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 4:
        print("Usage: python simple_umap_processor.py <h5ad_file> <metadata_file> <output_dir> [base_name]")
        print("")
        print("Two-stage process:")
        print("Stage 1: Merge CSV metadata into h5ad file")
        print("Stage 2: Generate UMAP plots using existing coordinates")
        sys.exit(1)
    
    h5ad_file = sys.argv[1]
    metadata_file = sys.argv[2]
    output_dir = sys.argv[3]
    base_name = sys.argv[4] if len(sys.argv) > 4 else None
    
    try:
        results = run_two_stage_workflow(h5ad_file, metadata_file, output_dir, base_name)
        print(f"\n✅ Two-stage UMAP workflow completed successfully!")
        print(f"📊 Generated {results['num_plots']} plots")
        print(f"📁 Enriched h5ad file: {results['enriched_h5ad_file']}")
        print(f"🖼️  Plot files saved to: {output_dir}")
    except Exception as e:
        logger.error(f"Two-stage UMAP workflow failed: {str(e)}")
        sys.exit(1)
