#!/usr/bin/env python3
"""
Test script for the simple two-stage UMAP processor
"""

import sys
from pathlib import Path
from simple_umap_processor import stage1_addmetadata, stage2_umapby, run_two_stage_workflow

def test_with_latest_results():
    """Test with the most recent results"""
    
    # Find the most recent results
    results_dir = Path("results")
    if not results_dir.exists():
        print("❌ No results directory found!")
        return False
    
    # Find the most recent job directory
    job_dirs = [d for d in results_dir.iterdir() if d.is_dir()]
    if not job_dirs:
        print("❌ No job directories found!")
        return False
    
    latest_job_dir = max(job_dirs, key=lambda x: x.stat().st_mtime)
    predictions_file = latest_job_dir / "predictions.csv"
    
    if not predictions_file.exists():
        print(f"❌ No predictions.csv found in {latest_job_dir}")
        return False
    
    # Find the original h5ad file
    uploads_dir = Path("uploads")
    original_h5ad = None
    if uploads_dir.exists():
        for job_dir in uploads_dir.iterdir():
            if job_dir.name == latest_job_dir.name:
                h5ad_files = list(job_dir.glob("*.h5ad"))
                if h5ad_files:
                    original_h5ad = h5ad_files[0]
                    break
    
    if not original_h5ad or not original_h5ad.exists():
        print("❌ Original h5ad file not found!")
        return False
    
    # Create output directory
    output_dir = latest_job_dir / "umap"
    output_dir.mkdir(exist_ok=True)
    
    print("🧪 Testing Simple Two-Stage UMAP Processor")
    print("=" * 50)
    print(f"📁 Job ID: {latest_job_dir.name}")
    print(f"📄 Input h5ad: {original_h5ad}")
    print(f"📊 Predictions: {predictions_file}")
    print(f"📂 Output directory: {output_dir}")
    print()
    
    try:
        # Run the complete workflow
        results = run_two_stage_workflow(
            h5ad_file=str(original_h5ad),
            metadata_file=str(predictions_file),
            output_dir=str(output_dir),
            base_name=original_h5ad.stem
        )
        
        print("✅ Two-stage UMAP workflow completed successfully!")
        print(f"📊 Generated {results['num_plots']} plots")
        print(f"📁 Enriched h5ad file: {results['enriched_h5ad_file']}")
        print(f"🖼️  Plot files: {len(results['plot_files'])}")
        
        return True
        
    except Exception as e:
        print(f"❌ Two-stage UMAP workflow failed: {str(e)}")
        return False

def test_individual_stages():
    """Test individual stages separately"""
    
    # Find the most recent results
    results_dir = Path("results")
    job_dirs = [d for d in results_dir.iterdir() if d.is_dir()]
    latest_job_dir = max(job_dirs, key=lambda x: x.stat().st_mtime)
    
    predictions_file = latest_job_dir / "predictions.csv"
    
    # Find the original h5ad file
    uploads_dir = Path("uploads")
    original_h5ad = None
    for job_dir in uploads_dir.iterdir():
        if job_dir.name == latest_job_dir.name:
            h5ad_files = list(job_dir.glob("*.h5ad"))
            if h5ad_files:
                original_h5ad = h5ad_files[0]
                break
    
    # Create output directory
    output_dir = latest_job_dir / "umap"
    output_dir.mkdir(exist_ok=True)
    
    print("🧪 Testing Individual Stages")
    print("=" * 50)
    
    try:
        # Stage 1: Data Preparation
        print("Stage 1: Data Preparation (addmetadata)")
        enriched_h5ad = output_dir / f"{original_h5ad.stem}_enriched.h5ad"
        stage1_addmetadata(str(original_h5ad), str(predictions_file), str(enriched_h5ad))
        print(f"✅ Stage 1 completed: {enriched_h5ad}")
        
        # Stage 2: Plot Generation
        print("\nStage 2: Plot Generation (umapby)")
        plot_files = stage2_umapby(str(enriched_h5ad), str(output_dir), original_h5ad.stem)
        print(f"✅ Stage 2 completed: {len(plot_files)} plots generated")
        
        return True
        
    except Exception as e:
        print(f"❌ Individual stages test failed: {str(e)}")
        return False

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "individual":
        success = test_individual_stages()
    else:
        success = test_with_latest_results()
    
    if not success:
        sys.exit(1)
