#!/usr/bin/env python3
"""
Script to update existing completed jobs with UMAP results
"""

import os
import sys
from pathlib import Path
import json
from datetime import datetime

# Add the backend directory to the Python path
backend_dir = Path(__file__).parent
sys.path.insert(0, str(backend_dir))

from simple_umap_processor import run_two_stage_workflow

def find_completed_jobs():
    """Find all completed jobs that have predictions.csv but no UMAP results"""
    results_dir = Path("results")
    completed_jobs = []
    
    for job_dir in results_dir.iterdir():
        if not job_dir.is_dir():
            continue
            
        predictions_file = job_dir / "predictions.csv"
        if not predictions_file.exists():
            continue
            
        # Check if UMAP results exist
        umap_dir = job_dir / "umap"
        simple_umap_dir = job_dir / "simple_umap"
        
        has_umap = umap_dir.exists() and len(list(umap_dir.glob("*.png"))) > 0
        has_simple_umap = simple_umap_dir.exists() and len(list(simple_umap_dir.glob("*.png"))) > 0
        
        if has_umap or has_simple_umap:
            # Find the original h5ad file
            uploads_dir = Path("uploads")
            original_h5ad = None
            for upload_job_dir in uploads_dir.iterdir():
                if upload_job_dir.name == job_dir.name:
                    h5ad_files = list(upload_job_dir.glob("*.h5ad"))
                    if h5ad_files:
                        original_h5ad = h5ad_files[0]
                        break
            
            if original_h5ad and original_h5ad.exists():
                completed_jobs.append({
                    'job_id': job_dir.name,
                    'predictions_file': predictions_file,
                    'original_h5ad': original_h5ad,
                    'umap_dir': umap_dir if has_umap else simple_umap_dir,
                    'has_umap': has_umap,
                    'has_simple_umap': has_simple_umap
                })
    
    return completed_jobs

def update_job_with_umap_results(job_info):
    """Update a job with UMAP results information"""
    job_id = job_info['job_id']
    umap_dir = job_info['umap_dir']
    
    print(f"📊 Processing job {job_id}")
    print(f"   UMAP directory: {umap_dir}")
    
    # Count files in UMAP directory
    # Look for PNG files with the "umap" prefix pattern
    png_files = []
    for pattern in ["*.png", "umap*.png"]:
        png_files.extend(umap_dir.glob(pattern))
    png_files = list(set(png_files))  # Remove duplicates
    
    h5ad_files = list(umap_dir.glob("*.h5ad"))
    
    print(f"   PNG files: {len(png_files)}")
    print(f"   H5AD files: {len(h5ad_files)}")
    
    if len(png_files) == 0:
        print(f"   ⚠️  No PNG files found, skipping")
        return False
    
    # Create UMAP results structure
    enriched_h5ad_file = h5ad_files[0] if h5ad_files else None
    plot_files = [str(f) for f in png_files]
    
    umap_results = {
        "enriched_h5ad_file": str(enriched_h5ad_file) if enriched_h5ad_file else None,
        "plot_files": plot_files,
        "num_plots": len(png_files)
    }
    
    # Create job data structure
    job_data = {
        "id": job_id,
        "status": "completed",
        "filename": "EVAL_snRNA_no_enriched.h5ad",  # Default filename
        "upload_time": datetime.now().isoformat(),
        "start_time": datetime.now().isoformat(),
        "end_time": datetime.now().isoformat(),
        "error": None,
        "results": {
            "predictions_file": str(job_info['predictions_file']),
            "umap_results": umap_results,
            "enriched_h5ad_file": str(enriched_h5ad_file) if enriched_h5ad_file else None,
            "plot_files": plot_files,
            "num_plots": len(png_files)
        }
    }
    
    # Save job data to a file (for the backend to load)
    job_file = Path("jobs") / f"{job_id}.json"
    job_file.parent.mkdir(exist_ok=True)
    
    with open(job_file, 'w') as f:
        json.dump(job_data, f, indent=2)
    
    print(f"   ✅ Updated job data saved to {job_file}")
    return True

def main():
    """Main function to update existing jobs"""
    print("🔧 Updating Existing Jobs with UMAP Results")
    print("=" * 50)
    
    # Find completed jobs
    completed_jobs = find_completed_jobs()
    
    if not completed_jobs:
        print("❌ No completed jobs with UMAP results found")
        return
    
    print(f"📊 Found {len(completed_jobs)} completed jobs with UMAP results")
    
    # Process each job
    updated_count = 0
    for job_info in completed_jobs:
        try:
            if update_job_with_umap_results(job_info):
                updated_count += 1
        except Exception as e:
            print(f"   ❌ Failed to update job {job_info['job_id']}: {e}")
    
    print(f"\n🎉 Updated {updated_count} jobs successfully!")
    print("💡 The backend should now be able to serve UMAP results for these jobs")

if __name__ == "__main__":
    main()
