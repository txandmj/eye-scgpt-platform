from fastapi import FastAPI, File, UploadFile, HTTPException, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
import os
import uuid
import shutil
import subprocess
import asyncio
from pathlib import Path
from typing import Optional
import json
import logging
from datetime import datetime
import threading
import tempfile
<<<<<<< HEAD
from app import auth_router
=======
from simple_umap_processor import run_two_stage_workflow
>>>>>>> 49c03b5 (add umap processor)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(
    title="Eye-scGPT Annotation Platform API",
    description="API for single-cell omics data annotation using scGPT",
    version="1.0.0"
)

app.include_router(auth_router.router)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, specify your frontend URL
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Create necessary directories
UPLOAD_DIR = Path("uploads")
RESULTS_DIR = Path("results")
MODELS_DIR = Path("models")
TEST_MODEL_DIR = MODELS_DIR / "test_model"

# Create directories if they don't exist
for directory in [UPLOAD_DIR, RESULTS_DIR, MODELS_DIR]:
    directory.mkdir(exist_ok=True)

# In-memory job tracking (in production, use Redis or database)
jobs = {}

def load_existing_jobs():
    """Load existing job data from files"""
    jobs_dir = Path("jobs")
    if not jobs_dir.exists():
        return
    
    for job_file in jobs_dir.glob("*.json"):
        try:
            with open(job_file, 'r') as f:
                job_data = json.load(f)
                job_id = job_data['id']
                jobs[job_id] = job_data
                logger.info(f"Loaded existing job: {job_id}")
        except Exception as e:
            logger.warning(f"Failed to load job file {job_file}: {e}")

# Load existing jobs on startup
load_existing_jobs()

# Job status enum
class JobStatus:
    PENDING = "pending"
    PROCESSING = "processing"
    COMPLETED = "completed"
    FAILED = "failed"

@app.get("/")
async def root():
    return {"message": "Eye-scGPT Annotation Platform API", "version": "1.0.0"}

@app.get("/health")
async def health_check():
    return {"status": "healthy", "timestamp": datetime.now().isoformat()}

@app.post("/api/upload")
async def upload_file(file: UploadFile = File(...)):
    """Upload a file and create a job for processing"""
    
    # Validate file type
    if not file.filename.lower().endswith('.h5ad'):
        raise HTTPException(status_code=400, detail="Only .h5ad files are supported")
    
    # Generate unique job ID
    job_id = str(uuid.uuid4())
    
    # Create job directory
    job_dir = UPLOAD_DIR / job_id
    job_dir.mkdir(exist_ok=True)
    
    # Save uploaded file
    file_path = job_dir / file.filename
    try:
        with open(file_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
    except Exception as e:
        logger.error(f"Error saving file: {e}")
        raise HTTPException(status_code=500, detail="Failed to save uploaded file")
    
    # Initialize job status
    jobs[job_id] = {
        "id": job_id,
        "status": JobStatus.PENDING,
        "filename": file.filename,
        "upload_time": datetime.now().isoformat(),
        "start_time": None,
        "end_time": None,
        "error": None,
        "results": None
    }
    
    logger.info(f"File uploaded successfully: {file.filename}, Job ID: {job_id}")
    
    return {
        "job_id": job_id,
        "status": JobStatus.PENDING,
        "message": "File uploaded successfully"
    }

@app.post("/api/annotate/{job_id}")
async def start_annotation(job_id: str, background_tasks: BackgroundTasks):
    """Start the annotation process for a specific job"""
    
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    if jobs[job_id]["status"] != JobStatus.PENDING:
        raise HTTPException(status_code=400, detail="Job is not in pending status")
    
    # Check if test model exists
    if not TEST_MODEL_DIR.exists():
        raise HTTPException(status_code=500, detail="Test model not found. Please ensure test_model directory exists in models/")
    
    # Start processing in background
    background_tasks.add_task(process_annotation, job_id)
    
    jobs[job_id]["status"] = JobStatus.PROCESSING
    jobs[job_id]["start_time"] = datetime.now().isoformat()
    
    logger.info(f"Started annotation process for job: {job_id}")
    
    return {
        "job_id": job_id,
        "status": JobStatus.PROCESSING,
        "message": "Annotation process started"
    }

async def process_annotation(job_id: str):
    """Process the annotation in background"""
    temp_dir = None
    try:
        job = jobs[job_id]
        job_dir = UPLOAD_DIR / job_id
        
        # Create temporary directory for processing
        temp_dir = Path(tempfile.mkdtemp(prefix=f"annotation_{job_id}_"))
        logger.info(f"Using temporary directory: {temp_dir}")
        
        # Get the uploaded file
        uploaded_files = list(job_dir.glob("*.h5ad"))
        if not uploaded_files:
            raise Exception("No .h5ad file found in job directory")
        
        uploaded_file = uploaded_files[0]
        
        # Step 1: Preprocessing
        logger.info(f"Starting step 1 (preprocessing) for job {job_id}")
        step1_dir = temp_dir / "step1"
        step1_dir.mkdir(exist_ok=True)
        
        # Set up environment for step1
        step1_env = os.environ.copy()
        step1_env.update({
            'DATASET_DIRECTORY': str(uploaded_file),
            'SAVE_DIR': str(step1_dir),
            'LOAD_MODEL': str(TEST_MODEL_DIR)
        })
        
        # Run step 1 directly
        result1 = subprocess.run(
            ["python", "step1_preprocess.py"],
            cwd=backend_dir,
            env=step1_env,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout
        )
        
        print("DEBUG STEP1 STDOUT:\n", result1.stdout)
        print("DEBUG STEP1 STDERR:\n", result1.stderr)
        
        if result1.returncode != 0:
            import traceback
            print("DEBUG STEP1 EXCEPTION:\n", traceback.format_exc())
            raise Exception(f"Step 1 failed: {result1.stderr}")
        
        logger.info(f"Step 1 completed for job {job_id}")
        
        # Step 2: Inference
        logger.info(f"Starting step 2 (inference) for job {job_id}")
        step2_dir = temp_dir / "step2"
        step2_dir.mkdir(exist_ok=True)
        
        # Set up environment for step2
        step2_env = os.environ.copy()
        step2_env.update({
            'LOAD_MODEL': str(TEST_MODEL_DIR),
            'TRAIN_ARGS': str(step1_dir / "train_args.yml"),
            'SAVE_DIR': str(step2_dir)
        })
        
        # Run step 2 directly
        result2 = subprocess.run(
            ["python", "step2_inference.py"],
            cwd=backend_dir,
            env=step2_env,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout
        )
        
        if result2.returncode != 0:
            raise Exception(f"Step 2 failed: {result2.stderr}")
        
        logger.info(f"Step 2 completed for job {job_id}")
        
        # Check if results exist
        temp_predictions_file = step2_dir / "predictions.csv"
        if not temp_predictions_file.exists():
            raise Exception("Predictions file not found after processing")
        
        # Create results directory and copy only the final results
        results_dir = RESULTS_DIR / job_id
        results_dir.mkdir(exist_ok=True)
        
        # Copy the final predictions file to results directory
        final_predictions_file = results_dir / "predictions.csv"
        shutil.copy2(temp_predictions_file, final_predictions_file)
        
        # Step 3: UMAP Processing
        logger.info(f"Starting UMAP processing for job {job_id}")
        umap_dir = results_dir / "umap"
        umap_dir.mkdir(exist_ok=True)
        
        try:
            # Run the simple two-stage UMAP workflow
            umap_results = run_two_stage_workflow(
                h5ad_file=str(uploaded_file),
                metadata_file=str(final_predictions_file),
                output_dir=str(umap_dir),
                base_name=uploaded_file.stem
            )
            
            logger.info(f"UMAP processing completed successfully for job {job_id}")
            logger.info(f"Generated {umap_results['num_plots']} UMAP plots")
            
            # Update job status with complete results
            jobs[job_id]["status"] = JobStatus.COMPLETED
            jobs[job_id]["end_time"] = datetime.now().isoformat()
            jobs[job_id]["results"] = {
                "predictions_file": str(final_predictions_file),
                "umap_results": umap_results,
                "enriched_h5ad_file": umap_results["enriched_h5ad_file"],
                "plot_files": umap_results["plot_files"],
                "num_plots": umap_results["num_plots"]
            }
            
            logger.info(f"Complete workflow finished successfully for job {job_id}")
            
        except Exception as umap_error:
            logger.warning(f"UMAP processing failed for job {job_id}: {str(umap_error)}")
            logger.warning("Continuing with annotation results only (no UMAP plots)")
            
            # Update job status with annotation results only
            jobs[job_id]["status"] = JobStatus.COMPLETED
            jobs[job_id]["end_time"] = datetime.now().isoformat()
            jobs[job_id]["results"] = {
                "predictions_file": str(final_predictions_file),
                "umap_error": str(umap_error)
            }
            
            logger.info(f"Annotation completed (without UMAP) for job {job_id}")
        
    except Exception as e:
        logger.error(f"Annotation failed for job {job_id}: {str(e)}")
        jobs[job_id]["status"] = JobStatus.FAILED
        jobs[job_id]["end_time"] = datetime.now().isoformat()
        jobs[job_id]["error"] = str(e)
    
    finally:
        # Clean up temporary directory
        if temp_dir and temp_dir.exists():
            try:
                shutil.rmtree(temp_dir)
                logger.info(f"Cleaned up temporary directory: {temp_dir}")
            except Exception as e:
                logger.warning(f"Failed to clean up temporary directory {temp_dir}: {e}")


@app.get("/api/job/{job_id}")
async def get_job_status(job_id: str):
    """Get the status of a specific job"""
    
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    return jobs[job_id]

@app.get("/api/jobs")
async def list_jobs():
    """List all jobs"""
    return {"jobs": list(jobs.values())}

@app.get("/api/download/{job_id}")
async def download_results(job_id: str):
    """Download the results of a completed job"""
    
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    if job["status"] != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job is not completed yet")
    
    results = job["results"]
    predictions_file = Path(results["predictions_file"])
    
    if not predictions_file.exists():
        raise HTTPException(status_code=404, detail="Results file not found")
    
    # Read the predictions file
    try:
        with open(predictions_file, "r") as f:
            content = f.read()
        
        response_data = {
            "job_id": job_id,
            "filename": predictions_file.name,
            "content": content,
            "download_url": f"/api/download/{job_id}/file",
            "umap_available": False
        }
        
        # Check if UMAP results are available
        if "umap_results" in results:
            response_data["umap_available"] = True
            response_data["num_plots"] = results.get("num_plots", 0)
            response_data["umap_download_url"] = f"/api/download/{job_id}/umap"
        
        return response_data
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error reading results: {str(e)}")

@app.get("/api/download/{job_id}/file")
async def download_file(job_id: str):
    """Download the actual results file"""
    
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    if job["status"] != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job is not completed yet")
    
    results = job["results"]
    predictions_file = Path(results["predictions_file"])
    
    if not predictions_file.exists():
        raise HTTPException(status_code=404, detail="Results file not found")
    
    from fastapi.responses import FileResponse
    return FileResponse(
        path=predictions_file,
        filename=predictions_file.name,
        media_type='text/csv'
    )

@app.get("/api/download/{job_id}/umap")
async def download_umap_results(job_id: str):
    """Download UMAP results as a zip file"""
    
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    if job["status"] != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job is not completed yet")
    
    results = job["results"]
    
    if "umap_results" not in results:
        raise HTTPException(status_code=404, detail="UMAP results not available")
    
    umap_results = results["umap_results"]
    
    # Create a zip file with all UMAP results
    import zipfile
    import io
    
    zip_buffer = io.BytesIO()
    
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        # Add enriched h5ad file
        enriched_h5ad_path = Path(umap_results["enriched_h5ad_file"])
        if enriched_h5ad_path.exists():
            zip_file.write(enriched_h5ad_path, enriched_h5ad_path.name)
        
        # Add plot files
        for plot_file in umap_results["plot_files"]:
            plot_path = Path(plot_file)
            if plot_path.exists():
                zip_file.write(plot_path, plot_path.name)
    
    zip_buffer.seek(0)
    
    from fastapi.responses import Response
    return Response(
        content=zip_buffer.getvalue(),
        media_type="application/zip",
        headers={"Content-Disposition": f"attachment; filename=umap_results_{job_id}.zip"}
    )


# Get the backend directory path
backend_dir = Path(__file__).parent

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
