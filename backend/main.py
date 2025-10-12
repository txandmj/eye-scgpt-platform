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
from app import auth_router

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
        
        if result1.returncode != 0:
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
        
        # Update job status
        jobs[job_id]["status"] = JobStatus.COMPLETED
        jobs[job_id]["end_time"] = datetime.now().isoformat()
        jobs[job_id]["results"] = {
            "predictions_file": str(final_predictions_file)
        }
        
        logger.info(f"Annotation completed successfully for job {job_id}")
        
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
        
        return {
            "job_id": job_id,
            "filename": predictions_file.name,
            "content": content,
            "download_url": f"/api/download/{job_id}/file"
        }
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

# Get the backend directory path
backend_dir = Path(__file__).parent

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
