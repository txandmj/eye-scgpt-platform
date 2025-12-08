from fastapi import FastAPI, File, UploadFile, HTTPException, BackgroundTasks, Depends, Header, Form
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, FileResponse, Response
import os
import uuid
import shutil
import subprocess
import asyncio
from pathlib import Path
from typing import Optional, Dict, Any
import json
import logging
from datetime import datetime
import threading
import tempfile
from pydantic import BaseModel
from simple_umap_processor import run_two_stage_workflow
from app import auth_router
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(
    title="Eye-scGPT Annotation Platform API",
    description="API for single-cell omics data annotation using scGPT",
    version="1.0.0"
)

# Include auth router
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
JOBS_DIR = Path("jobs")
for directory in [UPLOAD_DIR, RESULTS_DIR, MODELS_DIR, JOBS_DIR]:
    directory.mkdir(exist_ok=True)

# Get the backend directory path (needed for subprocess calls)
backend_dir = Path(__file__).parent

# In-memory job tracking (in production, use Redis or database)
jobs = {}

def load_existing_jobs():
    """Load existing job data from files"""
    if not JOBS_DIR.exists():
        return
    
    for job_file in JOBS_DIR.glob("*.json"):
        try:
            with open(job_file, 'r') as f:
                job_data = json.load(f)
                job_id = job_data['id']
                jobs[job_id] = job_data
                logger.info(f"Loaded existing job: {job_id}")
        except Exception as e:
            logger.warning(f"Failed to load job file {job_file}: {e}")

def save_job_to_file(job_id: str, job_data: dict):
    """Save job data to file for persistence"""
    try:
        job_file = JOBS_DIR / f"{job_id}.json"
        with open(job_file, 'w') as f:
            json.dump(job_data, f, indent=2)
    except Exception as e:
        logger.warning(f"Failed to save job file for {job_id}: {e}")

def send_completion_email(job_id: str, user_email: str, user_name: Optional[str], filename: str, frontend_url: str = None):
    """Send email notification when job completes"""
    # Get SMTP configuration from environment variables
    smtp_host = os.getenv("SMTP_HOST")
    smtp_port = int(os.getenv("SMTP_PORT", "587"))
    smtp_user = os.getenv("SMTP_USER")
    smtp_password = os.getenv("SMTP_PASSWORD")
    smtp_from_email = os.getenv("SMTP_FROM_EMAIL", smtp_user)
    frontend_url = frontend_url or os.getenv("FRONTEND_URL", "http://localhost:3000")
    
    # Skip email if SMTP is not configured
    if not smtp_host or not smtp_user or not smtp_password:
        logger.info(f"Email notification skipped for job {job_id}: SMTP not configured")
        return
    
    if not user_email:
        logger.warning(f"Email notification skipped for job {job_id}: No user email")
        return
    
    try:
        # Create message
        msg = MIMEMultipart('alternative')
        msg['Subject'] = f"Your analysis job is complete: {filename}"
        msg['From'] = smtp_from_email
        msg['To'] = user_email
        
        # Create email body
        user_display_name = user_name or "User"
        job_link = f"{frontend_url}/#/history"  # Link to history page
        download_link = f"{frontend_url}/#/download?jobId={job_id}"  # Link to download page
        
        text_body = f"""
Hello {user_display_name},

Your analysis job has completed successfully!

Job Details:
- Job ID: {job_id}
- Filename: {filename}
- Status: Completed

You can view and download your results at:
{download_link}

Or check your job history at:
{job_link}

Thank you for using Eye-scGPT Annotation Platform!
"""
        
        html_body = f"""
<!DOCTYPE html>
<html>
<head>
    <style>
        body {{ font-family: Arial, sans-serif; line-height: 1.6; color: #333; }}
        .container {{ max-width: 600px; margin: 0 auto; padding: 20px; }}
        .header {{ background-color: #667eea; color: white; padding: 20px; text-align: center; border-radius: 5px 5px 0 0; }}
        .content {{ background-color: #f9f9f9; padding: 30px; border-radius: 0 0 5px 5px; }}
        .button {{ display: inline-block; padding: 12px 24px; background-color: #667eea; color: white; text-decoration: none; border-radius: 5px; margin: 10px 5px; }}
        .button:hover {{ background-color: #5568d3; }}
        .job-details {{ background-color: white; padding: 15px; border-radius: 5px; margin: 20px 0; }}
        .job-detail {{ margin: 10px 0; }}
        .footer {{ text-align: center; color: #666; font-size: 12px; margin-top: 30px; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>✅ Analysis Complete!</h1>
        </div>
        <div class="content">
            <p>Hello {user_display_name},</p>
            <p>Your analysis job has completed successfully!</p>
            
            <div class="job-details">
                <div class="job-detail"><strong>Job ID:</strong> {job_id}</div>
                <div class="job-detail"><strong>Filename:</strong> {filename}</div>
                <div class="job-detail"><strong>Status:</strong> ✅ Completed</div>
            </div>
            
            <p>You can view and download your results:</p>
            <a href="{download_link}" class="button">📥 View Results</a>
            <a href="{job_link}" class="button">📊 Job History</a>
            
            <div class="footer">
                <p>Thank you for using Eye-scGPT Annotation Platform!</p>
            </div>
        </div>
    </div>
</body>
</html>
"""
        
        # Add parts
        part1 = MIMEText(text_body, 'plain')
        part2 = MIMEText(html_body, 'html')
        msg.attach(part1)
        msg.attach(part2)
        
        # Send email
        with smtplib.SMTP(smtp_host, smtp_port) as server:
            server.starttls()
            server.login(smtp_user, smtp_password)
            server.send_message(msg)
        
        logger.info(f"Email notification sent for job {job_id} to {user_email}")
        
    except Exception as e:
        logger.error(f"Failed to send email notification for job {job_id}: {str(e)}")
        # Don't raise exception - email failure shouldn't break the job completion

# Load existing jobs on startup
load_existing_jobs()

# Job status enum
class JobStatus:
    UPLOADING = "uploading"  # Job created, file is being uploaded
    PENDING = "pending"
    PROCESSING = "processing"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"  # Job was cancelled by user


# Authentication dependency - simplified version for header-based auth
async def require_auth(
    x_google_sub: Optional[str] = Header(default=None, alias="X-Google-Sub"),
    x_user_email: Optional[str] = Header(default=None, alias="X-User-Email"),
    x_user_name: Optional[str] = Header(default=None, alias="X-User-Name"),
) -> Dict[str, Any]:
    """Require authentication - returns user info or raises exception"""
    if not x_google_sub:
        raise HTTPException(
            status_code=401,
            detail="Authentication required. Please sign in with Google."
        )
    return {
        "google_sub": x_google_sub,
        "email": x_user_email,
        "name": x_user_name,
        "authenticated": True
    }

@app.get("/")
async def root():
    return {"message": "Eye-scGPT Annotation Platform API", "version": "1.0.0"}

@app.get("/health")
async def health_check():
    return {"status": "healthy", "timestamp": datetime.now().isoformat()}

class PrepareUploadRequest(BaseModel):
    filename: str

@app.post("/api/upload/prepare")
async def prepare_upload(
    request: PrepareUploadRequest,
    user: Dict[str, Any] = Depends(require_auth)
):
    """Create a job before upload starts - returns job_id immediately"""
    filename = request.filename
    # Generate unique job ID
    job_id = str(uuid.uuid4())
    
    # Create job directory
    job_dir = UPLOAD_DIR / job_id
    job_dir.mkdir(exist_ok=True)
    
    # Initialize job status as UPLOADING
    jobs[job_id] = {
        "id": job_id,
        "status": JobStatus.UPLOADING,
        "filename": filename,
        "upload_time": datetime.now().isoformat(),
        "start_time": None,
        "end_time": None,
        "error": None,
        "results": None,
        "user_id": user.get("google_sub"),
        "user_email": user.get("email"),
        "user_name": user.get("name")
    }
    
    # Save job to file immediately
    save_job_to_file(job_id, jobs[job_id])
    
    logger.info(f"Prepared job {job_id} for file {filename}")
    
    return {
        "job_id": job_id,
        "status": JobStatus.UPLOADING,
        "message": "Job created, ready for upload"
    }

@app.post("/api/upload")
async def upload_file(
    file: UploadFile = File(...),
    job_id: Optional[str] = Form(None),
    user: Dict[str, Any] = Depends(require_auth)
):
    """Upload a file for a prepared job or create new job (requires authentication)"""
    
    # Validate file type
    if not file.filename.lower().endswith('.h5ad'):
        raise HTTPException(status_code=400, detail="Only .h5ad files are supported")
    
    # If job_id provided, use existing job (from prepare_upload)
    # Otherwise create new job (backward compatibility)
    if job_id and job_id in jobs:
        # Verify user owns this job
        job_user_id = jobs[job_id].get("user_id")
        if job_user_id and job_user_id != user.get("google_sub"):
            raise HTTPException(status_code=403, detail="You don't have permission to access this job")
        
        # Verify job is in UPLOADING status
        if jobs[job_id]["status"] != JobStatus.UPLOADING:
            raise HTTPException(status_code=400, detail="Job is not in uploading status")
        
        job_dir = UPLOAD_DIR / job_id
    else:
        # Create new job (backward compatibility)
        job_id = str(uuid.uuid4())
        job_dir = UPLOAD_DIR / job_id
        job_dir.mkdir(exist_ok=True)
        
        jobs[job_id] = {
            "id": job_id,
            "status": JobStatus.UPLOADING,
            "filename": file.filename,
            "upload_time": datetime.now().isoformat(),
            "start_time": None,
            "end_time": None,
            "error": None,
            "results": None,
            "user_id": user.get("google_sub"),
            "user_email": user.get("email"),
            "user_name": user.get("name")
        }
        save_job_to_file(job_id, jobs[job_id])
    
    logger.info(f"Starting upload for job {job_id}, file: {file.filename}")
    
    # Save uploaded file
    file_path = job_dir / file.filename
    try:
        with open(file_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
        
        # Update job status to PENDING after successful upload
        jobs[job_id]["status"] = JobStatus.PENDING
        save_job_to_file(job_id, jobs[job_id])
        
        logger.info(f"File uploaded successfully: {file.filename}, Job ID: {job_id}")
        
        return {
            "job_id": job_id,
            "status": JobStatus.PENDING,
            "message": "File uploaded successfully"
        }
    except Exception as e:
        # Mark job as failed if upload fails
        jobs[job_id]["status"] = JobStatus.FAILED
        jobs[job_id]["error"] = f"Upload failed: {str(e)}"
        jobs[job_id]["end_time"] = datetime.now().isoformat()
        save_job_to_file(job_id, jobs[job_id])
        
        logger.error(f"Error saving file: {e}")
        raise HTTPException(status_code=500, detail="Failed to save uploaded file")

@app.post("/api/jobs/{job_id}/cancel")
async def cancel_job(
    job_id: str,
    user: Dict[str, Any] = Depends(require_auth)
):
    """Cancel a job (requires authentication)"""
    
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    # Verify user owns this job
    job_user_id = job.get("user_id")
    if job_user_id and job_user_id != user.get("google_sub"):
        raise HTTPException(status_code=403, detail="You don't have permission to access this job")
    
    # Check if job can be cancelled
    current_status = job["status"]
    if current_status in [JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELLED]:
        raise HTTPException(
            status_code=400, 
            detail=f"Job cannot be cancelled. Current status: {current_status}"
        )
    
    # Mark job as cancelled
    jobs[job_id]["status"] = JobStatus.CANCELLED
    jobs[job_id]["end_time"] = datetime.now().isoformat()
    jobs[job_id]["error"] = "Job cancelled by user"
    
    # Save updated job status
    save_job_to_file(job_id, jobs[job_id])
    
    logger.info(f"Job {job_id} cancelled by user")
    
    return {
        "job_id": job_id,
        "status": JobStatus.CANCELLED,
        "message": "Job cancelled successfully"
    }

@app.post("/api/annotate/{job_id}")
async def start_annotation(
    job_id: str,
    background_tasks: BackgroundTasks,
    user: Dict[str, Any] = Depends(require_auth)
):
    """Start the annotation process for a specific job (requires authentication)"""
    
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    # Verify user owns this job
    job_user_id = jobs[job_id].get("user_id")
    if job_user_id and job_user_id != user.get("google_sub"):
        raise HTTPException(status_code=403, detail="You don't have permission to access this job")
    
    if jobs[job_id]["status"] != JobStatus.PENDING:
        raise HTTPException(status_code=400, detail="Job is not in pending status")
    
    # Check if test model exists
    if not TEST_MODEL_DIR.exists():
        raise HTTPException(status_code=500, detail="Test model not found. Please ensure test_model directory exists in models/")
    
    # Start processing in background
    background_tasks.add_task(process_annotation, job_id)
    
    jobs[job_id]["status"] = JobStatus.PROCESSING
    jobs[job_id]["start_time"] = datetime.now().isoformat()
    
    # Save updated job status
    save_job_to_file(job_id, jobs[job_id])
    
    logger.info(f"Started annotation process for job: {job_id}")
    
    return {
        "job_id": job_id,
        "status": JobStatus.PROCESSING,
        "message": "Annotation process started"
    }

def _run_subprocess_blocking(cmd, cwd, env, timeout):
    """Blocking subprocess runner - to be executed in thread pool"""
    return subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        timeout=timeout
    )

async def process_annotation(job_id: str):
    """Process the annotation in background - runs blocking subprocess calls in thread pool"""
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
        
        # Run step 1 in thread pool to avoid blocking event loop
        loop = asyncio.get_event_loop()
        result1 = await loop.run_in_executor(
            None,
            _run_subprocess_blocking,
            ["python", "step1_preprocess.py"],
            backend_dir,
            step1_env,
            3600  # 1 hour timeout
        )
        
        if result1.returncode != 0:
            raise Exception(f"Step 1 failed: {result1.stderr}")
        
        logger.info(f"Step 1 completed for job {job_id}")
        
        # Check if job was cancelled after step 1
        if job_id not in jobs or jobs[job_id]["status"] == JobStatus.CANCELLED:
            logger.info(f"Job {job_id} was cancelled after step 1, aborting")
            return
        
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
        
        # Run step 2 in thread pool to avoid blocking event loop
        result2 = await loop.run_in_executor(
            None,
            _run_subprocess_blocking,
            ["python", "step2_inference.py"],
            backend_dir,
            step2_env,
            3600  # 1 hour timeout
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
            # Run the simple two-stage UMAP workflow (this is CPU-bound, run in thread pool)
            loop = asyncio.get_event_loop()
            umap_results = await loop.run_in_executor(
                None,
                run_two_stage_workflow,
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
            
            # Save updated job status
            save_job_to_file(job_id, jobs[job_id])
            
            # Send email notification
            user_email = jobs[job_id].get("user_email")
            user_name = jobs[job_id].get("user_name")
            filename = jobs[job_id].get("filename", "unknown")
            if user_email:
                try:
                    send_completion_email(job_id, user_email, user_name, filename)
                except Exception as email_error:
                    logger.warning(f"Failed to send completion email for job {job_id}: {email_error}")
            
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
            
            # Save updated job status
            save_job_to_file(job_id, jobs[job_id])
            
            # Send email notification
            user_email = jobs[job_id].get("user_email")
            user_name = jobs[job_id].get("user_name")
            filename = jobs[job_id].get("filename", "unknown")
            if user_email:
                try:
                    send_completion_email(job_id, user_email, user_name, filename)
                except Exception as email_error:
                    logger.warning(f"Failed to send completion email for job {job_id}: {email_error}")
            
            logger.info(f"Annotation completed (without UMAP) for job {job_id}")
        
    except Exception as e:
        logger.error(f"Annotation failed for job {job_id}: {str(e)}")
        jobs[job_id]["status"] = JobStatus.FAILED
        jobs[job_id]["end_time"] = datetime.now().isoformat()
        jobs[job_id]["error"] = str(e)
        
        # Save failed job status
        save_job_to_file(job_id, jobs[job_id])
    
    finally:
        # Clean up temporary directory
        if temp_dir and temp_dir.exists():
            try:
                shutil.rmtree(temp_dir)
                logger.info(f"Cleaned up temporary directory: {temp_dir}")
            except Exception as e:
                logger.warning(f"Failed to clean up temporary directory {temp_dir}: {e}")


@app.get("/api/job/{job_id}")
async def get_job_status(
    job_id: str,
    user: Dict[str, Any] = Depends(require_auth)
):
    """Get the status of a specific job (requires authentication)"""
    
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    # Verify user owns this job
    job_user_id = jobs[job_id].get("user_id")
    if job_user_id and job_user_id != user.get("google_sub"):
        raise HTTPException(status_code=403, detail="You don't have permission to access this job")
    
    return jobs[job_id]

@app.get("/api/jobs")
async def list_jobs(user: Dict[str, Any] = Depends(require_auth)):
    """List all jobs for the authenticated user (requires authentication)"""
    user_id = user.get("google_sub")
    user_jobs = [
        job for job in jobs.values()
        if not job.get("user_id") or job.get("user_id") == user_id
    ]
    return {"jobs": user_jobs}

@app.get("/api/download/{job_id}")
async def download_results(
    job_id: str,
    user: Dict[str, Any] = Depends(require_auth)
):
    """Download the results of a completed job (requires authentication)"""
    
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    # Verify user owns this job
    job_user_id = job.get("user_id")
    if job_user_id and job_user_id != user.get("google_sub"):
        raise HTTPException(status_code=403, detail="You don't have permission to access this job")
    
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
            plot_files = results.get("plot_files", [])
            preview_plot = next(
                (pf for pf in plot_files if pf.endswith("_umap_predictions_wolabel.png")),
                None
            )
            if preview_plot and Path(preview_plot).exists():
                response_data["umap_preview_url"] = f"/api/download/{job_id}/umap/preview"
                response_data["umap_preview_filename"] = Path(preview_plot).name
        
        return response_data
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error reading results: {str(e)}")

@app.get("/api/download/{job_id}/file")
async def download_file(
    job_id: str,
    user: Dict[str, Any] = Depends(require_auth)
):
    """Download the actual results file (requires authentication)"""
    
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    # Verify user owns this job
    job_user_id = job.get("user_id")
    if job_user_id and job_user_id != user.get("google_sub"):
        raise HTTPException(status_code=403, detail="You don't have permission to access this job")
    
    if job["status"] != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job is not completed yet")
    
    results = job["results"]
    predictions_file = Path(results["predictions_file"])
    
    if not predictions_file.exists():
        raise HTTPException(status_code=404, detail="Results file not found")
    
    return FileResponse(
        path=predictions_file,
        filename=predictions_file.name,
        media_type='text/csv'
    )

@app.get("/api/download/{job_id}/umap")
async def download_umap_results(
    job_id: str,
    user: Dict[str, Any] = Depends(require_auth)
):
    """Download UMAP results as a zip file (requires authentication)"""
    
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    # Verify user owns this job
    job_user_id = job.get("user_id")
    if job_user_id and job_user_id != user.get("google_sub"):
        raise HTTPException(status_code=403, detail="You don't have permission to access this job")
    
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
    
    return Response(
        content=zip_buffer.getvalue(),
        media_type="application/zip",
        headers={"Content-Disposition": f"attachment; filename=umap_results_{job_id}.zip"}
    )


@app.get("/api/download/{job_id}/umap/preview")
async def get_umap_preview(
    job_id: str,
    user: Dict[str, Any] = Depends(require_auth)
):
    """Return the UMAP predictions_wolabel plot for quick preview (requires authentication)"""

    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = jobs[job_id]
    
    # Verify user owns this job
    job_user_id = job.get("user_id")
    if job_user_id and job_user_id != user.get("google_sub"):
        raise HTTPException(status_code=403, detail="You don't have permission to access this job")
    
    if job["status"] != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job is not completed yet")

    results = job["results"]

    if "umap_results" not in results:
        raise HTTPException(status_code=404, detail="UMAP results not available")

    plot_files = results.get("plot_files", [])
    preview_plot = next(
        (pf for pf in plot_files if pf.endswith("_umap_predictions_wolabel.png")),
        None
    )

    if not preview_plot:
        raise HTTPException(status_code=404, detail="UMAP preview not available")

    preview_path = Path(preview_plot)
    if not preview_path.exists():
        raise HTTPException(status_code=404, detail="UMAP preview file missing")

    return FileResponse(
        path=preview_path,
        filename=preview_path.name,
        media_type="image/png"
    )


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
