#!/usr/bin/env python3
"""
Test script for the Eye-scGPT API endpoints.
This script tests the basic functionality of the API.
"""

import requests
import time
import os
from pathlib import Path

API_BASE_URL = "http://localhost:8000"

def test_health_check():
    """Test the health check endpoint"""
    print("🔍 Testing health check...")
    try:
        response = requests.get(f"{API_BASE_URL}/health")
        if response.status_code == 200:
            print("✅ Health check passed")
            return True
        else:
            print(f"❌ Health check failed: {response.status_code}")
            return False
    except Exception as e:
        print(f"❌ Health check error: {e}")
        return False

def test_upload():
    """Test file upload endpoint"""
    print("🔍 Testing file upload...")
    
    # Check if we have a test file
    test_file = Path("outdir_step1/input_data.h5ad")
    if not test_file.exists():
        print("⚠️  No test file found. Creating a dummy test...")
        # Create a minimal test file
        test_file.parent.mkdir(exist_ok=True)
        # This would normally be a real h5ad file
        print("📝 Note: You need to provide a real .h5ad file for testing")
        return None
    
    try:
        with open(test_file, "rb") as f:
            files = {"file": ("test_data.h5ad", f, "application/octet-stream")}
            response = requests.post(f"{API_BASE_URL}/api/upload", files=files)
        
        if response.status_code == 200:
            data = response.json()
            print(f"✅ Upload successful: Job ID {data['job_id']}")
            return data['job_id']
        else:
            print(f"❌ Upload failed: {response.status_code} - {response.text}")
            return None
    except Exception as e:
        print(f"❌ Upload error: {e}")
        return None

def test_annotation(job_id):
    """Test annotation endpoint"""
    print(f"🔍 Testing annotation for job {job_id}...")
    
    try:
        response = requests.post(f"{API_BASE_URL}/api/annotate/{job_id}")
        if response.status_code == 200:
            print("✅ Annotation started successfully")
            return True
        else:
            print(f"❌ Annotation failed: {response.status_code} - {response.text}")
            return False
    except Exception as e:
        print(f"❌ Annotation error: {e}")
        return False

def test_job_status(job_id):
    """Test job status endpoint"""
    print(f"🔍 Testing job status for {job_id}...")
    
    try:
        response = requests.get(f"{API_BASE_URL}/api/job/{job_id}")
        if response.status_code == 200:
            data = response.json()
            print(f"✅ Job status: {data['status']}")
            return data
        else:
            print(f"❌ Job status failed: {response.status_code} - {response.text}")
            return None
    except Exception as e:
        print(f"❌ Job status error: {e}")
        return None

def test_download(job_id):
    """Test download endpoint"""
    print(f"🔍 Testing download for job {job_id}...")
    
    try:
        response = requests.get(f"{API_BASE_URL}/api/download/{job_id}")
        if response.status_code == 200:
            print("✅ Download endpoint accessible")
            return True
        else:
            print(f"❌ Download failed: {response.status_code} - {response.text}")
            return False
    except Exception as e:
        print(f"❌ Download error: {e}")
        return False

def main():
    """Main test function"""
    print("🧪 Eye-scGPT API Test Suite")
    print("=" * 50)
    
    # Test 1: Health check
    if not test_health_check():
        print("❌ API server is not running. Please start the backend first.")
        print("   Run: python start.py")
        return
    
    # Test 2: Upload
    job_id = test_upload()
    if not job_id:
        print("⚠️  Skipping remaining tests due to upload failure")
        return
    
    # Test 3: Job status (before annotation)
    test_job_status(job_id)
    
    # Test 4: Start annotation
    if test_annotation(job_id):
        # Test 5: Monitor job status
        print("⏳ Monitoring job status...")
        for i in range(10):  # Check for up to 10 iterations
            time.sleep(2)
            job_data = test_job_status(job_id)
            if job_data and job_data['status'] in ['completed', 'failed']:
                break
        
        # Test 6: Download (if completed)
        if job_data and job_data['status'] == 'completed':
            test_download(job_id)
        else:
            print("⚠️  Job not completed yet. You can check status later.")
    
    print("=" * 50)
    print("🎉 Test suite completed!")

if __name__ == "__main__":
    main()
