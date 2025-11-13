#!/usr/bin/env python3
"""
Test script for the complete workflow: upload -> annotation -> UMAP processing -> download
"""

import requests
import time
import sys
from pathlib import Path

API_BASE_URL = "http://localhost:8000"

def test_health_check():
    """Test if the API is running"""
    print("🔍 Testing API health...")
    try:
        response = requests.get(f"{API_BASE_URL}/health")
        if response.status_code == 200:
            print("✅ API is healthy")
            return True
        else:
            print(f"❌ API health check failed: {response.status_code}")
            return False
    except requests.exceptions.RequestException as e:
        print(f"❌ API not reachable: {e}")
        print("💡 Make sure the backend is running: python start.py")
        return False

def test_upload_and_workflow():
    """Test the complete workflow"""
    print("\n🧪 Testing Complete Workflow")
    print("=" * 50)
    
    # Check if we have a test file
    test_file = Path("../EVAL_snRNA_no_enriched.h5ad")
    if not test_file.exists():
        print("❌ Test file not found. Please ensure EVAL_snRNA_no_enriched.h5ad exists in the project root.")
        return False
    
    try:
        # Step 1: Upload file
        print("📤 Step 1: Uploading .h5ad file...")
        with open(test_file, "rb") as f:
            files = {"file": (test_file.name, f, "application/octet-stream")}
            response = requests.post(f"{API_BASE_URL}/api/upload", files=files)
        
        if response.status_code != 200:
            print(f"❌ Upload failed: {response.status_code} - {response.text}")
            return False
        
        upload_data = response.json()
        job_id = upload_data['job_id']
        print(f"✅ Upload successful: Job ID {job_id}")
        
        # Step 2: Start annotation
        print("\n🔬 Step 2: Starting eye-scgpt annotation...")
        response = requests.post(f"{API_BASE_URL}/api/annotate/{job_id}")
        
        if response.status_code != 200:
            print(f"❌ Annotation start failed: {response.status_code} - {response.text}")
            return False
        
        annotation_data = response.json()
        print(f"✅ Annotation started: {annotation_data['message']}")
        
        # Step 3: Monitor progress
        print("\n⏳ Step 3: Monitoring progress...")
        print("This may take several minutes for annotation + UMAP processing...")
        
        while True:
            response = requests.get(f"{API_BASE_URL}/api/job/{job_id}")
            if response.status_code != 200:
                print(f"❌ Status check failed: {response.status_code}")
                return False
            
            job_data = response.json()
            status = job_data['status']
            
            print(f"📊 Status: {status}")
            
            if status == "completed":
                print("✅ Processing completed successfully!")
                
                # Check results
                results = job_data.get('results', {})
                print(f"📄 Predictions file: {results.get('predictions_file', 'N/A')}")
                
                if 'umap_results' in results:
                    print(f"🖼️  UMAP plots generated: {results.get('num_plots', 0)}")
                    print(f"📁 Enhanced h5ad file: {results.get('enriched_h5ad_file', 'N/A')}")
                else:
                    print("⚠️  UMAP processing failed or not available")
                
                return job_id
                
            elif status == "failed":
                error = job_data.get('error', 'Unknown error')
                print(f"❌ Processing failed: {error}")
                return False
            
            elif status == "processing":
                print("⏳ Still processing... (waiting 30 seconds)")
                time.sleep(30)
            
            else:
                print(f"⚠️  Unknown status: {status}")
                time.sleep(10)
    
    except Exception as e:
        print(f"❌ Workflow test failed: {e}")
        return False

def test_download_results(job_id):
    """Test downloading results"""
    print(f"\n📥 Step 4: Testing download for job {job_id}")
    
    try:
        # Test CSV download
        print("📄 Testing CSV download...")
        response = requests.get(f"{API_BASE_URL}/api/download/{job_id}")
        
        if response.status_code == 200:
            data = response.json()
            print(f"✅ CSV download info retrieved")
            print(f"   Filename: {data['filename']}")
            print(f"   UMAP available: {data.get('umap_available', False)}")
            if data.get('umap_available'):
                print(f"   Number of plots: {data.get('num_plots', 0)}")
        else:
            print(f"❌ CSV download failed: {response.status_code}")
        
        # Test UMAP download if available
        if data.get('umap_available'):
            print("\n🖼️  Testing UMAP download...")
            response = requests.get(f"{API_BASE_URL}/api/download/{job_id}/umap")
            
            if response.status_code == 200:
                # Save the zip file
                zip_file = f"umap_results_{job_id}.zip"
                with open(zip_file, 'wb') as f:
                    f.write(response.content)
                print(f"✅ UMAP results downloaded: {zip_file}")
            else:
                print(f"❌ UMAP download failed: {response.status_code}")
        
        return True
        
    except Exception as e:
        print(f"❌ Download test failed: {e}")
        return False

def main():
    """Main test function"""
    print("🔬 Complete Workflow Test")
    print("=" * 60)
    print("This test will:")
    print("1. Upload a .h5ad file")
    print("2. Run eye-scgpt annotation")
    print("3. Run UMAP processing")
    print("4. Download results (CSV + PNG plots)")
    print("=" * 60)
    
    # Check API health
    if not test_health_check():
        return
    
    # Run complete workflow
    job_id = test_upload_and_workflow()
    if not job_id:
        return
    
    # Test downloads
    test_download_results(job_id)
    
    print("\n🎉 Complete workflow test finished!")
    print(f"📊 Job ID: {job_id}")
    print("💡 Check the results directory for generated files")

if __name__ == "__main__":
    main()
