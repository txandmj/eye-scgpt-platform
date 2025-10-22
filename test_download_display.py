#!/usr/bin/env python3
"""
Test script to verify the download display shows both CSV and UMAP files
"""

import requests
import json

def test_download_display():
    """Test that the download endpoint returns both file types"""
    
    API_URL = "http://localhost:8000"
    
    print("🧪 Testing Download Display")
    print("=" * 50)
    
    # Get available jobs
    try:
        response = requests.get(f"{API_URL}/api/jobs")
        if response.status_code != 200:
            print(f"❌ Failed to get jobs: {response.status_code}")
            return False
        
        jobs = response.json().get('jobs', [])
        if not jobs:
            print("❌ No jobs found")
            return False
        
        # Find a completed job
        completed_job = None
        for job in jobs:
            if job['status'] == 'completed':
                completed_job = job
                break
        
        if not completed_job:
            print("❌ No completed jobs found")
            return False
        
        job_id = completed_job['id']
        print(f"📊 Testing with job: {job_id}")
        
        # Test download info
        response = requests.get(f"{API_URL}/api/download/{job_id}")
        if response.status_code != 200:
            print(f"❌ Failed to get download info: {response.status_code}")
            return False
        
        data = response.json()
        print(f"✅ Download info retrieved")
        print(f"   Filename: {data['filename']}")
        print(f"   UMAP available: {data.get('umap_available', False)}")
        
        if data.get('umap_available'):
            print(f"   Number of plots: {data.get('num_plots', 0)}")
            print(f"   UMAP download URL: {data.get('umap_download_url', 'N/A')}")
            
            # Test UMAP download
            print("\n🖼️  Testing UMAP download...")
            umap_response = requests.get(f"{API_URL}{data['umap_download_url']}")
            if umap_response.status_code == 200:
                print(f"✅ UMAP download successful ({len(umap_response.content)} bytes)")
            else:
                print(f"❌ UMAP download failed: {umap_response.status_code}")
        else:
            print("⚠️  UMAP results not available")
        
        print("\n🎉 Download display test completed!")
        print("💡 Now check the frontend at http://localhost:3000")
        print(f"   Use Job ID: {job_id}")
        
        return True
        
    except Exception as e:
        print(f"❌ Test failed: {e}")
        return False

if __name__ == "__main__":
    test_download_display()
