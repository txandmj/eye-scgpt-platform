#!/usr/bin/env python3
"""
Startup script for the Eye-scGPT backend API.
This script sets up the environment and starts the FastAPI server.
"""

import os
import sys
import subprocess
from pathlib import Path

def setup_environment():
    """Setup the environment for running the backend"""
    
    # Create necessary directories
    directories = ['uploads', 'results', 'models']
    for directory in directories:
        Path(directory).mkdir(exist_ok=True)
    
    # Check if test model exists
    test_model_dir = Path('models/test_model')
    if not test_model_dir.exists():
        print("⚠️  Test model directory not found. Creating minimal structure...")
        # Run the setup script
        try:
            subprocess.run([sys.executable, 'setup_model.py'], check=True)
            print("✅ Test model structure created")
        except subprocess.CalledProcessError as e:
            print(f"❌ Failed to create test model structure: {e}")
            return False
    
    return True

def main():
    """Main startup function"""
    
    print("🚀 Starting Eye-scGPT Backend API...")
    print("=" * 50)
    
    # Setup environment
    if not setup_environment():
        print("❌ Failed to setup environment. Exiting.")
        sys.exit(1)
    
    # Check if required files exist
    required_files = ['main.py', 'step1_preprocess.py', 'step2_inference.py']
    for file in required_files:
        if not Path(file).exists():
            print(f"❌ Required file not found: {file}")
            sys.exit(1)
    
    print("✅ Environment setup complete")
    print("🌐 Starting FastAPI server on http://localhost:8000")
    print("📚 API documentation available at http://localhost:8000/docs")
    print("=" * 50)
    
    # Start the server
    try:
        subprocess.run([
            sys.executable, '-m', 'uvicorn', 
            'main:app', 
            '--host', '0.0.0.0', 
            '--port', '8000', 
            '--reload'
        ], check=True)
    except KeyboardInterrupt:
        print("\n👋 Server stopped by user")
    except subprocess.CalledProcessError as e:
        print(f"❌ Server failed to start: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
