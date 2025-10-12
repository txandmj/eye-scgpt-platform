#!/usr/bin/env python3
"""
Startup script for the entire Eye-scGPT platform.
This script helps you start both frontend and backend services.
"""

import subprocess
import sys
import time
import os
from pathlib import Path

def check_requirements():
    """Check if required tools are available"""
    print("🔍 Checking requirements...")
    
    # Check Python
    try:
        result = subprocess.run([sys.executable, "--version"], capture_output=True, text=True)
        print(f"✅ Python: {result.stdout.strip()}")
    except Exception as e:
        print(f"❌ Python not found: {e}")
        return False
    
    # Check Node.js
    try:
        result = subprocess.run(["node", "--version"], capture_output=True, text=True)
        print(f"✅ Node.js: {result.stdout.strip()}")
    except Exception as e:
        print(f"❌ Node.js not found: {e}")
        print("   Please install Node.js from https://nodejs.org/")
        return False
    
    # Check npm
    try:
        result = subprocess.run(["npm", "--version"], capture_output=True, text=True)
        print(f"✅ npm: {result.stdout.strip()}")
    except Exception as e:
        print(f"❌ npm not found: {e}")
        return False
    
    return True

def setup_backend():
    """Setup and start the backend"""
    print("\n🔧 Setting up backend...")
    
    backend_dir = Path("backend")
    if not backend_dir.exists():
        print("❌ Backend directory not found")
        return False
    
    # Install Python dependencies
    print("📦 Installing Python dependencies...")
    try:
        subprocess.run([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"], 
                      cwd=backend_dir, check=True)
        print("✅ Python dependencies installed")
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to install Python dependencies: {e}")
        return False
    
    # Setup test model
    print("🤖 Setting up test model...")
    try:
        subprocess.run([sys.executable, "setup_model.py"], cwd=backend_dir, check=True)
        print("✅ Test model setup complete")
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to setup test model: {e}")
        return False
    
    return True

def setup_frontend():
    """Setup the frontend"""
    print("\n🔧 Setting up frontend...")
    
    frontend_dir = Path("frontend")
    if not frontend_dir.exists():
        print("❌ Frontend directory not found")
        return False
    
    # Install Node.js dependencies
    print("📦 Installing Node.js dependencies...")
    try:
        subprocess.run(["npm", "install"], cwd=frontend_dir, check=True)
        print("✅ Node.js dependencies installed")
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to install Node.js dependencies: {e}")
        return False
    
    return True

def start_services():
    """Start both backend and frontend services"""
    print("\n🚀 Starting services...")
    
    backend_dir = Path("backend")
    frontend_dir = Path("frontend")
    
    try:
        # Start backend
        print("🔧 Starting backend server...")
        backend_process = subprocess.Popen(
            [sys.executable, "start.py"],
            cwd=backend_dir
        )
        
        # Wait a bit for backend to start
        time.sleep(3)
        
        # Start frontend
        print("🎨 Starting frontend server...")
        frontend_process = subprocess.Popen(
            ["npm", "start"],
            cwd=frontend_dir
        )
        
        print("\n" + "="*60)
        print("🎉 Eye-scGPT Platform is starting!")
        print("="*60)
        print("📊 Frontend: http://localhost:3000")
        print("🔧 Backend API: http://localhost:8000")
        print("📚 API Docs: http://localhost:8000/docs")
        print("="*60)
        print("\n💡 Tips:")
        print("   - Wait for both services to fully start")
        print("   - Check the terminal output for any errors")
        print("   - Use Ctrl+C to stop both services")
        print("\n⏳ Starting services... Please wait...")
        
        # Wait for processes
        try:
            backend_process.wait()
        except KeyboardInterrupt:
            print("\n🛑 Stopping services...")
            backend_process.terminate()
            frontend_process.terminate()
            print("👋 Services stopped")
            
    except Exception as e:
        print(f"❌ Failed to start services: {e}")
        return False
    
    return True

def main():
    """Main startup function"""
    print("🔬 Eye-scGPT Platform Startup")
    print("=" * 50)
    
    # Check requirements
    if not check_requirements():
        print("\n❌ Requirements check failed. Please install missing tools.")
        return
    
    # Setup backend
    if not setup_backend():
        print("\n❌ Backend setup failed.")
        return
    
    # Setup frontend
    if not setup_frontend():
        print("\n❌ Frontend setup failed.")
        return
    
    # Start services
    if not start_services():
        print("\n❌ Failed to start services.")
        return

if __name__ == "__main__":
    main()
