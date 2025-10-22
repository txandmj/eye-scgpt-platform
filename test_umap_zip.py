#!/usr/bin/env python3
"""
Test script to verify UMAP ZIP functionality
"""

import zipfile
import io
from pathlib import Path

def test_umap_zip_creation():
    """Test creating a UMAP ZIP file with the existing results"""
    
    # Use the existing UMAP results
    umap_dir = Path("backend/results/ace89739-332e-44a4-a128-ed804b9d94e8/umap")
    
    if not umap_dir.exists():
        print(f"❌ UMAP directory not found: {umap_dir}")
        return False
    
    print(f"📁 Testing UMAP ZIP creation from: {umap_dir}")
    
    # List files in the directory
    files = list(umap_dir.iterdir())
    png_files = [f for f in files if f.suffix == '.png']
    h5ad_files = [f for f in files if f.suffix == '.h5ad']
    
    print(f"📊 Found {len(files)} total files")
    print(f"🖼️  PNG files: {len(png_files)}")
    print(f"🔬 H5AD files: {len(h5ad_files)}")
    
    if len(png_files) == 0:
        print("❌ No PNG files found!")
        return False
    
    # Create ZIP file
    zip_buffer = io.BytesIO()
    
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        # Add h5ad file
        for h5ad_file in h5ad_files:
            zip_file.write(h5ad_file, h5ad_file.name)
            print(f"✅ Added {h5ad_file.name}")
        
        # Add PNG files
        for png_file in png_files:
            zip_file.write(png_file, png_file.name)
            print(f"✅ Added {png_file.name}")
    
    zip_buffer.seek(0)
    zip_content = zip_buffer.getvalue()
    
    print(f"\n📦 ZIP file created: {len(zip_content)} bytes")
    
    # Save the ZIP file
    output_zip = "test_umap_results.zip"
    with open(output_zip, 'wb') as f:
        f.write(zip_content)
    
    print(f"💾 Saved to: {output_zip}")
    
    # Verify the ZIP file
    with zipfile.ZipFile(output_zip, 'r') as zip_file:
        file_list = zip_file.namelist()
        print(f"\n📋 ZIP contents ({len(file_list)} files):")
        for i, filename in enumerate(file_list[:10]):  # Show first 10 files
            print(f"  {i+1}. {filename}")
        if len(file_list) > 10:
            print(f"  ... and {len(file_list) - 10} more files")
    
    return True

if __name__ == "__main__":
    success = test_umap_zip_creation()
    if success:
        print("\n🎉 UMAP ZIP test completed successfully!")
    else:
        print("\n❌ UMAP ZIP test failed!")
