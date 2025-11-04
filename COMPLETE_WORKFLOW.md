# Complete Eye-scGPT + UMAP Workflow

This document describes the complete workflow that processes .h5ad files through eye-scgpt annotation and UMAP visualization.

## 🎯 Overview

The platform now provides a **complete end-to-end workflow**:

1. **Upload** .h5ad file
2. **Eye-scGPT Annotation** → generates predictions.csv
3. **UMAP Processing** → generates PNG plots
4. **Download Results** → CSV + PNG plots

## 🚀 Quick Start

### 1. Start the Platform

```bash
# Start both backend and frontend
python start_platform.py
```

Or start individually:
```bash
# Backend only
cd backend && python start.py

# Frontend only (in another terminal)
cd frontend && npm start
```

### 2. Access the Platform

- **Frontend**: http://localhost:3000
- **Backend API**: http://localhost:8000
- **API Documentation**: http://localhost:8000/docs

## 📋 Complete Workflow

### Step 1: Upload .h5ad File
```bash
curl -X POST "http://localhost:8000/api/upload" \
     -F "file=@your_data.h5ad"
```

**Response:**
```json
{
  "job_id": "uuid-here",
  "status": "pending",
  "message": "File uploaded successfully"
}
```

### Step 2: Start Annotation
```bash
curl -X POST "http://localhost:8000/api/annotate/{job_id}"
```

**Response:**
```json
{
  "job_id": "uuid-here",
  "status": "processing",
  "message": "Annotation process started"
}
```

### Step 3: Monitor Progress
```bash
curl "http://localhost:8000/api/job/{job_id}"
```

**Response (Processing):**
```json
{
  "id": "uuid-here",
  "status": "processing",
  "filename": "your_data.h5ad",
  "start_time": "2024-01-01T12:00:00"
}
```

**Response (Completed):**
```json
{
  "id": "uuid-here",
  "status": "completed",
  "results": {
    "predictions_file": "results/uuid/predictions.csv",
    "umap_results": {
      "enriched_h5ad_file": "results/uuid/umap/data_enriched.h5ad",
      "plot_files": ["plot1.png", "plot2.png", ...],
      "num_plots": 39
    }
  }
}
```

### Step 4: Download Results

#### Download CSV (Predictions)
```bash
curl "http://localhost:8000/api/download/{job_id}/file" -o predictions.csv
```

#### Download UMAP Results (ZIP with all plots)
```bash
curl "http://localhost:8000/api/download/{job_id}/umap" -o umap_results.zip
```

## 📊 Generated Files

### After Complete Processing:

```
results/{job_id}/
├── predictions.csv                    # Eye-scGPT predictions
└── umap/
    ├── data_enriched.h5ad            # Enhanced data with metadata
    ├── data_umap_celltype_wolabel.png
    ├── data_umap_celltype_ondata.png
    ├── data_umap_celltype_fontline.png
    ├── data_umap_predictions_wolabel.png
    ├── data_umap_predictions_ondata.png
    ├── data_umap_predictions_fontline.png
    └── ... (more plots for each metadata column)
```

### File Types:

1. **predictions.csv** - Eye-scGPT annotation results
2. **data_enriched.h5ad** - Original data + predictions metadata
3. **UMAP plots** - PNG files with different visualization styles:
   - `_wolabel.png` - Clean UMAP without labels
   - `_ondata.png` - UMAP with labels on data points
   - `_fontline.png` - UMAP with labels and text outline

## 🧪 Testing

### Test Complete Workflow
```bash
cd backend
python test_complete_workflow.py
```

This will:
1. Upload a test file
2. Run annotation + UMAP processing
3. Download and verify results

### Test Individual Components
```bash
# Test simple UMAP processor
python test_simple_umap.py

# Test API endpoints
python test_api.py
```

## 🔧 API Endpoints

### Core Endpoints
- `POST /api/upload` - Upload .h5ad file
- `POST /api/annotate/{job_id}` - Start annotation process
- `GET /api/job/{job_id}` - Get job status
- `GET /api/jobs` - List all jobs

### Download Endpoints
- `GET /api/download/{job_id}` - Get download info
- `GET /api/download/{job_id}/file` - Download CSV file
- `GET /api/download/{job_id}/umap` - Download UMAP results (ZIP)

### Health Check
- `GET /health` - API health status
- `GET /` - API info

## 🎨 Frontend Integration

The frontend at http://localhost:3000 provides a web interface for:

- File upload
- Progress monitoring
- Result visualization
- Download management

## 📝 Processing Details

### Eye-scGPT Annotation
- Runs step1_preprocess.py and step2_inference.py
- Generates predictions.csv with cell type predictions
- Takes 5-30 minutes depending on data size

### UMAP Processing
- Automatically runs after annotation completes
- Uses simple two-stage processor:
  - Stage 1: Merge CSV metadata into h5ad
  - Stage 2: Generate UMAP plots using existing coordinates
- Generates 2-3 plots per metadata column
- Takes 1-5 minutes depending on data size

### Error Handling
- If UMAP processing fails, annotation results are still available
- Detailed error messages in job status
- Graceful degradation (annotation without UMAP plots)

## 🚨 Requirements

### Backend Requirements
- Python 3.8+
- All dependencies in requirements.txt
- Test model in models/test_model/

### Frontend Requirements
- Node.js 16+
- All dependencies in package.json

### Data Requirements
- Input .h5ad file must contain UMAP coordinates (`X_umap` in .obsm)
- File should be properly formatted for scanpy

## 🎉 Success!

The platform now provides a **complete end-to-end solution**:

✅ Upload .h5ad files  
✅ Run eye-scgpt annotation  
✅ Generate UMAP visualizations  
✅ Download CSV + PNG results  
✅ Web interface for easy use  

Everything is integrated and ready to use! 🚀
