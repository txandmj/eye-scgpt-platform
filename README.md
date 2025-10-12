# Eye-scGPT Annotation Platform

A streamlined web-based platform for automated cell type annotation of single-cell omics data using scGPT (Single-Cell GPT). The platform provides a direct workflow: upload data → process → get results.

## 🚀 Features

- **Simple Web Interface**: Clean React frontend for easy data upload and result download
- **Direct Processing Pipeline**: Streamlined two-step workflow (preprocessing → inference)
- **RESTful API**: FastAPI backend with comprehensive endpoints
- **Docker Support**: One-command deployment with Docker Compose
- **Real-time Job Tracking**: Monitor processing status and progress
- **Clean Results Storage**: Only final predictions stored in results folder
- **Automatic Cleanup**: Temporary files automatically removed after processing

## 📁 Project Structure

```
eye-scgpt-platform/
├── backend/                    # FastAPI backend
│   ├── main.py                # Main API server with simplified workflow
│   ├── step1_preprocess.py    # Preprocessing script (uses env variables)
│   ├── step2_inference.py     # Inference script (uses env variables)
│   ├── setup_model.py         # Model setup utility
│   ├── utils/                 # Utility modules and helpers
│   ├── models/                # Model storage
│   │   └── test_model/        # Pre-trained scGPT model
│   ├── uploads/               # Temporary upload storage
│   └── results/               # Final results only (predictions.csv)
├── frontend/                  # React frontend
│   ├── src/
│   │   ├── components/        # React components
│   │   └── App.jsx           # Main application
│   └── public/               # Static files
├── docker-compose.yml         # Docker configuration
└── start_platform.py         # Platform startup script
```

## 🛠️ Quick Start

### Option 1: Docker (Recommended)

```bash
# Clone the repository
git clone <repository-url>
cd eye-scgpt-platform

# Start the entire platform
docker-compose up --build

# Access the application
# Frontend: http://localhost:80
# Backend API: http://localhost:8000
```

### Option 2: Manual Setup

#### Backend Setup

```bash
cd backend

# Install dependencies
pip install -r requirements.txt

# Setup test model
python setup_model.py

# Start backend server
python start.py
# OR: uvicorn main:app --host 0.0.0.0 --port 8000 --reload
```

#### Frontend Setup

```bash
cd frontend

# Install dependencies
npm install

# Start development server
npm start
# Frontend: http://localhost:3000
```

## 📊 How It Works

### Simplified Workflow

1. **Upload**: User uploads `.h5ad` file via web interface
2. **Process**: Backend runs original `step1_preprocess.py` and `step2_inference.py` directly
3. **Store**: Only final `predictions.csv` is saved to `results/{job_id}/`
4. **Clean**: Temporary processing files are automatically removed
5. **Download**: User downloads final results


## 🎯 Usage

### Web Interface

1. **Upload Data**:
   - Go to Upload tab
   - Select `.h5ad` file
   - Click "Upload & Annotate"
   - Note your Job ID

2. **Monitor Progress**:
   - Check job status in real-time
   - Processing typically takes 5-15 minutes

3. **Download Results**:
   - Go to Download tab
   - Enter Job ID
   - Download `predictions.csv`

### API Usage

```bash
# Upload file
curl -X POST "http://localhost:8000/api/upload" \
  -H "Content-Type: multipart/form-data" \
  -F "file=@your_data.h5ad"

# Start annotation
curl -X POST "http://localhost:8000/api/annotate/{job_id}"

# Check status
curl "http://localhost:8000/api/job/{job_id}"

# Download results
curl "http://localhost:8000/api/download/{job_id}/file" -o predictions.csv
```

## ⚙️ Configuration

### Environment Variables

The processing scripts now use environment variables for flexibility:

**Step 1 (Preprocessing)**:
- `DATASET_DIRECTORY`: Path to input `.h5ad` file
- `SAVE_DIR`: Output directory for preprocessing
- `LOAD_MODEL`: Path to trained model

**Step 2 (Inference)**:
- `LOAD_MODEL`: Path to trained model
- `TRAIN_ARGS`: Path to training arguments from step 1
- `SAVE_DIR`: Output directory for final results

### Model Requirements

Ensure `backend/models/test_model/` contains:
- `best_model.pt` - Trained model weights
- `vocab.json` - Gene vocabulary
- `id2type.json` - Cell type mappings
- `dev_train_args.yml` - Training configuration

## 🔧 API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/` | GET | Platform information |
| `/health` | GET | Health check |
| `/api/upload` | POST | Upload `.h5ad` file |
| `/api/annotate/{job_id}` | POST | Start annotation process |
| `/api/job/{job_id}` | GET | Get job status |
| `/api/jobs` | GET | List all jobs |
| `/api/download/{job_id}` | GET | Download results (JSON) |
| `/api/download/{job_id}/file` | GET | Download results (CSV file) |

## 🐛 Troubleshooting

### Common Issues

**Model Not Found**:
```bash
# Ensure model directory exists
ls backend/models/test_model/
# Run setup if missing
python backend/setup_model.py
```

**Memory Issues**:
- Reduce `batch_size` in scripts
- Increase Docker memory limits
- Process smaller datasets

**File Upload Errors**:
- Only `.h5ad` files supported
- Check file size (default 5GB limit)
- Verify file format

**Processing Timeouts**:
- Large datasets: 10+ minutes
- Check job status via API
- Monitor system resources

### Logs and Debugging

- **Backend logs**: Terminal output or `docker logs`
- **API docs**: http://localhost:8000/docs
- **Job status**: Use `/api/job/{job_id}` endpoint


## 📈 Performance

- **Typical processing time**: 20-45 minutes(CPU) or 2 minutes(GPU) per dataset
- **Supported file sizes**: Up to 5GB
- **Concurrent jobs**: Limited by system resources
- **Memory usage**: ~2-8GB depending on dataset size

## 🤝 Contributing

1. Fork the repository
2. Create feature branch: `git checkout -b feature/amazing-feature`
3. Make changes and test
4. Commit: `git commit -m 'Add amazing feature'`
5. Push: `git push origin feature/amazing-feature`
6. Open Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- [scGPT](https://github.com/bowang-lab/scGPT) - Core single-cell annotation model
- [Scanpy](https://scanpy.readthedocs.io/) - Single-cell analysis library
- [FastAPI](https://fastapi.tiangolo.com/) - Modern Python web framework
- [React](https://reactjs.org/) - Frontend framework

## 📞 Support

- **Issues**: Create GitHub issue
- **Documentation**: Check API docs at `/docs`
- **Troubleshooting**: See troubleshooting section above

---

**Note**: This platform is designed for research purposes. Ensure proper permissions and ethical approvals for data processing.