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


## 🛠️ Quick Start

### Option 1: Docker (Recommended)

```bash
# Clone the repository
git clone <repository-url>
cd eye-scgpt-platform

# Start the entire platform
docker-compose up --build

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