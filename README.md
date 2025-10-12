# Eye-scGPT Annotation Platform

A web-based platform for automated cell type annotation of single-cell omics data using scGPT (Single-Cell GPT).

## 🚀 Features

- **Web-based Interface**: Easy-to-use React frontend for data upload and result download
- **Automated Processing**: Two-step pipeline for preprocessing and inference
- **RESTful API**: FastAPI backend with comprehensive endpoints
- **Docker Support**: Containerized deployment with Docker Compose
- **Job Tracking**: Real-time job status monitoring
- **File Management**: Secure file upload and result download

## 📁 Project Structure

```
eye-scgpt-platform/
├── backend/                 # FastAPI backend
│   ├── main.py             # Main API server
│   ├── step1_preprocess.py # Preprocessing script
│   ├── step2_inference.py  # Inference script
│   ├── utils/              # Utility modules
│   ├── models/             # Model storage
│   ├── uploads/            # Uploaded files
│   └── results/            # Processing results
├── frontend/               # React frontend
│   ├── src/
│   │   ├── components/     # React components
│   │   └── App.jsx         # Main app component
│   └── public/             # Static files
└── docker-compose.yml      # Docker configuration
```

## 🛠️ Setup Instructions

### Prerequisites

- Python 3.11+
- Node.js 16+
- Docker & Docker Compose (optional)

### Backend Setup

1. **Navigate to backend directory:**
   ```bash
   cd backend
   ```

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Setup test model (for development):**
   ```bash
   python setup_model.py
   ```

4. **Start the backend server:**
   ```bash
   python start.py
   ```
   
   Or manually:
   ```bash
   uvicorn main:app --host 0.0.0.0 --port 8000 --reload
   ```

The API will be available at:
- **API**: http://localhost:8000
- **Documentation**: http://localhost:8000/docs
- **Health Check**: http://localhost:8000/health

### Frontend Setup

1. **Navigate to frontend directory:**
   ```bash
   cd frontend
   ```

2. **Install dependencies:**
   ```bash
   npm install
   ```

3. **Start the development server:**
   ```bash
   npm start
   ```

The frontend will be available at: http://localhost:3000

### Docker Setup (Recommended)

1. **Build and start all services:**
   ```bash
   docker-compose up --build
   ```

2. **Access the application:**
   - **Frontend**: http://localhost:80
   - **Backend API**: http://localhost:8000

## 📊 Usage

### 1. Upload Data

1. Navigate to the **Upload** tab
2. Select a `.h5ad` file (AnnData format)
3. Click **Upload & Annotate**
4. Wait for processing to complete

### 2. Download Results

1. Navigate to the **Download** tab
2. Enter your Job ID (provided after upload)
3. Click **Fetch Results**
4. Download the predictions CSV file

### 3. API Usage

#### Upload a file:
```bash
curl -X POST "http://localhost:8000/api/upload" \
  -H "Content-Type: multipart/form-data" \
  -F "file=@your_data.h5ad"
```

#### Start annotation:
```bash
curl -X POST "http://localhost:8000/api/annotate/{job_id}"
```

#### Check job status:
```bash
curl "http://localhost:8000/api/job/{job_id}"
```

#### Download results:
```bash
curl "http://localhost:8000/api/download/{job_id}/file" -o predictions.csv
```

## 🔧 Configuration

### Backend Configuration

Key parameters in `step1_preprocess.py` and `step2_inference.py`:

- `load_model`: Path to the trained scGPT model
- `batch_size`: Processing batch size
- `n_hvg`: Number of highly variable genes
- `max_seq_len`: Maximum sequence length

### Model Requirements

The `models/test_model/` directory should contain:

- `best_model.pt`: Trained model weights
- `args.json`: Model configuration
- `vocab.json`: Gene vocabulary
- `dev_train_args.yml`: Training arguments
- `id2type.json`: Cell type mappings

## 🐛 Troubleshooting

### Common Issues

1. **Model not found error:**
   - Ensure the `test_model` directory exists in `backend/models/`
   - Run `python setup_model.py` to create a minimal structure

2. **Memory issues:**
   - Reduce `batch_size` in the processing scripts
   - Increase Docker memory limits in `docker-compose.yml`

3. **File upload errors:**
   - Check file format (only `.h5ad` files supported)
   - Ensure file size is within limits (5GB default)

4. **Processing timeouts:**
   - Large datasets may take 10+ minutes
   - Check job status via API or frontend

### Logs

- **Backend logs**: Check terminal output or Docker logs
- **Processing logs**: Available in `results/{job_id}/step2/run.log`

## 🔒 Security Notes

- The current setup is for development/demo purposes
- For production deployment:
  - Implement proper authentication
  - Use HTTPS
  - Set up proper file storage (S3, etc.)
  - Implement rate limiting
  - Add input validation

## 📝 API Documentation

Full API documentation is available at:
- **Swagger UI**: http://localhost:8000/docs
- **ReDoc**: http://localhost:8000/redoc

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🙏 Acknowledgments

- [scGPT](https://github.com/bowang-lab/scGPT) - The core model for single-cell annotation
- [Scanpy](https://scanpy.readthedocs.io/) - Single-cell analysis library
- [FastAPI](https://fastapi.tiangolo.com/) - Modern Python web framework
- [React](https://reactjs.org/) - Frontend framework

## 📞 Support

For issues and questions:
- Create an issue in the GitHub repository
- Check the troubleshooting section above
- Review the API documentation

---

**Note**: This platform is designed for research purposes. Ensure you have proper permissions and ethical approvals for any data processing.
