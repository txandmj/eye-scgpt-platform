# 🚀 Quick Start Guide

Get your Eye-scGPT platform up and running in minutes!

## Option 1: Automated Setup (Recommended)

Run the automated setup script:

```bash
python start_platform.py
```

This will:
- ✅ Check all requirements
- 📦 Install dependencies
- 🤖 Setup test model
- 🚀 Start both frontend and backend

## Option 2: Manual Setup

### 1. Backend Setup
```bash
cd backend
pip install -r requirements.txt
python setup_model.py
python start.py
```

### 2. Frontend Setup (New Terminal)
```bash
cd frontend
npm install
npm start
```

## Option 3: Docker (Easiest)

```bash
docker-compose up --build
```

## 🌐 Access Your Platform

Once started, access:
- **Frontend**: http://localhost:3000 (or http://localhost:80 with Docker)
- **Backend API**: http://localhost:8000
- **API Docs**: http://localhost:8000/docs

## 📊 Test Your Setup

1. **Upload Tab**: Upload a `.h5ad` file
2. **Download Tab**: Enter Job ID to get results
3. **API Test**: Run `python backend/test_api.py`

## 🔧 Troubleshooting

### Backend won't start?
- Check if port 8000 is free
- Ensure Python 3.11+ is installed
- Run `pip install -r backend/requirements.txt`

### Frontend won't start?
- Check if port 3000 is free
- Ensure Node.js 16+ is installed
- Run `npm install` in frontend directory

### Model errors?
- Run `python backend/setup_model.py`
- Ensure `backend/models/test_model/` exists

## 📝 Next Steps

1. **Add your real model**: Replace files in `backend/models/test_model/`
2. **Configure parameters**: Edit `step1_preprocess.py` and `step2_inference.py`
3. **Test with real data**: Upload your `.h5ad` files
4. **Deploy**: Use Docker or cloud services

## 🆘 Need Help?

- Check the full [README.md](README.md)
- Review API docs at http://localhost:8000/docs
- Run the test script: `python backend/test_api.py`

---

**Ready to annotate some cells? 🧬**
