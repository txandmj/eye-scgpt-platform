import React, { useState } from 'react';
import axios from 'axios';
import { getAuthHeaders } from '../auth';

const rawBase = process.env.REACT_APP_API_BASE || process.env.REACT_APP_API_URL || 'http://localhost:8000';
const normalizeBase = (url) => {
  if (!url) return 'http://localhost:8000';
  const trimmed = url.trim();
  // Handle values like ":8000" → protocol + hostname + port
  if (trimmed.startsWith(':')) {
    return `${window.location.protocol}//${window.location.hostname}${trimmed}`;
  }
  // Handle values like "/api" → graft onto current origin
  if (trimmed.startsWith('/')) {
    return `${window.location.origin}${trimmed}`;
  }
  // If missing protocol, default to http://
  if (!/^https?:\/\//i.test(trimmed)) {
    return `http://${trimmed}`;
  }
  return trimmed.replace(/\/$/, '');
};
const API_BASE = normalizeBase(rawBase);

function Upload({ onUploadSuccess }) {
  const [file, setFile] = useState(null);
  const [uploading, setUploading] = useState(false);
  const [annotating, setAnnotating] = useState(false);
  const [message, setMessage] = useState('');
  const [error, setError] = useState('');
  const [jobId, setJobId] = useState(null);

  const handleFileChange = (e) => {
    const selectedFile = e.target.files[0];
    setFile(selectedFile);
    setMessage('');
    setError('');
  };

  const handleUpload = async () => {
    if (!file) {
      setError('Please select a file first');
      return;
    }

    const formData = new FormData();
    formData.append('file', file);

    setUploading(true);
    setError('');
    setMessage('Uploading file...');

    try {
      const response = await axios.post(`${API_BASE}/api/upload`, formData, {
        headers: {
          'Content-Type': 'multipart/form-data',
          ...getAuthHeaders(),
        },
        onUploadProgress: (progressEvent) => {
          const percentCompleted = Math.round(
            (progressEvent.loaded * 100) / progressEvent.total
          );
          setMessage(`Uploading... ${percentCompleted}%`);
        },
      });

      const uploadedJobId = response.data.job_id;
      setJobId(uploadedJobId);
      setMessage('File uploaded successfully! Starting annotation...');
      
      // Automatically start annotation
      await handleAnnotate(uploadedJobId);
      
    } catch (err) {
      setError(`Upload failed: ${err.response?.data?.detail || err.message}`);
      setUploading(false);
    }
  };

  const handleAnnotate = async (uploadedJobId) => {
    setAnnotating(true);
    setMessage('Running annotation... This may take a few minutes.');

    try {
      const response = await axios.post(`${API_BASE}/api/annotate/${uploadedJobId}`, null, { headers: getAuthHeaders() });
      
      setMessage('Annotation completed successfully!');
      setAnnotating(false);
      setUploading(false);
      
      // Notify parent component
      onUploadSuccess(uploadedJobId);
      
    } catch (err) {
      setError(`Annotation failed: ${err.response?.data?.detail || err.message}`);
      setAnnotating(false);
      setUploading(false);
    }
  };

  const formatFileSize = (bytes) => {
    if (bytes === 0) return '0 Bytes';
    const k = 1024;
    const sizes = ['Bytes', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return Math.round(bytes / Math.pow(k, i) * 100) / 100 + ' ' + sizes[i];
  };

  return (
    <div className="upload-container">
      <div className="card">
        <h2>📤 Upload Your Data</h2>
        <p className="subtitle">
          Upload single-cell omics data for automated cell type annotation
        </p>

        <div className="supported-formats">
          <h3>Supported Formats:</h3>
          <ul>
            <li><strong>.h5ad</strong> - Scanpy AnnData HDF5 format</li>
            <li><strong>.h5</strong> - HDF5 format</li>
            <li><strong>.rds</strong> - Seurat RDS format</li>
          </ul>
        </div>

        <div className="upload-section">
          <input
            type="file"
            accept=".h5ad,.h5,.rds"
            onChange={handleFileChange}
            disabled={uploading}
            className="file-input"
            id="file-upload"
          />
          <label htmlFor="file-upload" className="file-label">
            {file ? '📄 ' + file.name : '📁 Choose File'}
          </label>
          
          {file && (
            <div className="file-info">
              <p>Size: {formatFileSize(file.size)}</p>
            </div>
          )}

          <button
            onClick={handleUpload}
            disabled={!file || uploading || annotating}
            className="btn btn-primary"
          >
            {uploading || annotating ? '⏳ Processing...' : '🚀 Upload & Annotate'}
          </button>
        </div>

        {message && !error && (
          <div className="message success">
            ✅ {message}
          </div>
        )}

        {error && (
          <div className="message error">
            ❌ {error}
          </div>
        )}

        {(uploading || annotating) && (
          <div className="progress-indicator">
            <div className="spinner"></div>
            <p>Please wait, do not close this window...</p>
          </div>
        )}

        <div className="info-box">
          <h3>📝 Important Notes:</h3>
          <ul>
            <li>Maximum file size: 5GB</li>
            <li>Processing time varies based on dataset size (typically 2-10 minutes)</li>
            <li>Ensure your data is preprocessed and normalized</li>
            <li>Results will be available in the Download tab</li>
          </ul>
        </div>
      </div>
    </div>
  );
}

export default Upload;