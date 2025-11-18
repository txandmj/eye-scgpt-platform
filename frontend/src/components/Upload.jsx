import React, { useState } from 'react';
import axios from 'axios';
import { getAuthHeaders } from '../auth';

const API_URL = process.env.REACT_APP_API_URL || 'http://localhost:8000';

function Upload({ onUploadSuccess }) {
  const [file, setFile] = useState(null);
  const [uploading, setUploading] = useState(false);
  const [annotating, setAnnotating] = useState(false);
  const [message, setMessage] = useState('');
  const [error, setError] = useState('');
  const [jobId, setJobId] = useState(null);
  const [uploadProgress, setUploadProgress] = useState(0);
  const [cancelSource, setCancelSource] = useState(null);

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

    setUploading(true);
    setError('');
    setMessage('Preparing upload...');
    setUploadProgress(0);

    // Create cancel token for upload
    const source = axios.CancelToken.source();
    setCancelSource(source);

    try {
      // Step 1: Create job and get job_id IMMEDIATELY (before upload starts)
      const prepareResponse = await axios.post(
        `${API_URL}/api/upload/prepare`,
        { filename: file.name },
        { headers: getAuthHeaders() }
      );
      
      const uploadedJobId = prepareResponse.data.job_id;
      setJobId(uploadedJobId);
      setMessage(`Uploading file... Job ID: ${uploadedJobId}`);

      // Step 2: Upload file with job_id
      const formData = new FormData();
      formData.append('file', file);
      formData.append('job_id', uploadedJobId);

      const uploadResponse = await axios.post(`${API_URL}/api/upload`, formData, {
        headers: {
          'Content-Type': 'multipart/form-data',
          ...getAuthHeaders(),
        },
        cancelToken: source.token,
        onUploadProgress: (progressEvent) => {
          const percentCompleted = Math.round(
            (progressEvent.loaded * 100) / progressEvent.total
          );
          setUploadProgress(percentCompleted);
          setMessage(`Uploading... ${percentCompleted}%`);
        },
      });

      setUploadProgress(100);
      setMessage('File uploaded successfully! Starting annotation...');
      setCancelSource(null); // Clear cancel source after upload completes
      
      // Automatically start annotation
      await handleAnnotate(uploadedJobId);
      
    } catch (err) {
      if (axios.isCancel(err)) {
        setMessage('Upload cancelled');
        setError('');
        // Cancel the job on backend
        if (jobId) {
          try {
            await axios.post(
              `${API_URL}/api/jobs/${jobId}/cancel`,
              null,
              { headers: getAuthHeaders() }
            );
          } catch (cancelErr) {
            console.error('Failed to cancel job:', cancelErr);
          }
        }
      } else {
        setError(`Upload failed: ${err.response?.data?.detail || err.message}`);
      }
      setUploading(false);
      setUploadProgress(0);
      setJobId(null);
      setCancelSource(null);
    }
  };

  const handleCancelUpload = () => {
    if (cancelSource) {
      cancelSource.cancel('Upload cancelled by user');
      setCancelSource(null);
    }
  };

  const handleAnnotate = async (uploadedJobId) => {
    setAnnotating(true);
    setMessage('Running annotation... This may take a few minutes.');

    try {
      const response = await axios.post(
        `${API_URL}/api/annotate/${uploadedJobId}`,
        null,
        { headers: getAuthHeaders() }
      );
      
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

        {(uploading || annotating) && jobId && (
          <div className="info-box" style={{ marginTop: 20 }}>
            <p><strong>📋 Job ID:</strong> <code>{jobId}</code></p>
            <p style={{ fontSize: '0.9em', color: '#666', marginTop: 5 }}>
              You can navigate to other pages (History, Download) - your job will continue processing in the background.
            </p>
            {uploading && (
              <button
                onClick={handleCancelUpload}
                className="btn btn-secondary"
                style={{ 
                  marginTop: 15,
                  backgroundColor: '#ef4444', 
                  color: 'white',
                  border: 'none'
                }}
              >
                🚫 Cancel Upload
              </button>
            )}
          </div>
        )}

        {uploading && uploadProgress > 0 && (
          <div style={{ marginTop: 20 }}>
            <div style={{ 
              width: '100%', 
              backgroundColor: '#e9ecef', 
              borderRadius: '8px', 
              height: '30px',
              overflow: 'hidden',
              position: 'relative'
            }}>
              <div style={{
                width: `${uploadProgress}%`,
                backgroundColor: '#667eea',
                height: '100%',
                transition: 'width 0.3s ease',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                color: 'white',
                fontSize: '0.85em',
                fontWeight: 600
              }}>
                {uploadProgress < 100 ? `${uploadProgress}%` : 'Upload Complete!'}
              </div>
            </div>
          </div>
        )}

        {(uploading || annotating) && (
          <div className="progress-indicator">
            <div className="spinner"></div>
            <p>{uploading ? 'Uploading file...' : 'Processing annotation...'}</p>
            <p style={{ fontSize: '0.9em', color: '#666', marginTop: 10 }}>
              You can close this page - the job will continue in the background. Check History page for updates.
            </p>
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