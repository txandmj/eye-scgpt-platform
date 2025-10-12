import React, { useState, useEffect } from 'react';
import axios from 'axios';

const API_URL = process.env.REACT_APP_API_URL || 'http://localhost:8000';

function Download({ jobId }) {
  const [results, setResults] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [inputJobId, setInputJobId] = useState('');

  useEffect(() => {
    if (jobId) {
      setInputJobId(jobId);
      fetchResults(jobId);
    }
  }, [jobId]);

  const fetchResults = async (id) => {
    if (!id) {
      setError('Please enter a Job ID');
      return;
    }

    setLoading(true);
    setError('');

    try {
      const response = await axios.get(`${API_URL}/api/job/${id}`);
      const job = response.data;
      
      if (job.status === 'pending' || job.status === 'processing') {
        setError('Job is still processing. Please wait for completion.');
        setResults([]);
      } else if (job.status === 'failed') {
        setError(`Job failed: ${job.error || 'Unknown error'}`);
        setResults([]);
      } else if (job.status === 'completed' && job.results) {
        // Fetch the actual results
        const downloadResponse = await axios.get(`${API_URL}/api/download/${id}`);
        setResults([{
          filename: downloadResponse.data.filename,
          content: downloadResponse.data.content,
          downloadUrl: downloadResponse.data.download_url
        }]);
      } else {
        setError('No results found for this Job ID.');
        setResults([]);
      }
    } catch (err) {
      setError(`Failed to fetch results: ${err.response?.data?.detail || err.message}`);
      setResults([]);
    } finally {
      setLoading(false);
    }
  };

  const handleDownload = (file) => {
    if (file.downloadUrl) {
      window.open(`${API_URL}${file.downloadUrl}`, '_blank');
    } else {
      // Fallback to direct file download
      const downloadUrl = `${API_URL}/api/download/${inputJobId}/file`;
      window.open(downloadUrl, '_blank');
    }
  };

  const formatFileSize = (bytes) => {
    if (bytes === 0) return '0 Bytes';
    const k = 1024;
    const sizes = ['Bytes', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return Math.round(bytes / Math.pow(k, i) * 100) / 100 + ' ' + sizes[i];
  };

  const getFileIcon = (filename) => {
    if (filename.endsWith('.csv')) return '📊';
    if (filename.endsWith('.xlsx')) return '📈';
    if (filename.endsWith('.h5ad')) return '🔬';
    if (filename.endsWith('.png') || filename.endsWith('.jpg')) return '🖼️';
    return '📄';
  };

  return (
    <div className="download-container">
      <div className="card">
        <h2>📥 Download Results</h2>
        <p className="subtitle">
          Retrieve your annotated single-cell data and visualizations
        </p>

        <div className="job-id-section">
          <label htmlFor="job-id-input">Enter Job ID:</label>
          <div className="input-group">
            <input
              id="job-id-input"
              type="text"
              value={inputJobId}
              onChange={(e) => setInputJobId(e.target.value)}
              placeholder="e.g., 20241002_143025_a1b2c3d4"
              className="text-input"
            />
            <button
              onClick={() => fetchResults(inputJobId)}
              disabled={loading || !inputJobId}
              className="btn btn-secondary"
            >
              {loading ? '⏳ Loading...' : '🔍 Fetch Results'}
            </button>
          </div>
        </div>

        {error && (
          <div className="message error">
            ❌ {error}
          </div>
        )}

        {results.length > 0 && (
          <div className="results-section">
            <h3>Available Downloads:</h3>
            <div className="results-list">
              {results.map((file, index) => (
                <div key={index} className="result-item">
                  <div className="file-details">
                    <span className="file-icon">{getFileIcon(file.filename)}</span>
                    <div className="file-name-size">
                      <strong>{file.filename}</strong>
                      <span className="file-size">Ready for download</span>
                    </div>
                  </div>
                  <button
                    onClick={() => handleDownload(file)}
                    className="btn btn-small"
                  >
                    ⬇️ Download
                  </button>
                </div>
              ))}
            </div>
          </div>
        )}

        {results.length === 0 && !loading && !error && (
          <div className="empty-state">
            <p>📭 No results to display</p>
            <p>Enter a Job ID to retrieve your annotation results</p>
          </div>
        )}

        <div className="info-box">
          <h3>📋 Result Files:</h3>
          <ul>
            <li><strong>*_annotations.csv</strong> - Cell type annotations in CSV format</li>
            <li><strong>*_annotations.xlsx</strong> - Cell type annotations in Excel format</li>
            <li><strong>*_annotated.h5ad</strong> - Annotated AnnData object with predictions</li>
            <li><strong>*_umap.png</strong> - UMAP visualization with cell type labels</li>
          </ul>
          
          <h3>💡 Tips:</h3>
          <ul>
            <li>Save your Job ID to retrieve results later</li>
            <li>Results are retained for 30 days</li>
            <li>Download all files for complete analysis</li>
          </ul>
        </div>
      </div>
    </div>
  );
}

export default Download;