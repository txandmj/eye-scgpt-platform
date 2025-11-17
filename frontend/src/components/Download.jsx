import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { getAuthHeaders } from '../auth';

const API_URL = process.env.REACT_APP_API_URL || 'http://localhost:8000';

function Download({ jobId }) {
  const [results, setResults] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [inputJobId, setInputJobId] = useState('');
  const [umapPreviewUrl, setUmapPreviewUrl] = useState(null);
  const [umapPreviewFilename, setUmapPreviewFilename] = useState('');

  useEffect(() => {
    if (jobId) {
      setInputJobId(jobId);
      fetchResults(jobId);
    }
  }, [jobId]);

  const fetchResults = async (id) => {
    if (!id) {
      setError('Please enter a Job ID');
      setUmapPreviewUrl(null);
      setUmapPreviewFilename('');
      return;
    }

    setLoading(true);
    setError('');
    setUmapPreviewUrl(null);
    setUmapPreviewFilename('');

    try {
      const response = await axios.get(`${API_URL}/api/job/${id}`, { headers: getAuthHeaders() });
      const job = response.data;
      
      if (job.status === 'pending' || job.status === 'processing') {
        setError('Job is still processing. Please wait for completion.');
        setResults([]);
      } else if (job.status === 'failed') {
        setError(`Job failed: ${job.error || 'Unknown error'}`);
        setResults([]);
      } else if (job.status === 'completed' && job.results) {
        // Fetch the actual results
        const downloadResponse = await axios.get(`${API_URL}/api/download/${id}`, { headers: getAuthHeaders() });
        const downloadData = downloadResponse.data;
        
        const availableResults = [{
          filename: downloadData.filename,
          content: downloadData.content,
          downloadUrl: downloadData.download_url,
          type: 'csv',
          description: 'Eye-scGPT predictions'
        }];
        
        // Add UMAP results if available
        if (downloadData.umap_available) {
          availableResults.push({
            filename: `umap_results_${id}.zip`,
            downloadUrl: downloadData.umap_download_url,
            type: 'zip',
            description: `UMAP visualizations (${downloadData.num_plots} plots)`,
            numPlots: downloadData.num_plots
          });

          if (downloadData.umap_preview_url) {
            // Fetch preview as blob with auth headers and create object URL
            try {
              const previewResponse = await axios.get(
                `${API_URL}${downloadData.umap_preview_url}`,
                { responseType: 'blob', headers: getAuthHeaders() }
              );
              const blobUrl = window.URL.createObjectURL(previewResponse.data);
              setUmapPreviewUrl(blobUrl);
              setUmapPreviewFilename(downloadData.umap_preview_filename || 'UMAP predictions preview');
            } catch (err) {
              console.error('Failed to load UMAP preview:', err);
              // Preview not critical, continue without it
            }
          }
        }
        
        setResults(availableResults);
      } else {
        setError('No results found for this Job ID.');
        setResults([]);
      }
    } catch (err) {
      setError(`Failed to fetch results: ${err.response?.data?.detail || err.message}`);
      setResults([]);
      setUmapPreviewUrl(null);
      setUmapPreviewFilename('');
    } finally {
      setLoading(false);
    }
  };

  const handleDownload = async (file) => {
    try {
      let downloadUrl = file.downloadUrl || `/api/download/${inputJobId}/file`;
      if (!downloadUrl.startsWith('http')) {
        downloadUrl = `${API_URL}${downloadUrl}`;
      }
      
      // Use axios to download with auth headers
      const response = await axios.get(downloadUrl, {
        responseType: 'blob',
        headers: getAuthHeaders()
      });
      
      // Create a blob URL and trigger download
      const blob = new Blob([response.data]);
      const url = window.URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = file.filename || 'download';
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      window.URL.revokeObjectURL(url);
    } catch (err) {
      setError(`Download failed: ${err.response?.data?.detail || err.message}`);
    }
  };

  const formatFileSize = (bytes) => {
    if (bytes === 0) return '0 Bytes';
    const k = 1024;
    const sizes = ['Bytes', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return Math.round(bytes / Math.pow(k, i) * 100) / 100 + ' ' + sizes[i];
  };

  const getFileIcon = (filename, type) => {
    if (type === 'csv') return '📊';
    if (type === 'zip') return '📦';
    if (filename.endsWith('.csv')) return '📊';
    if (filename.endsWith('.xlsx')) return '📈';
    if (filename.endsWith('.h5ad')) return '🔬';
    if (filename.endsWith('.png') || filename.endsWith('.jpg')) return '🖼️';
    if (filename.endsWith('.zip')) return '📦';
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

        {umapPreviewUrl && (
          <div className="umap-preview">
            <h3>UMAP Preview (Predictions)</h3>
            <img
              src={umapPreviewUrl}
              alt="UMAP predictions visualization"
              loading="lazy"
            />
            <p className="preview-caption">
              {umapPreviewFilename} — download the UMAP archive for the complete collection.
            </p>
          </div>
        )}

        {results.length > 0 && (
          <div className="results-section">
            <h3>Available Downloads:</h3>
            <div className="results-list">
              {results.map((file, index) => (
                <div key={index} className="result-item">
                  <div className="file-details">
                    <span className="file-icon">{getFileIcon(file.filename, file.type)}</span>
                    <div className="file-name-size">
                      <strong>{file.filename}</strong>
                      <span className="file-description">{file.description}</span>
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
            <li><strong>predictions.csv</strong> - Eye-scGPT cell type predictions in CSV format</li>
            <li><strong>umap_results_*.zip</strong> - UMAP visualizations including:
              <ul>
                <li>Enhanced h5ad file with merged metadata</li>
                <li>Multiple PNG plots for each metadata column</li>
                <li>Different visualization styles (wolabel, ondata, fontline)</li>
              </ul>
            </li>
          </ul>
          
          <h3>💡 Tips:</h3>
          <ul>
            <li>Save your Job ID to retrieve results later</li>
            <li>Download both CSV and UMAP files for complete analysis</li>
            <li>UMAP plots show cell type predictions on 2D projections</li>
            <li>Results are retained for 30 days</li>
          </ul>
        </div>
      </div>
    </div>
  );
}

export default Download;