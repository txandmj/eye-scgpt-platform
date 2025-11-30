import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { getAuthHeaders } from '../auth';

const API_URL = 'http://eye.som.uci.edu:8100';

function History({ setActiveTab }) {
  const [jobs, setJobs] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState('');
  const [message, setMessage] = useState('');
  const [selectedJob, setSelectedJob] = useState(null);
  const [filterStatus, setFilterStatus] = useState('all'); // all, pending, processing, completed, failed
  const [sortBy, setSortBy] = useState('date'); // date, status, filename

  useEffect(() => {
    fetchHistory();
  }, [filterStatus, sortBy]);

  const handleCancelJob = async (jobId) => {
    if (!window.confirm('Are you sure you want to cancel this job? This action cannot be undone.')) {
      return;
    }

    try {
      const response = await axios.post(
        `${API_URL}/api/jobs/${jobId}/cancel`,
        null,
        { headers: getAuthHeaders() }
      );
      
      const successMsg = `Job ${jobId.substring(0, 8)}... cancelled successfully`;
      setMessage(successMsg);
      setError('');
      
      // Refresh history to show updated status
      setTimeout(() => {
        fetchHistory();
        setSelectedJob(null);
        setMessage(''); // Clear message after refresh
      }, 1000);
    } catch (err) {
      setError(`Failed to cancel job: ${err.response?.data?.detail || err.message}`);
    }
  };

  // Auto-refresh effect for active jobs
  useEffect(() => {
    // Auto-refresh every 30 seconds if there are active jobs (uploading, processing, or pending)
    const hasActiveJobs = jobs.some(job => 
      job.status === 'uploading' || job.status === 'processing' || job.status === 'pending'
    );
    
    let intervalId;
    if (hasActiveJobs && !loading) {
      intervalId = setInterval(() => {
        fetchHistory();
      }, 30000); // Refresh every 30 seconds
    }
    
    return () => {
      if (intervalId) {
        clearInterval(intervalId);
      }
    };
  }, [jobs, loading]); // Re-run when jobs change or loading state changes

  const fetchHistory = async () => {
    setLoading(true);
    setError('');
    
    try {
      const response = await axios.get(`${API_URL}/api/jobs`, {
        headers: getAuthHeaders()
      });
      
      let fetchedJobs = response.data.jobs || [];
      
      // Filter by status
      if (filterStatus !== 'all') {
        fetchedJobs = fetchedJobs.filter(job => job.status === filterStatus);
      }
      
      // Sort jobs
      fetchedJobs = [...fetchedJobs].sort((a, b) => {
        if (sortBy === 'date') {
          const dateA = new Date(a.upload_time || 0);
          const dateB = new Date(b.upload_time || 0);
          return dateB - dateA; // Newest first
        } else if (sortBy === 'status') {
          const statusOrder = { 
            uploading: 0,  // Show uploading jobs first
            processing: 1, 
            pending: 2, 
            completed: 3, 
            cancelled: 4,
            failed: 5 
          };
          return (statusOrder[a.status] || 99) - (statusOrder[b.status] || 99);
        } else if (sortBy === 'filename') {
          return (a.filename || '').localeCompare(b.filename || '');
        }
        return 0;
      });
      
      setJobs(fetchedJobs);
    } catch (err) {
      setError(`Failed to load history: ${err.response?.data?.detail || err.message}`);
      setJobs([]);
    } finally {
      setLoading(false);
    }
  };

  const handleJobClick = async (job) => {
    setSelectedJob(job);
    
    // If job is completed, fetch detailed results
    if (job.status === 'completed' && job.results) {
      try {
        const response = await axios.get(`${API_URL}/api/download/${job.id}`, {
          headers: getAuthHeaders()
        });
        setSelectedJob({
          ...job,
          downloadData: response.data
        });
      } catch (err) {
        console.error('Failed to fetch job details:', err);
      }
    }
  };

  const handleDownload = async (jobId, fileType = 'csv') => {
    try {
      let url;
      if (fileType === 'csv') {
        url = `${API_URL}/api/download/${jobId}/file`;
      } else if (fileType === 'umap') {
        url = `${API_URL}/api/download/${jobId}/umap`;
      }
      
      const response = await axios.get(url, {
        responseType: 'blob',
        headers: getAuthHeaders()
      });
      
      // Create blob URL and trigger download
      const blob = new Blob([response.data]);
      const blobUrl = window.URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = blobUrl;
      
      if (fileType === 'csv') {
        link.download = 'predictions.csv';
      } else {
        link.download = `umap_results_${jobId}.zip`;
      }
      
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      window.URL.revokeObjectURL(blobUrl);
    } catch (err) {
      setError(`Download failed: ${err.response?.data?.detail || err.message}`);
    }
  };

  const formatDate = (dateString) => {
    if (!dateString) return 'N/A';
    try {
      const date = new Date(dateString);
      return date.toLocaleString('en-US', {
        year: 'numeric',
        month: 'short',
        day: 'numeric',
        hour: '2-digit',
        minute: '2-digit'
      });
    } catch {
      return dateString;
    }
  };

  const getStatusBadge = (status) => {
    const badges = {
      completed: { emoji: '✅', color: '#10b981', label: 'Completed' },
      processing: { emoji: '⏳', color: '#3b82f6', label: 'Processing' },
      uploading: { emoji: '📤', color: '#8b5cf6', label: 'Uploading' },
      pending: { emoji: '⏸️', color: '#f59e0b', label: 'Pending' },
      failed: { emoji: '❌', color: '#ef4444', label: 'Failed' },
      cancelled: { emoji: '🚫', color: '#6b7280', label: 'Cancelled' }
    };
    const badge = badges[status] || badges.pending;
    return (
      <span style={{ 
        color: badge.color, 
        fontWeight: 600,
        display: 'inline-flex',
        alignItems: 'center',
        gap: '4px'
      }}>
        {badge.emoji} {badge.label}
      </span>
    );
  };

  const formatFileSize = (bytes) => {
    if (!bytes) return 'N/A';
    const k = 1024;
    const sizes = ['Bytes', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return Math.round(bytes / Math.pow(k, i) * 100) / 100 + ' ' + sizes[i];
  };

  const getDuration = (startTime, endTime) => {
    if (!startTime) return 'N/A';
    const start = new Date(startTime);
    const end = endTime ? new Date(endTime) : new Date();
    const diffMs = end - start;
    const diffMins = Math.floor(diffMs / 60000);
    const diffHours = Math.floor(diffMins / 60);
    
    if (diffHours > 0) {
      return `${diffHours}h ${diffMins % 60}m`;
    }
    return `${diffMins}m`;
  };

  return (
    <div className="history-container">
      <div className="card">
        <h2>📊 Analysis History</h2>
        <p className="subtitle">
          Review your past analyses and retrieve previous work
        </p>
        
        <div className="info-box" style={{ marginTop: 20, marginBottom: 20 }}>
          <p>
            <strong>💡 Tip:</strong> You can close this webpage at any time. Your analysis will continue running in the background. 
            When you return, check this page to see updated job statuses. Processing jobs auto-refresh every 30 seconds.
          </p>
        </div>

        {/* Filters and Sort */}
        <div className="history-controls">
          <div className="filter-group">
            <label htmlFor="status-filter">Filter by Status:</label>
            <select
              id="status-filter"
              value={filterStatus}
              onChange={(e) => setFilterStatus(e.target.value)}
              className="select-input"
            >
              <option value="all">All Jobs</option>
              <option value="completed">✅ Completed</option>
              <option value="processing">⏳ Processing</option>
              <option value="uploading">📤 Uploading</option>
              <option value="pending">⏸️ Pending</option>
              <option value="failed">❌ Failed</option>
              <option value="cancelled">🚫 Cancelled</option>
            </select>
          </div>

          <div className="filter-group">
            <label htmlFor="sort-filter">Sort by:</label>
            <select
              id="sort-filter"
              value={sortBy}
              onChange={(e) => setSortBy(e.target.value)}
              className="select-input"
            >
              <option value="date">Date (Newest First)</option>
              <option value="status">Status</option>
              <option value="filename">Filename</option>
            </select>
          </div>

          <button
            onClick={fetchHistory}
            className="btn btn-secondary"
            style={{ marginLeft: 'auto' }}
          >
            🔄 Refresh
          </button>
        </div>

        {error && (
          <div className="message error" style={{ marginTop: 20 }}>
            ❌ {error}
          </div>
        )}

        {message && !error && (
          <div className="message success" style={{ marginTop: 20 }}>
            ✅ {message}
          </div>
        )}

        {loading ? (
          <div className="loading-container" style={{ marginTop: 40 }}>
            <div className="spinner"></div>
            <p>Loading history...</p>
          </div>
        ) : jobs.length === 0 ? (
          <div className="empty-state" style={{ marginTop: 40 }}>
            <p>📭 No jobs found</p>
            <p>
              {filterStatus !== 'all'
                ? `No jobs with status "${filterStatus}"`
                : 'You haven\'t uploaded any files yet'}
            </p>
            {filterStatus !== 'all' && (
              <button
                className="btn btn-primary"
                onClick={() => setFilterStatus('all')}
                style={{ marginTop: 10 }}
              >
                Show All Jobs
              </button>
            )}
            {jobs.length === 0 && filterStatus === 'all' && (
              <button
                className="btn btn-primary"
                onClick={() => setActiveTab('upload')}
                style={{ marginTop: 10 }}
              >
                📤 Upload Your First File
              </button>
            )}
          </div>
        ) : (
          <div className="history-content">
            {/* Jobs List */}
            <div className="jobs-list">
              {jobs.map((job) => (
                <div
                  key={job.id}
                  className={`job-item ${selectedJob?.id === job.id ? 'selected' : ''}`}
                  onClick={() => handleJobClick(job)}
                >
                  <div className="job-item-header">
                    <div className="job-item-title">
                      <span className="job-id">{job.id.substring(0, 8)}...</span>
                      <span className="job-filename">{job.filename || 'Unknown file'}</span>
                    </div>
                    {getStatusBadge(job.status)}
                  </div>
                  
                  <div className="job-item-details">
                    <div className="job-detail">
                      <span className="detail-label">Uploaded:</span>
                      <span>{formatDate(job.upload_time)}</span>
                    </div>
                    {job.start_time && (
                      <div className="job-detail">
                        <span className="detail-label">Started:</span>
                        <span>{formatDate(job.start_time)}</span>
                      </div>
                    )}
                    {job.end_time && (
                      <div className="job-detail">
                        <span className="detail-label">Completed:</span>
                        <span>{formatDate(job.end_time)}</span>
                      </div>
                    )}
                    {job.start_time && (
                      <div className="job-detail">
                        <span className="detail-label">Duration:</span>
                        <span>{getDuration(job.start_time, job.end_time)}</span>
                      </div>
                    )}
                  </div>

                  {job.status === 'failed' && job.error && (
                    <div className="job-error">
                      <strong>Error:</strong> {job.error}
                    </div>
                  )}

                  {job.status === 'cancelled' && (
                    <div className="job-error" style={{ background: '#f3f4f6', borderLeftColor: '#6b7280' }}>
                      <strong>Cancelled:</strong> This job was cancelled by the user
                    </div>
                  )}

                  {/* Cancel button for active jobs */}
                  {(job.status === 'uploading' || job.status === 'processing' || job.status === 'pending') && (
                    <div style={{ marginTop: 15 }}>
                      <button
                        className="btn btn-secondary"
                        onClick={(e) => {
                          e.stopPropagation();
                          handleCancelJob(job.id);
                        }}
                        style={{ 
                          backgroundColor: '#ef4444', 
                          color: 'white',
                          border: 'none',
                          fontSize: '0.9em',
                          padding: '8px 16px'
                        }}
                      >
                        🚫 Cancel Job
                      </button>
                    </div>
                  )}
                </div>
              ))}
            </div>

            {/* Job Details Panel */}
            {selectedJob && (
              <div className="job-details-panel">
                <div className="panel-header">
                  <h3>Job Details</h3>
                  <button
                    className="close-button"
                    onClick={() => setSelectedJob(null)}
                  >
                    ✕
                  </button>
                </div>

                <div className="panel-content">
                  <div className="detail-section">
                    <h4>Basic Information</h4>
                    <div className="detail-grid">
                      <div className="detail-item">
                        <strong>Job ID:</strong>
                        <code>{selectedJob.id}</code>
                      </div>
                      <div className="detail-item">
                        <strong>Filename:</strong>
                        <span>{selectedJob.filename || 'N/A'}</span>
                      </div>
                      <div className="detail-item">
                        <strong>Status:</strong>
                        {getStatusBadge(selectedJob.status)}
                      </div>
                      <div className="detail-item">
                        <strong>Upload Time:</strong>
                        <span>{formatDate(selectedJob.upload_time)}</span>
                      </div>
                      {selectedJob.start_time && (
                        <div className="detail-item">
                          <strong>Start Time:</strong>
                          <span>{formatDate(selectedJob.start_time)}</span>
                        </div>
                      )}
                      {selectedJob.end_time && (
                        <div className="detail-item">
                          <strong>End Time:</strong>
                          <span>{formatDate(selectedJob.end_time)}</span>
                        </div>
                      )}
                      {selectedJob.start_time && (
                        <div className="detail-item">
                          <strong>Duration:</strong>
                          <span>{getDuration(selectedJob.start_time, selectedJob.end_time)}</span>
                        </div>
                      )}
                    </div>
                  </div>

                  {selectedJob.status === 'failed' && selectedJob.error && (
                    <div className="detail-section">
                      <h4>Error Information</h4>
                      <div className="error-box">
                        <pre>{selectedJob.error}</pre>
                      </div>
                    </div>
                  )}

                  {selectedJob.status === 'completed' && (
                    <div className="detail-section">
                      <h4>Results & Downloads</h4>
                      <div className="download-buttons">
                        <button
                          className="btn btn-primary"
                          onClick={() => handleDownload(selectedJob.id, 'csv')}
                        >
                          📊 Download Predictions CSV
                        </button>
                        
                        {selectedJob.downloadData?.umap_available && (
                          <button
                            className="btn btn-primary"
                            onClick={() => handleDownload(selectedJob.id, 'umap')}
                          >
                            📦 Download UMAP Results
                            {selectedJob.downloadData.num_plots && (
                              <span style={{ fontSize: '0.9em', opacity: 0.8 }}>
                                {' '}({selectedJob.downloadData.num_plots} plots)
                              </span>
                            )}
                          </button>
                        )}

                        {selectedJob.downloadData?.umap_preview_url && (
                          <div className="umap-preview-small" style={{ marginTop: 15 }}>
                            <h5>UMAP Preview</h5>
                            <img
                              src={`${API_URL}${selectedJob.downloadData.umap_preview_url}`}
                              alt="UMAP preview"
                              style={{ maxWidth: '100%', borderRadius: '8px' }}
                              onError={(e) => {
                                // Try to load with auth headers
                                const img = e.target;
                                axios.get(`${API_URL}${selectedJob.downloadData.umap_preview_url}`, {
                                  responseType: 'blob',
                                  headers: getAuthHeaders()
                                }).then(response => {
                                  const blobUrl = window.URL.createObjectURL(response.data);
                                  img.src = blobUrl;
                                }).catch(() => {
                                  img.style.display = 'none';
                                });
                              }}
                            />
                          </div>
                        )}

                        <button
                          className="btn btn-secondary"
                          onClick={() => {
                            setActiveTab('download');
                            // Pass job ID to download tab
                            window.history.pushState({ jobId: selectedJob.id }, '');
                          }}
                          style={{ marginTop: 10 }}
                        >
                          📥 Open in Download Tab
                        </button>
                      </div>
                    </div>
                  )}

                  {selectedJob.status === 'uploading' && (
                    <div className="detail-section">
                      <div className="info-box">
                        <p>📤 This job is currently uploading. Please wait for the upload to complete.</p>
                        <p>Once uploaded, the annotation process will start automatically.</p>
                        <button
                          className="btn btn-secondary"
                          onClick={() => handleCancelJob(selectedJob.id)}
                          style={{ 
                            marginTop: 15,
                            backgroundColor: '#ef4444', 
                            color: 'white',
                            border: 'none'
                          }}
                        >
                          🚫 Cancel Job
                        </button>
                      </div>
                    </div>
                  )}

                  {selectedJob.status === 'processing' && (
                    <div className="detail-section">
                      <div className="info-box">
                        <p>⏳ This job is currently processing. Please check back later.</p>
                        <p>Processing typically takes 5-15 minutes depending on dataset size.</p>
                        <button
                          className="btn btn-secondary"
                          onClick={() => handleCancelJob(selectedJob.id)}
                          style={{ 
                            marginTop: 15,
                            backgroundColor: '#ef4444', 
                            color: 'white',
                            border: 'none'
                          }}
                        >
                          🚫 Cancel Job
                        </button>
                      </div>
                    </div>
                  )}

                  {selectedJob.status === 'pending' && (
                    <div className="detail-section">
                      <div className="info-box">
                        <p>⏸️ This job is pending. It will start processing shortly.</p>
                        <button
                          className="btn btn-secondary"
                          onClick={() => handleCancelJob(selectedJob.id)}
                          style={{ 
                            marginTop: 15,
                            backgroundColor: '#ef4444', 
                            color: 'white',
                            border: 'none'
                          }}
                        >
                          🚫 Cancel Job
                        </button>
                      </div>
                    </div>
                  )}

                  {selectedJob.status === 'cancelled' && (
                    <div className="detail-section">
                      <div className="info-box" style={{ background: '#f3f4f6', borderLeftColor: '#6b7280' }}>
                        <p>🚫 This job was cancelled by the user.</p>
                        {selectedJob.error && (
                          <p style={{ marginTop: 10, fontStyle: 'italic' }}>{selectedJob.error}</p>
                        )}
                      </div>
                    </div>
                  )}
                </div>
              </div>
            )}
          </div>
        )}

        <div className="info-box" style={{ marginTop: 30 }}>
          <h3>ℹ️ Tips</h3>
          <ul>
            <li>Click on any job to view detailed information and download results</li>
            <li>Use filters to find specific jobs by status</li>
            <li>Completed jobs can download both predictions CSV and UMAP visualizations</li>
            <li>All jobs are stored and accessible as long as your account is active</li>
          </ul>
        </div>
      </div>
    </div>
  );
}

export default History;

