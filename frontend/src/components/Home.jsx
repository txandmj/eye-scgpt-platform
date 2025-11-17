import React from 'react';
import { isAuthenticated } from '../auth';

function Home({ setActiveTab }) {
  const authenticated = isAuthenticated();

  return (
    <div className="card">
      <h2>🔬 Welcome to Eye-scGPT Annotation Platform</h2>
      <p className="subtitle">
        Automated Cell Type Annotation for Single-Cell Omics Data
      </p>

      <div className="home-content">
        <div className="feature-section">
          <h3>✨ What is Eye-scGPT?</h3>
          <p>
            Eye-scGPT is a powerful platform for automated cell type annotation 
            of single-cell omics data using state-of-the-art GPT-based models. 
            Upload your data, get accurate predictions, and visualize results 
            with beautiful UMAP plots.
          </p>
        </div>

        <div className="feature-section">
          <h3>🚀 Key Features</h3>
          <ul className="feature-list">
            <li>📊 <strong>Automated Annotation</strong> - Fast and accurate cell type predictions</li>
            <li>📈 <strong>UMAP Visualization</strong> - Beautiful visualizations of your data</li>
            <li>💾 <strong>Easy Download</strong> - Export results in CSV format</li>
            <li>📚 <strong>Comprehensive Tutorials</strong> - Learn how to use the platform</li>
            <li>🔐 <strong>Secure</strong> - Your data is protected with Google authentication</li>
          </ul>
        </div>

        <div className="feature-section">
          <h3>📖 How to Get Started</h3>
          <ol className="steps-list">
            <li>Click <strong>Login</strong> to sign in with your Google account</li>
            <li>Upload your <code>.h5ad</code> file using the Upload tab</li>
            <li>Wait for processing to complete (typically 5-15 minutes)</li>
            <li>Download your results and UMAP visualizations</li>
          </ol>
        </div>

        {!authenticated && (
          <div className="cta-section">
            <h3>🔐 Ready to Start?</h3>
            <p>Sign in with your Google account to begin annotating your single-cell data.</p>
            <button 
              className="btn btn-primary"
              onClick={() => setActiveTab('login')}
              style={{ marginTop: '20px', fontSize: '18px', padding: '12px 30px' }}
            >
              🔐 Sign In with Google
            </button>
          </div>
        )}

        {authenticated && (
          <div className="cta-section">
            <h3>✅ You're Signed In!</h3>
            <p>Ready to analyze your data? Navigate to Upload to get started.</p>
            <button 
              className="btn btn-primary"
              onClick={() => setActiveTab('upload')}
              style={{ marginTop: '20px', fontSize: '18px', padding: '12px 30px' }}
            >
              📤 Upload Your Data
            </button>
          </div>
        )}

        <div className="info-box" style={{ marginTop: '30px' }}>
          <h3>📝 Important Notes</h3>
          <ul>
            <li>Supported file format: <code>.h5ad</code> (Scanpy AnnData HDF5 format)</li>
            <li>Maximum file size: 5GB</li>
            <li>Processing time: 5-15 minutes for typical datasets</li>
            <li>For research purposes only</li>
          </ul>
        </div>
      </div>
    </div>
  );
}

export default Home;

