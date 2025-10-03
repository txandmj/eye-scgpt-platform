import React, { useState } from 'react';
import './App.css';
import Upload from './components/Upload';
import Download from './components/Download';
import Tutorial from './components/Tutorial';
import Login from './components/Login';
import Account from './components/Account';
import Contact from './components/Contact';
import Navigation from './components/Navigation';

function App() {
  const [activeTab, setActiveTab] = useState('upload');
  const [currentJobId, setCurrentJobId] = useState(null);

  const handleUploadSuccess = (jobId) => {
    setCurrentJobId(jobId);
    setActiveTab('download');
  };

  return (
    <div className="App">
      <header className="App-header">
        <h1>🔬 Eye-scGPT Annotation Platform</h1>
        <p>Automated Cell Type Annotation for Single-Cell Omics Data</p>
      </header>

      <Navigation activeTab={activeTab} setActiveTab={setActiveTab} />

      <main className="App-main">
        {activeTab === 'upload' && (
          <Upload onUploadSuccess={handleUploadSuccess} />
        )}
        
        {activeTab === 'download' && (
          <Download jobId={currentJobId} />
        )}
        
        {activeTab === 'tutorial' && (
          <Tutorial />
        )}
        
        {activeTab === 'login' && (
          <Login setActiveTab={setActiveTab} />
        )}
        
        {activeTab === 'account' && (
          <Account setActiveTab={setActiveTab} />
        )}
        
        {activeTab === 'contact' && (
          <Contact />
        )}
      </main>

      <footer className="App-footer">
        <p>Eye-scGPT Annotation Platform v1.0.0 | Powered by scGPT</p>
        <p>For research use only</p>
      </footer>
    </div>
  );
}

export default App;