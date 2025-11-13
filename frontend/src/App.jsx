import React, { useEffect, useRef, useState } from 'react';
import './App.css';
import Upload from './components/Upload';
import Download from './components/Download';
import Tutorial from './components/Tutorial';
import Login from './components/Login';
import Contact from './components/Contact';
import Navigation from './components/Navigation';
import HistoryPage from './components/HistoryPage';
import { isAuthenticated as getAuth } from './auth';

function App() {
  const [activeTab, setActiveTab] = useState('upload');
  const [currentJobId, setCurrentJobId] = useState(null);
  const [authed, setAuthed] = useState(getAuth());
  const lastAttemptRef = useRef(null);

  // Lightweight hash-based routing and guards
  useEffect(() => {
    const readHash = () => {
      const h = window.location.hash.replace(/^#\/?/, '');
      return h || 'upload';
    };
    const onHash = () => {
      const next = readHash();
      setActiveTab(next);
    };
    onHash();
    window.addEventListener('hashchange', onHash);
    return () => window.removeEventListener('hashchange', onHash);
  }, []);

  useEffect(() => {
    const ok = getAuth();
    setAuthed(ok);
    const protectedPages = ['upload', 'download', 'history'];
    if (!ok && protectedPages.includes(activeTab)) {
      lastAttemptRef.current = activeTab;
      if (activeTab !== 'login') {
        window.location.hash = '#/login';
        setActiveTab('login');
      }
    }
    if (ok && activeTab === 'login') {
      const dest = lastAttemptRef.current || 'upload';
      window.location.hash = `#/${dest}`;
      setActiveTab(dest);
      lastAttemptRef.current = null;
    }
  }, [activeTab]);

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

      {authed && activeTab !== 'login' && (
        <Navigation activeTab={activeTab} setActiveTab={(t)=>{ window.location.hash = `#/${t}`; setActiveTab(t); }} />
      )}

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
        {activeTab === 'history' && (
          <HistoryPage />
        )}
        
        {activeTab === 'login' && (
          <Login setActiveTab={(t)=>{ window.location.hash = `#/${t}`; setActiveTab(t); }} onSuccess={()=>{ setAuthed(true); const dest = lastAttemptRef.current || 'upload'; window.location.hash = `#/${dest}`; setActiveTab(dest); }} />
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