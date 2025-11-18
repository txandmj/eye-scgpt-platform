import React, { useState, useEffect } from 'react';
import './App.css';
import Home from './components/Home';
import Upload from './components/Upload';
import Download from './components/Download';
import History from './components/History';
import Tutorial from './components/Tutorial';
import Login from './components/Login';
import Contact from './components/Contact';
import Navigation from './components/Navigation';
import { isAuthenticated, getUserInfo } from './auth';

function App() {
  const [activeTab, setActiveTab] = useState('home');
  const [currentJobId, setCurrentJobId] = useState(null);
  const [authenticated, setAuthenticated] = useState(isAuthenticated());
  const [userInfo, setUserInfo] = useState(getUserInfo());

  useEffect(() => {
    // Check auth status on mount and when tab changes
    const authStatus = isAuthenticated();
    setAuthenticated(authStatus);
    setUserInfo(getUserInfo());
    
    // Redirect to login if trying to access protected routes
    const protectedTabs = ['upload', 'download', 'history'];
    if (!authStatus && protectedTabs.includes(activeTab)) {
      setActiveTab('login');
    }
  }, [activeTab]);

  const handleUploadSuccess = (jobId) => {
    setCurrentJobId(jobId);
    setActiveTab('download');
  };

  const handleLoginSuccess = () => {
    setAuthenticated(true);
    setUserInfo(getUserInfo());
    // Redirect to upload after login
    setActiveTab('upload');
  };

  const handleLogout = () => {
    setAuthenticated(false);
    setUserInfo({ email: null, name: null });
    setActiveTab('home');
  };

  return (
    <div className="App">
      <header className="App-header">
        <div className="header-content">
          <div className="header-title">
            <h1>🔬 Eye-scGPT Annotation Platform</h1>
            <p>Automated Cell Type Annotation for Single-Cell Omics Data</p>
          </div>
          {authenticated && userInfo.email && (
            <div className="user-info">
              <span className="user-name">{userInfo.name || userInfo.email}</span>
              <span className="user-email">{userInfo.email}</span>
            </div>
          )}
        </div>
      </header>

      <Navigation 
        activeTab={activeTab} 
        setActiveTab={setActiveTab}
        authenticated={authenticated}
        onLogout={handleLogout}
      />

      <main className="App-main">
        {activeTab === 'home' && (
          <Home setActiveTab={setActiveTab} />
        )}
        
        {activeTab === 'upload' && (
          authenticated ? (
            <Upload onUploadSuccess={handleUploadSuccess} />
          ) : (
            <Login setActiveTab={setActiveTab} onSuccess={handleLoginSuccess} />
          )
        )}
        
        {activeTab === 'download' && (
          authenticated ? (
            <Download jobId={currentJobId} />
          ) : (
            <Login setActiveTab={setActiveTab} onSuccess={handleLoginSuccess} />
          )
        )}
        
        {activeTab === 'history' && (
          authenticated ? (
            <History setActiveTab={setActiveTab} />
          ) : (
            <Login setActiveTab={setActiveTab} onSuccess={handleLoginSuccess} />
          )
        )}
        
        {activeTab === 'tutorial' && (
          <Tutorial />
        )}
        
        {activeTab === 'login' && (
          <Login setActiveTab={setActiveTab} onSuccess={handleLoginSuccess} />
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