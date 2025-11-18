import React from 'react';
import { logout } from '../auth';

function Navigation({ activeTab, setActiveTab, authenticated, onLogout }) {
  const handleLogout = async () => {
    await logout();
    if (onLogout) {
      onLogout();
    }
  };

  return (
    <nav className="navigation">
      <button
        className={`nav-button ${activeTab === 'home' ? 'active' : ''}`}
        onClick={() => setActiveTab('home')}
      >
        🏠 Home
      </button>
      <button
        className={`nav-button ${activeTab === 'tutorial' ? 'active' : ''}`}
        onClick={() => setActiveTab('tutorial')}
      >
        📚 Tutorial
      </button>
      {authenticated ? (
        <>
          <button
            className={`nav-button ${activeTab === 'upload' ? 'active' : ''}`}
            onClick={() => setActiveTab('upload')}
          >
            📤 Upload
          </button>
          <button
            className={`nav-button ${activeTab === 'download' ? 'active' : ''}`}
            onClick={() => setActiveTab('download')}
          >
            📥 Download
          </button>
          <button
            className={`nav-button ${activeTab === 'history' ? 'active' : ''}`}
            onClick={() => setActiveTab('history')}
          >
            📊 History
          </button>
          <button
            className="nav-button"
            onClick={handleLogout}
            style={{ color: '#ff6b6b' }}
          >
            🚪 Logout
          </button>
        </>
      ) : (
        <button
          className={`nav-button ${activeTab === 'login' ? 'active' : ''}`}
          onClick={() => setActiveTab('login')}
        >
          🔐 Sign In
        </button>
      )}
      <button
        className={`nav-button ${activeTab === 'contact' ? 'active' : ''}`}
        onClick={() => setActiveTab('contact')}
      >
        📞 Contact
      </button>
    </nav>
  );
}

export default Navigation;