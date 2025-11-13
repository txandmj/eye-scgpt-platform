import React from 'react';
import { isAuthenticated, logout } from '../auth';

function Navigation({ activeTab, setActiveTab }) {
  return (
    <nav className="navigation">
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
        className={`nav-button ${activeTab === 'tutorial' ? 'active' : ''}`}
        onClick={() => setActiveTab('tutorial')}
      >
        📚 Tutorial
      </button>
      {/* Login and Account removed from global nav in Google-only mode */}
      <button
        className={`nav-button ${activeTab === 'history' ? 'active' : ''}`}
        onClick={() => setActiveTab('history')}
      >
        🕘 History
      </button>
      <button
        className={`nav-button ${activeTab === 'contact' ? 'active' : ''}`}
        onClick={() => setActiveTab('contact')}
      >
        📞 Contact Us
      </button>
      {isAuthenticated() && (
        <div style={{ marginLeft: 'auto' }}>
          <button
            className="nav-button"
            onClick={async () => {
              if (!window.confirm('Are you sure you want to log out?')) return;
              await logout();
              alert("You've been logged out.");
              window.location.hash = '#/login';
            }}
          >
            🚪 Logout
          </button>
        </div>
      )}
    </nav>
  );
}

export default Navigation;