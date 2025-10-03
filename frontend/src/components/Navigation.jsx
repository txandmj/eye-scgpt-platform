import React from 'react';

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
      <button
        className={`nav-button ${activeTab === 'login' ? 'active' : ''}`}
        onClick={() => setActiveTab('login')}
      >
        🔐 Login
      </button>
      <button
        className={`nav-button ${activeTab === 'account' ? 'active' : ''}`}
        onClick={() => setActiveTab('account')}
      >
        👤 Account
      </button>
      <button
        className={`nav-button ${activeTab === 'contact' ? 'active' : ''}`}
        onClick={() => setActiveTab('contact')}
      >
        📞 Contact Us
      </button>
    </nav>
  );
}

export default Navigation;