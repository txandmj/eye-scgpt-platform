import React, { useState } from 'react';

function Login({ setActiveTab }) {
  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [message, setMessage] = useState('');

  const handleSubmit = async (e) => {
    e.preventDefault();
    setLoading(true);
    setError('');
    setMessage('');

    try {
      // For now, this is a placeholder login
      // In production, you would implement actual authentication
      if (email && password) {
        setMessage('Login successful! (Demo mode)');
        setTimeout(() => {
          setActiveTab('account');
        }, 1000);
      } else {
        setError('Please enter both email and password');
      }
    } catch (err) {
      setError('Login failed. Please try again.');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="login-container">
      <div className="card">
        <h2>🔐 Login</h2>
        <p className="subtitle">
          Access your account to manage your annotation jobs
        </p>

        <form onSubmit={handleSubmit} className="login-form">
          <div className="form-group">
            <label htmlFor="email">Email:</label>
            <input
              id="email"
              type="email"
              value={email}
              onChange={(e) => setEmail(e.target.value)}
              placeholder="your.email@example.com"
              className="text-input"
              required
            />
          </div>

          <div className="form-group">
            <label htmlFor="password">Password:</label>
            <input
              id="password"
              type="password"
              value={password}
              onChange={(e) => setPassword(e.target.value)}
              placeholder="Enter your password"
              className="text-input"
              required
            />
          </div>

          <button
            type="submit"
            disabled={loading}
            className="btn btn-primary"
          >
            {loading ? '⏳ Logging in...' : '🚀 Login'}
          </button>
        </form>

        {message && !error && (
          <div className="message success">
            ✅ {message}
          </div>
        )}

        {error && (
          <div className="message error">
            ❌ {error}
          </div>
        )}

        <div className="login-links">
          <p>
            Don't have an account?{' '}
            <button 
              className="link-button" 
              onClick={() => setActiveTab('account')}
            >
              Create Account
            </button>
          </p>
          <p>
            <button className="link-button">Forgot Password?</button>
          </p>
        </div>

        <div className="info-box">
          <h3>📝 Demo Mode:</h3>
          <p>This is a demonstration version. In production, this would connect to a secure authentication system.</p>
          <p>You can use any email and password to "login" for demo purposes.</p>
        </div>
      </div>
    </div>
  );
}

export default Login;
