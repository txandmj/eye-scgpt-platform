import React, { useState, useEffect, useRef } from 'react';

function Login({ setActiveTab }) {
  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [message, setMessage] = useState('');
  const googleBtnRef = useRef(null);

  const GOOGLE_CLIENT_ID = process.env.REACT_APP_GOOGLE_CLIENT_ID || '978070226745-cpevojet4d29ep33vm1vm88mr9st4osm.apps.googleusercontent.com';
  const API_BASE = process.env.REACT_APP_API_BASE || 'http://localhost:8000';

  useEffect(() => {
    let retryTimer;

    const initGoogle = () => {
      const google = window.google;
      if (!google || !google.accounts || !google.accounts.id) {
        retryTimer = setTimeout(initGoogle, 300);
        return;
      }

      // Reduce noisy SDK logs in console
      try {
        google.accounts.id.setLogLevel('warn');
      } catch (_) {
        // ignore if not supported
      }

      google.accounts.id.initialize({
        client_id: GOOGLE_CLIENT_ID,
        use_fedcm_for_prompt_api: false,
        callback: async (response) => {
          try {
            setLoading(true);
            setError('');
            setMessage('');

            const credential = response.credential;
            if (!credential) {
              console.error('[Google Sign-In] No credential returned. Likely origin mismatch.', {
                origin: window.location.origin,
              });
              setError(`Google sign-in failed: no credential returned. Add ${window.location.origin} to Authorized JavaScript origins in Google Cloud Console.`);
              return;
            }
            const res = await fetch(`${API_BASE}/google-login`, {
              method: 'POST',
              headers: { 'Content-Type': 'application/json' },
              body: JSON.stringify({ token: credential })
            });

            if (!res.ok) {
              let detail = 'Google login failed';
              try {
                const data = await res.json();
                detail = data?.detail || detail;
              } catch (_) {
                // ignore json parse error
              }
              throw new Error(detail);
            }

            const data = await res.json();
            setMessage(`Welcome ${data.name || data.email || 'user'}!`);
          } catch (err) {
            setError(err?.message || 'Google login failed. Please try again.');
          } finally {
            setLoading(false);
          }
        }
      });

      if (googleBtnRef.current) {
        google.accounts.id.renderButton(googleBtnRef.current, {
          type: 'standard',
          theme: 'outline',
          size: 'large',
          shape: 'rectangular',
          text: 'signin_with',
        });
      }

      // Expose client id for quick verification in console (read-only)
      try {
        Object.defineProperty(window, '__GOOGLE_CLIENT_ID', {
          value: GOOGLE_CLIENT_ID,
          writable: false,
          configurable: false,
          enumerable: false
        });
      } catch (_) {
        // ignore if already defined
      }
    };

    initGoogle();

    return () => {
      if (retryTimer) clearTimeout(retryTimer);
    };
  }, [GOOGLE_CLIENT_ID, setActiveTab]);

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

        <div style={{ marginBottom: '16px' }}>
          <div ref={googleBtnRef} style={{ display: 'inline-block' }} />
        </div>

        <div style={{ display: 'flex', alignItems: 'center', gap: 12, margin: '12px 0' }}>
          <div style={{ height: 1, background: '#eee', flex: 1 }} />
          <div style={{ color: '#999', fontSize: 12 }}>or</div>
          <div style={{ height: 1, background: '#eee', flex: 1 }} />
        </div>

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
              autoComplete="email"
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
              autoComplete="current-password"
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
