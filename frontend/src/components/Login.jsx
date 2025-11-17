import React, { useEffect, useRef } from 'react';
import { setAuthFromGoogle } from '../auth';

function Login({ setActiveTab, onSuccess }) {
  const [error, setError] = React.useState('');
  const [message, setMessage] = React.useState('');
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
            // Persist auth headers for subsequent API calls
            try {
              setAuthFromGoogle(credential, { email: data.email, name: data.name });
            } catch (_) {
              // ignore store errors
            }
            if (onSuccess) onSuccess();
          } catch (err) {
            setError(err?.message || 'Google login failed. Please try again.');
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
      if (retryTimer) {
        clearTimeout(retryTimer);
      }
    };
  }, [GOOGLE_CLIENT_ID, setActiveTab, onSuccess]);

  return (
    <div className="login-container">
      <div className="card">
        <h2>🔐 Sign In</h2>
        <p className="subtitle">
          Sign in with your Google account to access upload, download, and history features
        </p>
        <div ref={googleBtnRef} style={{ display: 'inline-block', marginTop: 12 }} />
        <p style={{ color: '#777', fontSize: 12, marginTop: 8 }}>Use your institutional Google account</p>

        {message && !error && (
          <div className="message success" style={{ marginTop: 16 }}>
            ✅ {message}
          </div>
        )}

        {error && (
          <div className="message error" style={{ marginTop: 16 }}>
            ❌ {error}
          </div>
        )}

        <div className="info-box" style={{ marginTop: 24 }}>
          <h3>📝 Why Sign In?</h3>
          <ul>
            <li>Upload and process your single-cell data</li>
            <li>Download your results securely</li>
            <li>Track your annotation history</li>
            <li>Your data is associated with your account</li>
          </ul>
        </div>
      </div>
    </div>
  );
}

export default Login;
