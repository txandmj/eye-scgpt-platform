// Simple auth header management using Google ID token payload

function decodeJwtPayload(token) {
  try {
    const [, payload] = token.split('.');
    const json = atob(payload.replace(/-/g, '+').replace(/_/g, '/'));
    return JSON.parse(decodeURIComponent(escape(json)));
  } catch {
    return null;
  }
}

export function setAuthFromGoogle(credential, fallback = {}) {
  const payload = decodeJwtPayload(credential) || {};
  const data = {
    google_sub: payload.sub || fallback.sub || 'test-sub',
    email: payload.email || fallback.email || 'test@example.com',
    name: payload.name || fallback.name || 'Test User',
  };
  try {
    localStorage.setItem('auth.google_sub', data.google_sub);
    localStorage.setItem('auth.email', data.email);
    localStorage.setItem('auth.name', data.name);
  } catch (_) {
    // ignore storage errors
  }
  return data;
}

export function getAuthHeaders() {
  try {
    const sub = localStorage.getItem('auth.google_sub');
    const email = localStorage.getItem('auth.email');
    const name = localStorage.getItem('auth.name');
    const headers = {};
    if (sub) headers['X-Google-Sub'] = sub;
    if (email) headers['X-User-Email'] = email;
    if (name) headers['X-User-Name'] = name;
    return headers;
  } catch {
    return {};
  }
}

export function isAuthenticated() {
  try {
    return !!localStorage.getItem('auth.google_sub');
  } catch {
    return false;
  }
}

export function clearAuth() {
  try {
    localStorage.removeItem('auth.google_sub');
    localStorage.removeItem('auth.email');
    localStorage.removeItem('auth.name');
  } catch {
    // ignore
  }
}

export async function logout() {
  // Try to log out from Google if the SDK is available (safe no-op otherwise)
  try {
    // Lazy import; if the package is not installed, it will throw and we ignore
    const mod = await import(/* webpackIgnore: true */ '@react-oauth/google').catch(() => null);
    if (mod && typeof mod.googleLogout === 'function') {
      try {
        mod.googleLogout();
      } catch (_) {
        // ignore runtime errors from SDK
      }
    }
  } catch (_) {
    // ignore import errors
  }
  // Always clear local auth
  clearAuth();
}


