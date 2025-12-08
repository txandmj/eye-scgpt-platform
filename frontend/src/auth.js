// Simple auth header management using Google ID token payload
// Using sessionStorage instead of localStorage to support multiple users in different browser tabs/windows
// sessionStorage is isolated per tab, so each user can have their own session

function decodeJwtPayload(token) {
  try {
    const [, payload] = token.split('.');
    const json = atob(payload.replace(/-/g, '+').replace(/_/g, '/'));
    return JSON.parse(decodeURIComponent(escape(json)));
  } catch {
    return null;
  }
}

// Helper to get storage (sessionStorage for per-tab isolation)
function getStorage() {
  try {
    return window.sessionStorage;
  } catch {
    // Fallback to localStorage if sessionStorage is not available
    try {
      return window.localStorage;
    } catch {
      return null;
    }
  }
}

export function setAuthFromGoogle(credential, fallback = {}) {
  const payload = decodeJwtPayload(credential) || {};
  const data = {
    google_sub: payload.sub || fallback.sub || 'test-sub',
    email: payload.email || fallback.email || 'test@example.com',
    name: payload.name || fallback.name || 'Test User',
  };
  const storage = getStorage();
  if (storage) {
    try {
      storage.setItem('auth.google_sub', data.google_sub);
      storage.setItem('auth.email', data.email);
      storage.setItem('auth.name', data.name);
    } catch (_) {
      // ignore storage errors
    }
  }
  return data;
}

export function getAuthHeaders() {
  const storage = getStorage();
  if (!storage) return {};
  
  try {
    const sub = storage.getItem('auth.google_sub');
    const email = storage.getItem('auth.email');
    const name = storage.getItem('auth.name');
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
  const storage = getStorage();
  if (!storage) return false;
  
  try {
    return !!storage.getItem('auth.google_sub');
  } catch {
    return false;
  }
}

export function getUserInfo() {
  const storage = getStorage();
  if (!storage) return { email: null, name: null };
  
  try {
    const email = storage.getItem('auth.email');
    const name = storage.getItem('auth.name');
    return {
      email: email || null,
      name: name || null,
    };
  } catch {
    return { email: null, name: null };
  }
}

export function clearAuth() {
  const storage = getStorage();
  if (!storage) return;
  
  try {
    storage.removeItem('auth.google_sub');
    storage.removeItem('auth.email');
    storage.removeItem('auth.name');
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

