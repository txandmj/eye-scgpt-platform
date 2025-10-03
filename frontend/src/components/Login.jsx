import React, { useState } from 'react';

function Login({ setActiveTab }) {
  const [formData, setFormData] = useState({
    email: '',
    password: ''
  });

  const handleInputChange = (e) => {
    const { name, value } = e.target;
    setFormData(prev => ({
      ...prev,
      [name]: value
    }));
  };

  const handleSubmit = (e) => {
    e.preventDefault();
    // TODO: Implement login logic
    console.log('Login attempt:', formData);
    alert('Login functionality will be implemented later');
  };

  return (
    <div className="card">
      <h2>🔐 Login</h2>
      <p className="subtitle">Access your account to manage your annotations</p>
      
      <form className="login-form" onSubmit={handleSubmit}>
        <div className="input-group">
          <input 
            type="email" 
            name="email"
            className="text-input" 
            placeholder="Email address"
            value={formData.email}
            onChange={handleInputChange}
            required
          />
          <input 
            type="password" 
            name="password"
            className="text-input" 
            placeholder="Password"
            value={formData.password}
            onChange={handleInputChange}
            required
          />
        </div>
        <button type="submit" className="btn btn-primary">
          Login
        </button>
        <p style={{marginTop: '20px', color: '#666'}}>
          Don't have an account?{' '}
          <a href="#" onClick={(e) => {
            e.preventDefault();
            setActiveTab('account');
          }}>
            Sign up here
          </a>
        </p>
      </form>
    </div>
  );
}

export default Login;
