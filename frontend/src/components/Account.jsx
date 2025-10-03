import React, { useState } from 'react';

function Account({ setActiveTab }) {
  const [formData, setFormData] = useState({
    fullName: '',
    email: '',
    password: '',
    confirmPassword: ''
  });

  const [errors, setErrors] = useState({});

  const handleInputChange = (e) => {
    const { name, value } = e.target;
    setFormData(prev => ({
      ...prev,
      [name]: value
    }));
    
    // Clear error when user starts typing
    if (errors[name]) {
      setErrors(prev => ({
        ...prev,
        [name]: ''
      }));
    }
  };

  const validateForm = () => {
    const newErrors = {};

    if (!formData.fullName.trim()) {
      newErrors.fullName = 'Full name is required';
    }

    if (!formData.email.trim()) {
      newErrors.email = 'Email is required';
    } else if (!/\S+@\S+\.\S+/.test(formData.email)) {
      newErrors.email = 'Email is invalid';
    }

    if (!formData.password) {
      newErrors.password = 'Password is required';
    } else if (formData.password.length < 6) {
      newErrors.password = 'Password must be at least 6 characters';
    }

    if (formData.password !== formData.confirmPassword) {
      newErrors.confirmPassword = 'Passwords do not match';
    }

    setErrors(newErrors);
    return Object.keys(newErrors).length === 0;
  };

  const handleSubmit = (e) => {
    e.preventDefault();
    
    if (validateForm()) {
      // TODO: Implement account creation logic
      console.log('Account creation:', formData);
      alert('Account creation functionality will be implemented later');
    }
  };

  return (
    <div className="card">
      <h2>👤 Account</h2>
      <p className="subtitle">Create your account to get started</p>
      
      <form className="account-form" onSubmit={handleSubmit}>
        <div className="input-group">
          <input 
            type="text" 
            name="fullName"
            className={`text-input ${errors.fullName ? 'error' : ''}`}
            placeholder="Full Name"
            value={formData.fullName}
            onChange={handleInputChange}
          />
          <input 
            type="email" 
            name="email"
            className={`text-input ${errors.email ? 'error' : ''}`}
            placeholder="Email address"
            value={formData.email}
            onChange={handleInputChange}
          />
        </div>
        {errors.fullName && <span className="error-message">{errors.fullName}</span>}
        {errors.email && <span className="error-message">{errors.email}</span>}
        
        <div className="input-group">
          <input 
            type="password" 
            name="password"
            className={`text-input ${errors.password ? 'error' : ''}`}
            placeholder="Password"
            value={formData.password}
            onChange={handleInputChange}
          />
          <input 
            type="password" 
            name="confirmPassword"
            className={`text-input ${errors.confirmPassword ? 'error' : ''}`}
            placeholder="Confirm Password"
            value={formData.confirmPassword}
            onChange={handleInputChange}
          />
        </div>
        {errors.password && <span className="error-message">{errors.password}</span>}
        {errors.confirmPassword && <span className="error-message">{errors.confirmPassword}</span>}
        
        <button type="submit" className="btn btn-primary">
          Create Account
        </button>
        <p style={{marginTop: '20px', color: '#666'}}>
          Already have an account?{' '}
          <a href="#" onClick={(e) => {
            e.preventDefault();
            setActiveTab('login');
          }}>
            Login here
          </a>
        </p>
      </form>
    </div>
  );
}

export default Account;
