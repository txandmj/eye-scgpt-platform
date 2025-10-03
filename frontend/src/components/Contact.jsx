import React, { useState } from 'react';

function Contact() {
  const [formData, setFormData] = useState({
    name: '',
    email: '',
    message: ''
  });

  const [isSubmitting, setIsSubmitting] = useState(false);
  const [submitStatus, setSubmitStatus] = useState(null);

  const handleInputChange = (e) => {
    const { name, value } = e.target;
    setFormData(prev => ({
      ...prev,
      [name]: value
    }));
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    setIsSubmitting(true);
    
    try {
      // TODO: Implement actual form submission
      console.log('Contact form submission:', formData);
      
      // Simulate API call
      await new Promise(resolve => setTimeout(resolve, 1000));
      
      setSubmitStatus('success');
      setFormData({ name: '', email: '', message: '' });
    } catch (error) {
      setSubmitStatus('error');
    } finally {
      setIsSubmitting(false);
    }
  };

  return (
    <div className="card">
      <h2>📞 Contact Us</h2>
      <p className="subtitle">Get in touch with our team</p>
      
      <div className="contact-info">
        <div className="contact-section">
          <h3>📧 Email Support</h3>
          <p>support@eye-scgpt.com</p>
        </div>
        
        <div className="contact-section">
          <h3>💬 General Inquiries</h3>
          <p>info@eye-scgpt.com</p>
        </div>
        
        <div className="contact-section">
          <h3>🔬 Technical Support</h3>
          <p>tech@eye-scgpt.com</p>
        </div>
        
        <div className="contact-section">
          <h3>📋 Feedback Form</h3>
          <form onSubmit={handleSubmit}>
            <div className="input-group">
              <input 
                type="text" 
                name="name"
                className="text-input" 
                placeholder="Your Name"
                value={formData.name}
                onChange={handleInputChange}
                required
              />
              <input 
                type="email" 
                name="email"
                className="text-input" 
                placeholder="Your Email"
                value={formData.email}
                onChange={handleInputChange}
                required
              />
            </div>
            <textarea 
              name="message"
              className="text-input" 
              rows="4" 
              placeholder="Your Message"
              style={{width: '100%', marginTop: '10px'}}
              value={formData.message}
              onChange={handleInputChange}
              required
            />
            <button 
              type="submit" 
              className="btn btn-primary" 
              style={{marginTop: '15px'}}
              disabled={isSubmitting}
            >
              {isSubmitting ? 'Sending...' : 'Send Message'}
            </button>
          </form>
          
          {submitStatus === 'success' && (
            <div className="message success" style={{marginTop: '15px'}}>
              Thank you for your message! We'll get back to you soon.
            </div>
          )}
          
          {submitStatus === 'error' && (
            <div className="message error" style={{marginTop: '15px'}}>
              Sorry, there was an error sending your message. Please try again.
            </div>
          )}
        </div>
      </div>
    </div>
  );
}

export default Contact;
