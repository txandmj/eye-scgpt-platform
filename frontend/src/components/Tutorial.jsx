import React from 'react';

function Tutorial() {
  return (
    <div className="tutorial-container">
      <div className="card">
        <h2>📚 Tutorial & Documentation</h2>
        <p className="subtitle">
          Learn how to use the Eye-scGPT Annotation Platform
        </p>

        <section className="tutorial-section">
          <h3>🎯 Overview</h3>
          <p>
            The Eye-scGPT Annotation Platform provides automated cell type annotation
            for single-cell omics data using our fine-tuned eye-scGPT model. This platform
            streamlines the analysis workflow for researchers studying eye tissue at single-cell resolution.
          </p>
        </section>

        <section className="tutorial-section">
          <h3>🚀 Quick Start Guide</h3>
          
          <div className="step">
            <h4>Step 1: Prepare Your Data</h4>
            <p>Ensure your single-cell data meets these requirements:</p>
            <ul>
              <li><strong>Format:</strong> .h5ad (Scanpy), .h5, or .rds (Seurat)</li>
              <li><strong>Preprocessing:</strong> Data should be normalized and log-transformed</li>
              <li><strong>Size:</strong> Maximum 5GB per file</li>
              <li><strong>Quality Control:</strong> Remove low-quality cells and doublets</li>
            </ul>
          </div>

          <div className="step">
            <h4>Step 2: Upload Your File</h4>
            <ol>
              <li>Navigate to the <strong>Upload</strong> tab</li>
              <li>Click "Choose File" and select your data file</li>
              <li>Click "Upload & Annotate" to begin processing</li>
              <li>Save the generated Job ID for later retrieval</li>
            </ol>
          </div>

          <div className="step">
            <h4>Step 3: Wait for Processing</h4>
            <p>
              The annotation process typically takes 2-10 minutes depending on dataset size.
              You can close the browser and return later using your Job ID.
            </p>
          </div>

          <div className="step">
            <h4>Step 4: Download Results</h4>
            <ol>
              <li>Go to the <strong>Download</strong> tab</li>
              <li>Enter your Job ID</li>
              <li>Click "Fetch Results" to view available files</li>
              <li>Download your annotated data and visualizations</li>
            </ol>
          </div>
        </section>

        <section className="tutorial-section">
          <h3>📊 Data Preparation Guidelines</h3>
          
          <h4>For Scanpy Users (Python):</h4>
          <div className="code-block">
            <pre>{`import scanpy as sc

# Load your data
adata = sc.read_10x_mtx('path/to/data/')

# Quality control
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Save for upload
adata.write('preprocessed_data.h5ad')`}</pre>
          </div>

          <h4>For Seurat Users (R):</h4>
          <div className="code-block">
            <pre>{`library(Seurat)

# Load and preprocess data
seurat_obj <- CreateSeuratObject(counts = raw_data)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# Save for upload
saveRDS(seurat_obj, "preprocessed_data.rds")`}</pre>
          </div>
        </section>

        <section className="tutorial-section">
          <h3>📈 Understanding Your Results</h3>
          
          <div className="result-explanation">
            <h4>1. Annotations CSV/Excel</h4>
            <p>
              Contains predicted cell types for each cell in your dataset.
              Includes confidence scores and metadata.
            </p>
          </div>

          <div className="result-explanation">
            <h4>2. Annotated H5AD</h4>
            <p>
              Complete AnnData object with annotations embedded. Can be loaded directly
              into Scanpy for further analysis.
            </p>
          </div>

          <div className="result-explanation">
            <h4>3. UMAP Visualization</h4>
            <p>
              Color-coded UMAP plot showing the spatial distribution of predicted
              cell types in your dataset.
            </p>
          </div>
        </section>

        <section className="tutorial-section">
          <h3>💡 Best Practices</h3>
          <ul>
            <li>Always perform quality control before uploading</li>
            <li>Normalize and log-transform your data</li>
            <li>Keep your Job ID in a safe place</li>
            <li>Download all result files for comprehensive analysis</li>
            <li>Verify annotations with known marker genes</li>
            <li>Consider batch effects if combining multiple samples</li>
          </ul>
        </section>

        <section className="tutorial-section">
          <h3>❓ FAQ</h3>
          
          <div className="faq-item">
            <h4>How long are results stored?</h4>
            <p>Results are retained for 30 days after processing.</p>
          </div>

          <div className="faq-item">
            <h4>Can I upload multiple files?</h4>
            <p>Currently, each upload processes one file at a time. For batch processing, upload files sequentially.</p>
          </div>

          <div className="faq-item">
            <h4>What if my file is larger than 5GB?</h4>
            <p>Consider downsampling your dataset or splitting it into smaller batches.</p>
          </div>

          <div className="faq-item">
            <h4>Which cell types can the model identify?</h4>
            <p>The eye-scGPT model is trained on eye tissue data and can identify major retinal cell types including photoreceptors, bipolar cells, ganglion cells, amacrine cells, horizontal cells, Müller glia, and more.</p>
          </div>

          <div className="faq-item">
            <h4>Can I use this for non-eye tissues?</h4>
            <p>This model is specifically fine-tuned for eye/retinal tissue. Performance on other tissues may vary.</p>
          </div>
        </section>

        <section className="tutorial-section">
          <h3>🔗 Additional Resources</h3>
          <ul>
            <li><a href="https://github.com/RCHENLAB/scGPT_fineTune_protocol" target="_blank" rel="noopener noreferrer">scGPT Fine-Tuning Protocol</a></li>
            <li><a href="https://scanpy.readthedocs.io/" target="_blank" rel="noopener noreferrer">Scanpy Documentation</a></li>
            <li><a href="https://satijalab.org/seurat/" target="_blank" rel="noopener noreferrer">Seurat Documentation</a></li>
          </ul>
        </section>

        <section className="tutorial-section">
          <h3>📧 Support</h3>
          <p>
            For technical support, questions, or feedback, please contact our research team
            or open an issue on our GitHub repository.
          </p>
        </section>
      </div>
    </div>
  );
}

export default Tutorial;