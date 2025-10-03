import sys
import os

# Add the scGPT protocol to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../scGPT_fineTune_protocol'))

# Import from the GitHub repo
from scGPT_fineTune_protocol.protocol_inference import run_inference
# OR import specific modules
from scgpt.model import TransformerModel
from scgpt.tokenizer import GeneVocab

class scGPTService:
    def __init__(self, model_path: str = "./models/eye-scgpt-model"):
        self.model_path = Path(model_path)
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        
        # Load model using the GitHub repo's method
        self.model = self._load_model()
        
    def _load_model(self):
        # Use the loading method from the GitHub repository
        # This depends on how the repository structures their code
        model = TransformerModel.load_from_checkpoint(
            self.model_path / "model.pt"
        )
        model.to(self.device)
        model.eval()
        return model
    
    def predict_cell_types(self, adata: sc.AnnData):
        # Use the inference method from GitHub repo
        # Replace the mock implementation with actual inference
        
        # Example (adjust based on actual GitHub code):
        predictions = self.model.predict(adata)
        confidence_scores = self.model.get_confidence_scores(adata)
        
        return {
            "predictions": predictions,
            "confidence": confidence_scores
        }