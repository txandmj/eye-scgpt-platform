#!/usr/bin/env python3
"""
Setup script to create a minimal test model structure for development.
This creates placeholder files that match the expected structure.
"""

import json
import yaml
from pathlib import Path

def create_test_model_structure():
    """Create a minimal test model structure"""
    
    model_dir = Path("models/test_model")
    model_dir.mkdir(parents=True, exist_ok=True)
    
    # Create basic vocab.json (this would normally be much larger)
    vocab_data = {
        "<pad>": 0,
        "<cls>": 1,
        "<eoc>": 2,
        "GENE1": 3,
        "GENE2": 4,
        "GENE3": 5
    }
    
    with open(model_dir / "vocab.json", "w") as f:
        json.dump(vocab_data, f, indent=2)
    
    # Create basic args.json
    model_args = {
        "ntokens": len(vocab_data),
        "embsize": 512,
        "nhead": 8,
        "d_hid": 512,
        "nlayers": 6,
        "nlayers_cls": 2,
        "n_cls": 10,
        "dropout": 0.1,
        "pad_token": "<pad>",
        "pad_value": 0,
        "do_mvc": True,
        "do_dab": False,
        "use_batch_labels": False,
        "domain_spec_batchnorm": False,
        "input_emb_style": "continuous",
        "n_input_bins": 5,
        "cell_emb_style": "cls",
        "ecs_threshold": 0.8,
        "use_fast_transformer": False,
        "fast_transformer_backend": "flash",
        "pre_norm": False
    }
    
    with open(model_dir / "args.json", "w") as f:
        json.dump(model_args, f, indent=2)
    
    # Create basic train_args.yml
    train_args = {
        "task_info": {
            "raw_dataset_directory": "",
            "raw_dataset_name": "test_dataset",
            "save_dir": "",
            "input_dataset_directory": "",
            "id2type_json": ""
        },
        "preprocess": {
            "do_norm": True,
            "filter_gene_by_counts": 0,
            "filter_cell_by_counts": 0,
            "n_hvg": 5000,
            "hvg_flavor": "seurat_v3",
            "dataset_cell_type_col": None,
            "dataset_batch_id_col": None,
            "normalize_total": 10000,
            "_FIXED_GENE_COL": "gene_name",
            "_FIXED_CELL_TYPE_COL": "cell_type",
            "_FIXED_CELL_TYPE_ID_COL": "cell_type_id",
            "_FIXED_BATCH_ID_COL": "batch_id"
        },
        "model_parameters": {
            "load_model": str(model_dir),
            "do_train": False,
            "epochs": 1,
            "batch_size": 32,
            "lr": 0.001,
            "dropout": 0.1,
            "nlayers": 6,
            "nheads": 8,
            "embsize": 512,
            "d_hid": 512,
            "nlayers_cls": 2,
            "freeze_predecoder": False,
            "use_batch_labels": False,
            "DSBN": False,
            "pre_norm": False,
            "use_flash_attn": False,
            "fast_transformer_backend": "flash",
            "schedule_interval": 1,
            "log_interval": 10,
            "save_eval_interval": 1,
            "seed": 42
        },
        "task_configs": {
            "pad_token": "<pad>",
            "cls_token": "<cls>",
            "eoc_token": "<eoc>",
            "mask_value": -1,
            "pad_value": 0,
            "append_cls": True,
            "include_zero_gene": False,
            "max_seq_len": 2000,
            "n_input_bins": 5,
            "CLS": True,
            "input_emb_style": "continuous",
            "cell_emb_style": "cls",
            "ecs_threshold": 0.8,
            "MVC": True,
            "DAB": False
        },
        "wandb_configs": {
            "project": "test_project",
            "name": "test_run",
            "reinit": True,
            "mode": "offline"
        }
    }
    
    with open(model_dir / "dev_train_args.yml", "w") as f:
        yaml.dump(train_args, f, default_flow_style=False)
    
    # Create basic id2type.json
    id2type = {
        "0": "CellType1",
        "1": "CellType2",
        "2": "CellType3",
        "3": "CellType4",
        "4": "CellType5"
    }
    
    with open(model_dir / "id2type.json", "w") as f:
        json.dump(id2type, f, indent=2)
    
    # Create a placeholder model file (this would normally be a PyTorch model)
    with open(model_dir / "best_model.pt", "w") as f:
        f.write("# Placeholder model file\n")
        f.write("# In production, this would be a PyTorch model file\n")
        f.write("# Created by setup_model.py for development\n")
    
    print(f"✅ Created test model structure in {model_dir}")
    print("📝 Note: This is a minimal setup for development.")
    print("   For production, replace these files with your actual trained model.")

if __name__ == "__main__":
    create_test_model_structure()
