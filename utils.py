import logging
import json
import os
from config import PROGRESS_FILE

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)

def save_progress(completed_genes):
    with open(PROGRESS_FILE, 'w') as f:
        json.dump({"completed_genes": completed_genes}, f)

def load_progress():
    if os.path.exists(PROGRESS_FILE):
        with open(PROGRESS_FILE, 'r') as f:
            data = json.load(f)
            return data.get("completed_genes", [])
    return []

def flatten_json(obj, parent_key='', sep='_'):
    items = []
    for k, v in obj.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_json(v, new_key, sep=sep).items())
        elif isinstance(v, list):
            if len(v) > 0 and isinstance(v[0], dict):
                items.extend(flatten_json(v[0], new_key, sep=sep).items())
            else:
                items.append((new_key, str(v)))
        else:
            items.append((new_key, v))
    return dict(items)
