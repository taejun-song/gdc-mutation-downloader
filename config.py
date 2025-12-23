PRIMARY_SITE = "Breast"
TOP_N_GENES = 100
MIN_AFFECTED_PERCENTAGE = 0.0
OUTPUT_DIR = "Breast Cancer"
OUTPUT_FILE = "frequent-mutations"
DATE_FORMAT = "%Y-%m-%d"

GDC_API_BASE = "https://api.gdc.cancer.gov"
CASES_ENDPOINT = f"{GDC_API_BASE}/cases"
SSMS_ENDPOINT = f"{GDC_API_BASE}/ssms"
GENES_ENDPOINT = f"{GDC_API_BASE}/analysis/top_mutated_genes_by_project"

RATE_LIMIT_REQUESTS = 5
RATE_LIMIT_WINDOW = 1.0
PAGE_SIZE = 1000
MAX_RETRIES = 5

PROGRESS_FILE = "progress.json"
