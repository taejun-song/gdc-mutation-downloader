# GDC Mutation Data Downloader

Automated tool to download mutation data from the [GDC Data Portal](https://portal.gdc.cancer.gov/) that exactly matches the Mutation Frequency app output.

## Features

- Downloads mutation data for top genes ranked by frequency in a specified cohort
- Filters by **Cancer Gene Census** genes (excludes large genes with high background mutation rates)
- Uses **Open Access MAF cohort** (matches GDC Portal Mutation Frequency app)
- Includes functional impact annotations (SIFT, PolyPhen, VEP)
- Outputs TSV file with 17 columns matching GDC Portal format
- Supports progress tracking and resume capability
- Rate limiting and automatic retry for API reliability

## Requirements

- Python 3.8+
- [uv](https://github.com/astral-sh/uv) package manager (recommended) or pip

## Installation

### Using uv (recommended)

```bash
# Install uv if not already installed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone the repository
git clone https://github.com/yourusername/download_ssm.git
cd download_ssm

# Install dependencies (uv will handle this automatically)
uv sync
```

### Using pip

```bash
# Clone the repository
git clone https://github.com/yourusername/download_ssm.git
cd download_ssm

# Create virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Configuration

Edit `config.py` to customize settings:

```python
PRIMARY_SITE = "Breast"              # Target primary site
TOP_N_GENES = 100                    # Number of top genes to process
MIN_AFFECTED_PERCENTAGE = 0.0        # Minimum affected percentage (0 = all mutations)
OUTPUT_DIR = "Breast Cancer"         # Output directory
```

## Usage

### Run the full download

```bash
uv run python main.py
```

### Run a test with 3 genes

```bash
uv run python test_download.py
```

### Resume interrupted download

The script automatically saves progress. Simply run the command again to resume from where it stopped.

## Output Format

The output TSV file contains 17 columns:

| Column | Description |
|--------|-------------|
| `ssm_id` | Mutation UUID |
| `gene` | Gene symbol |
| `dna_change` | Genomic change (HGVS notation) |
| `protein_change` | Protein change (e.g., "TP53 R175H") |
| `type` | Mutation type (e.g., "Single base substitution") |
| `consequence` | Variant consequence (e.g., "missense_variant") |
| `num_cohort_ssm_affected_cases` | Cases in cohort with this mutation |
| `num_cohort_ssm_cases` | Total cases in cohort with mutations in this gene |
| `cohort_ssm_affected_cases_percentage` | Percentage in cohort |
| `num_gdc_ssm_affected_cases` | Cases across GDC with this mutation |
| `num_gdc_ssm_cases` | Total cases across GDC with mutations in this gene |
| `gdc_ssm_affected_cases_percentage` | Percentage across GDC |
| `vep_impact` | VEP impact prediction |
| `sift_impact` | SIFT impact prediction |
| `sift_score` | SIFT score |
| `polyphen_impact` | PolyPhen impact prediction |
| `polyphen_score` | PolyPhen score |

## How It Works

### 1. Cohort Selection

The script uses the **Open Access MAF cohort**, which matches the GDC Portal's Mutation Frequency app:

```
Filters: Primary Site = "Breast" + Files.Access = "open" + Files.Data Format = "MAF"
Result: 1,307 cases for Breast cancer
```

### 2. Gene Ranking

Genes are ranked by **# SSM Affected Cases in Cohort** (number of unique cases with mutations):

```
Example for Breast cancer:
1. TP53: 447/1,307 cases (34.20%)
2. PIK3CA: 436/1,307 cases (33.36%)
3. CDH1: 165/1,307 cases (12.62%)
```

### 3. Cancer Gene Census Filter

Only genes in the [Cancer Gene Census](https://cancer.sanger.ac.uk/census) are included. This excludes:
- Large genes with high background mutation rates (TTN, MUC16, OBSCN, RYR2, etc.)
- Genes not causally implicated in cancer

This matches the GDC Portal's default filter setting.

### 4. Mutation Collection

For each gene, the script:
1. Fetches all mutations in the cohort (no frequency filter by default)
2. Gets functional impact annotations (SIFT, PolyPhen, VEP)
3. Calculates both cohort-specific and GDC-wide statistics

### 5. Output

All mutations from all genes are combined into a single TSV file, sorted by cohort frequency.

## Project Structure

```
download_ssm/
├── main.py                  # Main script
├── config.py                # Configuration settings
├── gdc_api_client.py        # GDC API client with rate limiting
├── file_manager.py          # File output handling
├── mutation_formatter.py    # Mutation data formatting
├── utils.py                 # Logging and progress tracking
├── requirements.txt         # Python dependencies
├── pyproject.toml          # uv project configuration
└── README.md               # This file
```

## API Reference

### GDC API Endpoints Used

- `/cases` - Case data and cohort filtering
- `/ssms` - Simple somatic mutation data
- `/ssm_occurrences` - Mutation occurrences (for gene aggregation)
- `/genes` - Gene information and Cancer Gene Census status

### Rate Limiting

The script implements rate limiting (5 requests per 10 seconds by default) to avoid overloading the GDC API.

## Troubleshooting

### Connection errors

If you encounter connection errors, the script will automatically retry with exponential backoff.

### Wrong case count

Ensure you're using the correct primary site name. Check available primary sites at:
https://portal.gdc.cancer.gov/

### Missing mutations

By default, `MIN_AFFECTED_PERCENTAGE = 0.0` includes all mutations. Increase this value to filter by frequency.

## Example Output

```bash
$ uv run python main.py

============================================================
GDC Mutation Data Downloader
============================================================
Target: Primary site 'Breast'
Top genes: 100
Min affected percentage: 0.0%

Fetching Open Access MAF cohort for Breast...
Found 1307 cases with Open Access MAF files

Getting gene IDs and calculating # SSM Affected Cases in Cohort...
Processing genes: 100%|██████████| 200/200 [01:45<00:00]

Top 10 genes by # SSM Affected Cases in Cohort:
  1. TP53: 447 cases
  2. PIK3CA: 436 cases
  3. CDH1: 165 cases
  4. GATA3: 158 cases
  5. KMT2C: 107 cases
  ...

Processing 46 genes
Processing genes: 100%|██████████| 46/46 [07:53<00:00]

============================================================
Total mutations collected: 1813
Output file: Breast Cancer/frequent-mutations.2025-12-23.tsv
============================================================
```

## References

- [GDC Data Portal](https://portal.gdc.cancer.gov/)
- [GDC API Documentation](https://docs.gdc.cancer.gov/API/)
- [Mutation Frequency App Documentation](https://docs.gdc.cancer.gov/Data_Portal/Users_Guide/mutation_frequency/)
- [Cancer Gene Census](https://cancer.sanger.ac.uk/census)

## License

MIT License

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
