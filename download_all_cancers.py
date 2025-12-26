import sys
import os
import json
import argparse
from datetime import datetime
from tqdm import tqdm
from gdc_api_client import GDCApiClient
from file_manager import FileManager
from mutation_formatter import format_mutation_for_output
from utils import setup_logging, save_progress, load_progress
from config import TOP_N_GENES, MIN_AFFECTED_PERCENTAGE, CASES_ENDPOINT

def get_downloaded_sites(base_dir="."):
    """Get list of primary sites that have already been downloaded"""
    downloaded = set()
    for item in os.listdir(base_dir):
        if item.endswith(" Cancer") and os.path.isdir(os.path.join(base_dir, item)):
            dir_path = os.path.join(base_dir, item)
            tsv_files = [f for f in os.listdir(dir_path) if f.endswith(".tsv")]
            if tsv_files:
                site_name = item.replace(" Cancer", "").replace("-", "/")
                downloaded.add(site_name.lower())
    return downloaded

def get_all_primary_sites(api_client):
    """Get all available primary sites from GDC API"""
    filters = {
        "op": "and",
        "content": [
            {
                "op": "=",
                "content": {
                    "field": "files.access",
                    "value": "open"
                }
            },
            {
                "op": "=",
                "content": {
                    "field": "files.data_format",
                    "value": "MAF"
                }
            }
        ]
    }

    params = {
        "filters": json.dumps(filters),
        "facets": "primary_site",
        "size": 0
    }

    response = api_client._make_request(CASES_ENDPOINT, params)

    primary_sites = []
    aggregations = response.get("data", {}).get("aggregations", {})
    site_agg = aggregations.get("primary_site", {})
    buckets = site_agg.get("buckets", [])

    for bucket in buckets:
        site = bucket.get("key", "")
        count = bucket.get("doc_count", 0)
        if site and count > 0:
            primary_sites.append({
                "site": site,
                "case_count": count
            })

    return primary_sites

def download_cancer_type(api_client, primary_site, top_n_genes, min_affected_pct, logger):
    """Download mutation data for a single cancer type"""
    output_dir = f"{primary_site.replace('/', '-')} Cancer"
    file_manager = FileManager(output_dir)

    try:
        logger.info("=" * 60)
        logger.info(f"Processing: {primary_site}")
        logger.info("=" * 60)

        cohort_case_ids = api_client.get_open_access_maf_cohort_cases(primary_site)
        total_cohort_cases = len(cohort_case_ids)

        if total_cohort_cases == 0:
            logger.warning(f"No Open Access MAF cases found for {primary_site}. Skipping.")
            return

        logger.info(f"Fetching all Cancer Gene Census genes...")
        all_genes = api_client.get_all_cancer_gene_census_genes()
        logger.info(f"Retrieved {len(all_genes)} Cancer Gene Census genes")
        logger.info("")

        logger.info("Calculating # SSM Affected Cases in Cohort for each gene...")
        genes_with_counts = []
        for gene in tqdm(all_genes, desc=f"Processing genes for {primary_site}"):
            gene_symbol = gene["symbol"].upper()
            gene_id = gene["gene_id"]

            cohort_gene_cases = api_client.get_gene_case_count(gene_id, primary_site, cohort_case_ids)

            if cohort_gene_cases == 0:
                continue

            genes_with_counts.append({
                "symbol": gene_symbol,
                "gene_id": gene_id,
                "cohort_affected_cases": cohort_gene_cases
            })

        genes_with_counts.sort(key=lambda x: x["cohort_affected_cases"], reverse=True)
        top_genes = genes_with_counts[:top_n_genes]

        logger.info("")
        logger.info(f"Top 10 genes by # SSM Affected Cases in Cohort:")
        for i, gene in enumerate(top_genes[:10], 1):
            logger.info(f"  {i}. {gene['symbol']}: {gene['cohort_affected_cases']} cases")
        logger.info("")

        genes_to_process = top_genes
        logger.info(f"Processing {len(genes_to_process)} genes")
        logger.info("")

        all_mutations = []

        for idx, gene in enumerate(tqdm(genes_to_process, desc=f"Processing mutations for {primary_site}"), 1):
            gene_symbol = gene["symbol"]
            gene_id = gene["gene_id"]
            cohort_gene_cases = gene["cohort_affected_cases"]

            logger.info(f"[{idx}/{len(genes_to_process)}] Processing {gene_symbol}...")

            gdc_gene_cases = api_client.get_gene_case_count(gene_id, None)

            logger.info(f"Gene case counts - Cohort: {cohort_gene_cases}, GDC: {gdc_gene_cases}")

            mutations = api_client.get_gene_mutations(
                gene_id=gene_id,
                primary_site=primary_site,
                total_cases=total_cohort_cases,
                min_affected_pct=min_affected_pct
            )

            if not mutations:
                logger.warning(f"No mutations found for {gene_symbol}")
                continue

            logger.info(f"Found {len(mutations)} mutations for {gene_symbol}")

            for mut in mutations:
                gdc_affected = api_client.get_gdc_wide_mutation_stats(mut["ssm_id"])

                formatted_mut = format_mutation_for_output(
                    mut, gene_symbol, cohort_gene_cases, gdc_gene_cases, gdc_affected
                )
                all_mutations.append(formatted_mut)

        logger.info("")
        logger.info("=" * 60)
        logger.info(f"Total mutations collected: {len(all_mutations)}")

        date_str = datetime.now().strftime("%Y-%m-%d")
        filename = f"frequent-mutations.{date_str}.tsv"

        file_manager.save_all_mutations(all_mutations, filename)

        logger.info("=" * 60)
        logger.info(f"Download complete for {primary_site}!")
        logger.info(f"Output file: {output_dir}/{filename}")
        logger.info("=" * 60)
        logger.info("")

    except Exception as e:
        logger.error(f"Error processing {primary_site}: {e}", exc_info=True)

def main():
    parser = argparse.ArgumentParser(
        description="Download mutation data for all cancer types from GDC Data Portal"
    )
    parser.add_argument(
        "--top-genes",
        type=int,
        default=TOP_N_GENES,
        help=f"Number of top genes to download per cancer type (default: {TOP_N_GENES})"
    )
    parser.add_argument(
        "--min-affected-pct",
        type=float,
        default=MIN_AFFECTED_PERCENTAGE,
        help=f"Minimum affected percentage (default: {MIN_AFFECTED_PERCENTAGE})"
    )
    parser.add_argument(
        "--min-cases",
        type=int,
        default=0,
        help="Minimum number of cases required for a cancer type (default: 0)"
    )
    parser.add_argument(
        "--skip-downloaded",
        action="store_true",
        default=True,
        help="Skip already downloaded primary sites (default: True)"
    )
    parser.add_argument(
        "--no-skip-downloaded",
        action="store_true",
        help="Re-download all primary sites including already downloaded ones"
    )
    args = parser.parse_args()

    top_n_genes = args.top_genes
    min_affected_pct = args.min_affected_pct
    min_cases = args.min_cases
    skip_downloaded = args.skip_downloaded and not args.no_skip_downloaded

    logger = setup_logging()
    logger.info("=" * 60)
    logger.info("GDC Mutation Data Downloader - ALL CANCER TYPES")
    logger.info("=" * 60)

    api_client = GDCApiClient()

    try:
        downloaded_sites = set()
        if skip_downloaded:
            downloaded_sites = get_downloaded_sites()
            if downloaded_sites:
                logger.info(f"Found {len(downloaded_sites)} already downloaded sites (will skip)")

        logger.info("Fetching all primary sites from GDC API...")
        primary_sites = get_all_primary_sites(api_client)
        logger.info(f"Found {len(primary_sites)} primary sites")
        logger.info("")

        filtered_sites = [site for site in primary_sites if site["case_count"] >= min_cases]
        if skip_downloaded:
            before_skip = len(filtered_sites)
            filtered_sites = [site for site in filtered_sites if site["site"].lower() not in downloaded_sites]
            skipped_count = before_skip - len(filtered_sites)
            if skipped_count > 0:
                logger.info(f"Skipping {skipped_count} already downloaded sites")

        logger.info(f"Processing {len(filtered_sites)} sites with >= {min_cases} cases")
        logger.info("")

        logger.info("Primary sites to process:")
        for i, site in enumerate(filtered_sites, 1):
            logger.info(f"  {i}. {site['site']}: {site['case_count']} cases")
        logger.info("")

        for idx, site_info in enumerate(filtered_sites, 1):
            primary_site = site_info["site"]
            logger.info(f"\n[{idx}/{len(filtered_sites)}] Starting download for: {primary_site}")

            download_cancer_type(
                api_client,
                primary_site,
                top_n_genes,
                min_affected_pct,
                logger
            )

        logger.info("=" * 60)
        logger.info("ALL DOWNLOADS COMPLETE!")
        logger.info(f"Processed {len(filtered_sites)} cancer types")
        logger.info("=" * 60)

    except KeyboardInterrupt:
        logger.warning("\nDownload interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error occurred: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
