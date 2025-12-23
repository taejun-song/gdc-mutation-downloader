import sys
import argparse
from datetime import datetime
from gdc_api_client import GDCApiClient
from file_manager import FileManager
from mutation_formatter import format_mutation_for_output
from utils import setup_logging
from config import PRIMARY_SITE, MIN_AFFECTED_PERCENTAGE

def main():
    parser = argparse.ArgumentParser(
        description="Test GDC mutation data downloader with 3 genes"
    )
    parser.add_argument(
        "--primary-site",
        default=PRIMARY_SITE,
        help=f"Primary site to query (default: {PRIMARY_SITE})"
    )
    parser.add_argument(
        "--output-dir",
        default="Test_Output",
        help="Output directory (default: Test_Output)"
    )

    args = parser.parse_args()

    primary_site = args.primary_site
    output_dir = args.output_dir

    logger = setup_logging()
    logger.info("=" * 60)
    logger.info("GDC Mutation Data Downloader - TEST MODE")
    logger.info("Testing with 3 genes")
    logger.info("=" * 60)

    api_client = GDCApiClient()
    file_manager = FileManager(output_dir)

    try:
        logger.info(f"Target: Primary site '{primary_site}'")
        logger.info(f"Output directory: {output_dir}")
        logger.info("")

        cohort_case_ids = api_client.get_open_access_maf_cohort_cases(primary_site)
        total_cohort_cases = len(cohort_case_ids)
        logger.info("")

        logger.info(f"Fetching all genes with mutations in cohort...")
        all_genes = api_client.get_all_mutated_genes_in_cohort(primary_site)
        logger.info(f"Retrieved {len(all_genes)} genes from aggregation")
        logger.info("")

        logger.info("Getting gene IDs for top 3 genes...")
        genes_with_counts = []
        for gene in all_genes[:10]:
            gene_symbol = gene["symbol"].upper()

            gene_data = api_client.get_gene_id_from_symbol(gene_symbol)
            if not gene_data:
                logger.warning(f"Could not find gene_id for {gene_symbol}, skipping")
                continue

            if not gene_data.get("is_cancer_gene_census", False):
                continue

            gene_id = gene_data["gene_id"]

            cohort_gene_cases = api_client.get_gene_case_count(gene_id, primary_site, cohort_case_ids)
            genes_with_counts.append({
                "symbol": gene_symbol,
                "gene_id": gene_id,
                "cohort_affected_cases": cohort_gene_cases
            })

            if len(genes_with_counts) >= 3:
                break

        genes_with_counts.sort(key=lambda x: x["cohort_affected_cases"], reverse=True)

        logger.info("")
        logger.info(f"Top 3 genes by # SSM Affected Cases in Cohort:")
        for i, gene in enumerate(genes_with_counts, 1):
            logger.info(f"  {i}. {gene['symbol']}: {gene['cohort_affected_cases']} cases")
        logger.info("")

        all_mutations = []

        for idx, gene in enumerate(genes_with_counts, 1):
            gene_symbol = gene["symbol"]
            gene_id = gene["gene_id"]
            cohort_gene_cases = gene["cohort_affected_cases"]

            logger.info(f"[{idx}/{len(genes_with_counts)}] Processing {gene_symbol}...")

            gdc_gene_cases = api_client.get_gene_case_count(gene_id, None)

            logger.info(f"Gene case counts - Cohort: {cohort_gene_cases}, GDC: {gdc_gene_cases}")

            mutations = api_client.get_gene_mutations(
                gene_id=gene_id,
                primary_site=primary_site,
                total_cases=total_cohort_cases,
                min_affected_pct=MIN_AFFECTED_PERCENTAGE
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
        filename = f"test-mutations.{date_str}.tsv"

        file_manager.save_all_mutations(all_mutations, filename)

        logger.info("=" * 60)
        logger.info("Test complete!")
        logger.info(f"Output file: {output_dir}/{filename}")
        logger.info("=" * 60)

    except Exception as e:
        logger.error(f"Error occurred: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
