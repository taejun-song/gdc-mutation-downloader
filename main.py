import sys
from datetime import datetime
from tqdm import tqdm
from gdc_api_client import GDCApiClient
from file_manager import FileManager
from mutation_formatter import format_mutation_for_output
from utils import setup_logging, save_progress, load_progress
from config import PRIMARY_SITE, OUTPUT_DIR, TOP_N_GENES, MIN_AFFECTED_PERCENTAGE

def main():
    logger = setup_logging()
    logger.info("=" * 60)
    logger.info("GDC Mutation Data Downloader")
    logger.info("=" * 60)

    api_client = GDCApiClient()
    file_manager = FileManager(OUTPUT_DIR)

    try:
        completed_genes = load_progress()
        logger.info(f"Loaded progress: {len(completed_genes)} genes already completed")

        logger.info(f"Target: Primary site '{PRIMARY_SITE}'")
        logger.info(f"Top genes: {TOP_N_GENES}")
        logger.info(f"Min affected percentage: {MIN_AFFECTED_PERCENTAGE}%")
        logger.info("")

        # Get Open Access MAF cohort (matches GDC Portal Mutation Frequency app)
        cohort_case_ids = api_client.get_open_access_maf_cohort_cases(PRIMARY_SITE)
        total_cohort_cases = len(cohort_case_ids)
        logger.info("")

        logger.info(f"Fetching all genes with mutations in cohort...")
        all_genes = api_client.get_all_mutated_genes_in_cohort(PRIMARY_SITE)
        logger.info(f"Retrieved {len(all_genes)} genes from aggregation")
        logger.info("")

        logger.info("Getting gene IDs and calculating # SSM Affected Cases in Cohort...")
        genes_with_counts = []
        for gene in tqdm(all_genes, desc="Processing genes"):
            gene_symbol = gene["symbol"].upper()

            if gene_symbol in completed_genes:
                continue

            # Get gene_id and cancer gene census status from symbol
            gene_data = api_client.get_gene_id_from_symbol(gene_symbol)
            if not gene_data:
                logger.warning(f"Could not find gene_id for {gene_symbol}, skipping")
                continue

            # Filter by Cancer Gene Census (matches GDC Portal Mutation Frequency default filter)
            if not gene_data.get("is_cancer_gene_census", False):
                continue

            gene_id = gene_data["gene_id"]

            # Get unique case count (filtered by Open Access MAF cohort)
            cohort_gene_cases = api_client.get_gene_case_count(gene_id, PRIMARY_SITE, cohort_case_ids)
            genes_with_counts.append({
                "symbol": gene_symbol,
                "gene_id": gene_id,
                "cohort_affected_cases": cohort_gene_cases
            })

        genes_with_counts.sort(key=lambda x: x["cohort_affected_cases"], reverse=True)
        top_genes = genes_with_counts[:TOP_N_GENES]

        logger.info("")
        logger.info(f"Top 10 genes by # SSM Affected Cases in Cohort:")
        for i, gene in enumerate(top_genes[:10], 1):
            logger.info(f"  {i}. {gene['symbol']}: {gene['cohort_affected_cases']} cases")
        logger.info("")

        genes_to_process = top_genes
        logger.info(f"Processing {len(genes_to_process)} genes")
        logger.info("")

        all_mutations = []

        for idx, gene in enumerate(tqdm(genes_to_process, desc="Processing genes"), 1):
            gene_symbol = gene["symbol"]
            gene_id = gene["gene_id"]
            cohort_gene_cases = gene["cohort_affected_cases"]

            logger.info(f"[{idx}/{len(genes_to_process)}] Processing {gene_symbol}...")

            gdc_gene_cases = api_client.get_gene_case_count(gene_id, None)

            logger.info(f"Gene case counts - Cohort: {cohort_gene_cases}, GDC: {gdc_gene_cases}")

            mutations = api_client.get_gene_mutations(
                gene_id=gene_id,
                primary_site=PRIMARY_SITE,
                total_cases=total_cohort_cases,
                min_affected_pct=MIN_AFFECTED_PERCENTAGE
            )

            if not mutations:
                logger.warning(f"No mutations found for {gene_symbol}")
                completed_genes.append(gene_symbol)
                save_progress(completed_genes)
                continue

            logger.info(f"Found {len(mutations)} mutations for {gene_symbol}")

            for mut in mutations:
                gdc_affected = api_client.get_gdc_wide_mutation_stats(mut["ssm_id"])

                formatted_mut = format_mutation_for_output(
                    mut, gene_symbol, cohort_gene_cases, gdc_gene_cases, gdc_affected
                )
                all_mutations.append(formatted_mut)

            completed_genes.append(gene_symbol)
            save_progress(completed_genes)

        logger.info("")
        logger.info("=" * 60)
        logger.info(f"Total mutations collected: {len(all_mutations)}")

        date_str = datetime.now().strftime("%Y-%m-%d")
        filename = f"frequent-mutations.{date_str}.tsv"

        file_manager.save_all_mutations(all_mutations, filename)

        logger.info("=" * 60)
        logger.info("Download complete!")
        logger.info(f"Output file: {OUTPUT_DIR}/{filename}")
        logger.info("=" * 60)

    except KeyboardInterrupt:
        logger.warning("\nDownload interrupted by user")
        logger.info(f"Progress saved. Run again to resume")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error occurred: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
