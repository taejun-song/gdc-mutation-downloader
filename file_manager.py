import os
import pandas as pd
import logging

class FileManager:
    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.logger = logging.getLogger(__name__)

    def save_all_mutations(self, all_mutations, filename):
        if not all_mutations:
            self.logger.warning("No mutations to save")
            return

        seen_dna_changes = {}
        unique_mutations = []

        for mut in all_mutations:
            dna_change = mut["dna_change"]
            if dna_change not in seen_dna_changes:
                seen_dna_changes[dna_change] = mut
                unique_mutations.append(mut)

        self.logger.info(f"Unique mutations: {len(unique_mutations)} (from {len(all_mutations)} total)")

        sorted_mutations = sorted(
            unique_mutations,
            key=lambda x: (
                -x["num_cohort_ssm_affected_cases"],
                -x["num_gdc_ssm_affected_cases"]
            )
        )

        df = pd.DataFrame(sorted_mutations)

        column_order = [
            "ssm_id", "gene", "dna_change", "protein_change", "type", "consequence",
            "num_cohort_ssm_affected_cases", "num_cohort_ssm_cases",
            "cohort_ssm_affected_cases_percentage",
            "num_gdc_ssm_affected_cases", "num_gdc_ssm_cases",
            "gdc_ssm_affected_cases_percentage",
            "vep_impact", "sift_impact", "sift_score",
            "polyphen_impact", "polyphen_score"
        ]

        df = df[column_order]

        output_file = os.path.join(self.output_dir, filename)
        os.makedirs(self.output_dir, exist_ok=True)
        df.to_csv(output_file, sep='\t', index=False)

        self.logger.info(f"Saved {len(sorted_mutations)} mutations to {output_file}")
        return output_file
