def format_mutation_for_output(mutation, gene_symbol, cohort_gene_cases, gdc_gene_cases, gdc_affected):
    consequence = mutation.get("consequence", [{}])[0] if mutation.get("consequence") else {}
    transcript = consequence.get("transcript", {})
    annotation = transcript.get("annotation", {})

    aa_change = transcript.get("aa_change", "")
    protein_change = f"{gene_symbol} {aa_change}" if aa_change else ""

    vep_impact = annotation.get("vep_impact", "")
    sift_impact = annotation.get("sift_impact", "")
    sift_score = annotation.get("sift_score", "")
    polyphen_impact = annotation.get("polyphen_impact", "")
    polyphen_score = annotation.get("polyphen_score", "")

    cohort_affected = mutation.get("num_affected_cases", 0)
    cohort_pct = (cohort_affected / cohort_gene_cases * 100) if cohort_gene_cases > 0 else 0
    gdc_pct = (gdc_affected / gdc_gene_cases * 100) if gdc_gene_cases > 0 else 0

    return {
        "ssm_id": mutation.get("ssm_id", ""),
        "gene": gene_symbol,
        "dna_change": mutation.get("genomic_dna_change", ""),
        "protein_change": protein_change,
        "type": mutation.get("mutation_subtype", ""),
        "consequence": transcript.get("consequence_type", ""),
        "num_cohort_ssm_affected_cases": cohort_affected,
        "num_cohort_ssm_cases": cohort_gene_cases,
        "cohort_ssm_affected_cases_percentage": round(cohort_pct, 2),
        "num_gdc_ssm_affected_cases": gdc_affected,
        "num_gdc_ssm_cases": gdc_gene_cases,
        "gdc_ssm_affected_cases_percentage": round(gdc_pct, 2),
        "vep_impact": vep_impact,
        "sift_impact": sift_impact,
        "sift_score": sift_score,
        "polyphen_impact": polyphen_impact,
        "polyphen_score": polyphen_score
    }
