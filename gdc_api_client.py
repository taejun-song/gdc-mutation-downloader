import requests
import json
import time
import logging
from datetime import datetime, timedelta
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from config import (
    CASES_ENDPOINT, SSMS_ENDPOINT, GENES_ENDPOINT, GDC_API_BASE,
    RATE_LIMIT_REQUESTS, RATE_LIMIT_WINDOW, PAGE_SIZE, MAX_RETRIES
)

class RateLimiter:
    def __init__(self, max_requests=RATE_LIMIT_REQUESTS, time_window=RATE_LIMIT_WINDOW):
        self.max_requests = max_requests
        self.time_window = time_window
        self.requests = []

    def wait_if_needed(self):
        now = datetime.now()
        self.requests = [
            req_time for req_time in self.requests
            if now - req_time < timedelta(seconds=self.time_window)
        ]

        if len(self.requests) >= self.max_requests:
            sleep_time = self.time_window - (now - self.requests[0]).total_seconds()
            if sleep_time > 0:
                time.sleep(sleep_time)

        self.requests.append(now)

class GDCApiClient:
    def __init__(self):
        self.session = self._create_session_with_retries()
        self.rate_limiter = RateLimiter()
        self.logger = logging.getLogger(__name__)

    def _create_session_with_retries(self):
        session = requests.Session()
        retry = Retry(
            total=MAX_RETRIES,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504]
        )
        adapter = HTTPAdapter(max_retries=retry)
        session.mount('https://', adapter)
        session.mount('http://', adapter)
        return session

    def _make_request(self, endpoint, params=None):
        self.rate_limiter.wait_if_needed()
        try:
            response = self.session.get(endpoint, params=params)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            self.logger.error(f"Request failed: {e}")
            raise

    def get_primary_site_case_count(self, primary_site):
        filters = {
            "op": "=",
            "content": {
                "field": "primary_site",
                "value": primary_site
            }
        }

        params = {
            "filters": json.dumps(filters),
            "size": 0
        }

        self.logger.info(f"Getting case count for primary site: {primary_site}...")
        response = self._make_request(CASES_ENDPOINT, params)
        total = response["data"]["pagination"]["total"]
        self.logger.info(f"Total cases: {total}")
        return total

    def get_open_access_maf_cohort_cases(self, primary_site):
        """Get case IDs for cases with Open Access MAF files (matches GDC Portal Mutation Frequency app)"""
        filters = {
            "op": "and",
            "content": [
                {
                    "op": "=",
                    "content": {
                        "field": "primary_site",
                        "value": primary_site
                    }
                },
                {
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
            ]
        }

        params = {
            "filters": json.dumps(filters),
            "fields": "case_id",
            "size": 10000
        }

        self.logger.info(f"Fetching Open Access MAF cohort for {primary_site}...")
        response = self._make_request(CASES_ENDPOINT, params)

        case_ids = set()
        for hit in response.get("data", {}).get("hits", []):
            case_id = hit.get("case_id")
            if case_id:
                case_ids.add(case_id)

        self.logger.info(f"Found {len(case_ids)} cases with Open Access MAF files")
        return case_ids

    def get_all_mutated_genes_in_cohort(self, primary_site):
        """Get all genes with mutations in the cohort using aggregation"""
        ssm_occ_url = f"{GDC_API_BASE}/ssm_occurrences"

        filters = {
            "op": "=",
            "content": {
                "field": "case.primary_site",
                "value": primary_site
            }
        }

        params = {
            "filters": json.dumps(filters),
            "facets": "ssm.consequence.transcript.gene.symbol",
            "size": 0
        }

        self.logger.info(f"Fetching all genes with mutations in primary site: {primary_site}...")
        response = self._make_request(ssm_occ_url, params)

        genes = []
        aggregations = response.get("data", {}).get("aggregations", {})
        gene_agg = aggregations.get("ssm.consequence.transcript.gene.symbol", {})
        buckets = gene_agg.get("buckets", [])

        for bucket in buckets:
            gene_symbol = bucket.get("key", "")
            if gene_symbol:
                genes.append({
                    "symbol": gene_symbol,
                    "occurrence_count": bucket.get("doc_count", 0)
                })

        self.logger.info(f"Retrieved {len(genes)} genes from aggregation")
        return genes

    def get_gene_id_from_symbol(self, gene_symbol):
        """Get gene_id and cancer gene census status from gene symbol"""
        genes_url = f"{GDC_API_BASE}/genes"

        filters = {
            "op": "=",
            "content": {
                "field": "symbol",
                "value": gene_symbol
            }
        }

        params = {
            "filters": json.dumps(filters),
            "fields": "gene_id,is_cancer_gene_census",
            "size": 1
        }

        response = self._make_request(genes_url, params)
        hits = response.get("data", {}).get("hits", [])
        if hits:
            gene_data = hits[0]
            return {
                "gene_id": gene_data.get("gene_id", ""),
                "is_cancer_gene_census": gene_data.get("is_cancer_gene_census", False)
            }
        return None

    def get_gene_mutations(self, gene_id, primary_site, total_cases, min_affected_pct=0.1):
        filters = {
            "op": "and",
            "content": [
                {
                    "op": "=",
                    "content": {
                        "field": "consequence.transcript.gene.gene_id",
                        "value": gene_id
                    }
                },
                {
                    "op": "=",
                    "content": {
                        "field": "occurrence.case.primary_site",
                        "value": primary_site
                    }
                }
            ]
        }

        all_mutations = []
        from_offset = 0

        fields = "ssm_id,genomic_dna_change,mutation_subtype,consequence.transcript.gene.symbol,consequence.transcript.aa_change,consequence.transcript.consequence_type,consequence.transcript.annotation.vep_impact,consequence.transcript.annotation.sift_impact,consequence.transcript.annotation.sift_score,consequence.transcript.annotation.polyphen_impact,consequence.transcript.annotation.polyphen_score,occurrence.case.case_id"

        while True:
            params = {
                "filters": json.dumps(filters),
                "fields": fields,
                "expand": "consequence.transcript.annotation,occurrence.case",
                "format": "JSON",
                "size": PAGE_SIZE,
                "from": from_offset
            }

            response = self._make_request(SSMS_ENDPOINT, params)
            data = response.get("data", {})
            mutations = data.get("hits", [])

            if not mutations:
                break

            for mut in mutations:
                num_affected = len(mut.get("occurrence", []))
                affected_pct = (num_affected / total_cases) * 100 if total_cases > 0 else 0

                if affected_pct > min_affected_pct:
                    mut["num_affected_cases"] = num_affected
                    mut["affected_percentage"] = affected_pct
                    all_mutations.append(mut)

            pagination = data.get("pagination", {})
            total = pagination.get("total", 0)
            from_offset += PAGE_SIZE

            if from_offset >= total:
                break

        return all_mutations

    def get_gene_case_count(self, gene_id, primary_site=None, cohort_case_ids=None):
        """
        Count unique cases with mutations in a gene.
        If cohort_case_ids is provided, only count cases in that cohort.
        """
        filters_content = [{
            "op": "=",
            "content": {
                "field": "consequence.transcript.gene.gene_id",
                "value": gene_id
            }
        }]

        if primary_site:
            filters_content.append({
                "op": "=",
                "content": {
                    "field": "occurrence.case.primary_site",
                    "value": primary_site
                }
            })

        filters = {
            "op": "and",
            "content": filters_content
        }

        params = {
            "filters": json.dumps(filters),
            "fields": "occurrence.case.case_id",
            "expand": "occurrence.case",
            "size": 10000
        }

        response = self._make_request(SSMS_ENDPOINT, params)
        unique_cases = set()
        for hit in response.get("data", {}).get("hits", []):
            for occ in hit.get("occurrence", []):
                case_id = occ.get("case", {}).get("case_id")
                if case_id:
                    if cohort_case_ids is None or case_id in cohort_case_ids:
                        unique_cases.add(case_id)

        return len(unique_cases)

    def get_gdc_wide_mutation_stats(self, ssm_id):
        filters = {
            "op": "=",
            "content": {
                "field": "ssm_id",
                "value": ssm_id
            }
        }

        params = {
            "filters": json.dumps(filters),
            "fields": "ssm_id,occurrence.case.case_id",
            "size": 1
        }

        response = self._make_request(SSMS_ENDPOINT, params)
        if response["data"]["hits"]:
            mut = response["data"]["hits"][0]
            return len(mut.get("occurrence", []))
        return 0
