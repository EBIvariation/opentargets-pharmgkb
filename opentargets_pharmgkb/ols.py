import logging
import os
from functools import lru_cache

import requests
from requests import RequestException
from retry import retry

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

OLS_API_ROOT = 'https://www.ebi.ac.uk/ols/api'


@lru_cache
@retry(exceptions=(ConnectionError, RequestException), tries=4, delay=2, backoff=1.2, jitter=(1, 3))
def get_chebi_iri(drug_name):
    chebi_search_url = os.path.join(OLS_API_ROOT, f'search?ontology=chebi&q={drug_name}')
    response = requests.get(chebi_search_url)
    response.raise_for_status()
    data = response.json()
    if 'response' in data:
        results = data['response']['docs']
        candidates = set()
        for result in results:
            # Check that we've found the drug exactly
            # TODO this is too strict for e.g. fluorouracil, which is in drugs.tsv as 5-fluorouracil (CHEBI:46345)
            if result['label'].lower() == drug_name.lower():
                candidates.add(result['iri'])
        # Only return a result if we can find it unambiguously
        if len(candidates) == 1:
            return candidates.pop()
    logger.warning(f'Could not find a CHEBI IRI for {drug_name}')
    return None
