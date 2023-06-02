import logging
import os
from functools import lru_cache

import requests
from cmat.trait_mapping.zooma import get_zooma_results
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
            # Check that we've found the drug exactly (strict case-insensitive string match)
            if result['label'].lower() == drug_name.lower():
                candidates.add(result['iri'])
        # Only return a result if we can find it unambiguously
        if len(candidates) == 1:
            return candidates.pop()
    logger.warning(f'Could not find a CHEBI IRI for {drug_name}')
    return None


@lru_cache
@retry(exceptions=(ConnectionError, RequestException), tries=4, delay=2, backoff=1.2, jitter=(1, 3))
def get_efo_iri(phenotype_name):
    if not phenotype_name:
        return ''
    filters = {'ontologies': 'efo,ordo,hp,mondo',
               'required': 'cttv,eva-clinvar,clinvar-xrefs,gwas',  # TODO do we want all of these?
               'preferred': 'eva-clinvar,cttv,gwas,clinvar-xrefs'}
    zooma_results = get_zooma_results(phenotype_name, filters=filters, zooma_host='https://www.ebi.ac.uk')
    current_efo_uris = []
    for result in zooma_results:
        if result.confidence.lower() != "high":
            continue
        for mapping in result.mapping_list:
            if mapping.in_efo and mapping.is_current:
                current_efo_uris.append(mapping.uri)
    if len(current_efo_uris) > 1:
        raise ValueError(f'Found multiple high-confidence mappings for {phenotype_name}: {",".join(current_efo_uris)}')
    if len(current_efo_uris) == 0:
        return ''
    return current_efo_uris[0]
