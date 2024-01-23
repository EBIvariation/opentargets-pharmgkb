import logging
import os
from functools import lru_cache

import requests
from cmat.trait_mapping.main import process_trait
from cmat.trait_mapping.trait import Trait
from requests import RequestException
from retry import retry

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

HOST = 'https://www.ebi.ac.uk'
OLS_API_ROOT = f'{HOST}/ols4/api'

# Defaults from CMAT for getting EFO mappings that don't require manual curation
zooma_filters = {'ontologies': 'efo,ordo,hp,mondo',
                 'required': 'cttv,eva-clinvar,clinvar-xrefs,gwas',
                 'preferred': 'eva-clinvar,cttv,gwas,clinvar-xrefs'}
oxo_targets = ['Orphanet', 'efo', 'hp', 'mondo']
oxo_distance = 1


@lru_cache
@retry(exceptions=(ConnectionError, RequestException), tries=4, delay=2, backoff=1.2, jitter=(1, 3))
def get_chebi_iri(drug_name):
    chebi_search_url = os.path.join(OLS_API_ROOT, f'search?ontology=chebi&q={drug_name}&queryFields=label&exact=true')
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
        return None

    # Trait to store Zooma/OxO results - other attributes not used
    trait = Trait(phenotype_name, None, None)
    processed_trait = process_trait(trait, zooma_filters, HOST, oxo_targets, oxo_distance)
    if processed_trait.is_finished:
        efo_uris = [ontology_entry.uri for ontology_entry in processed_trait.finished_mapping_set]
        if len(efo_uris) > 1:
            # Don't expect multiple mappings for PharmGKB phenotypes
            logger.warning(f'Found multiple mappings for {phenotype_name}: {",".join(efo_uris)}')
            return None
        return efo_uris[0]
    else:
        logger.warning(f'Found no mappings for {phenotype_name}')
        return None
