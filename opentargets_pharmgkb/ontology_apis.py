import logging
from functools import lru_cache

from cmat.trait_mapping.main import process_trait
from cmat.trait_mapping.ols import get_uri_from_exact_match
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
# Only used for finished mappings, so no need to search synonyms or other ontologies
ols_query_fields = 'label'
ols_field_list = 'iri,label,ontology_name'
target_ontology = 'EFO'
preferred_ontologies = [target_ontology]


def get_chebi_iri(drug_name):
    return get_uri_from_exact_match(drug_name, ontology='chebi')


@lru_cache
@retry(exceptions=(ConnectionError, RequestException), tries=4, delay=2, backoff=1.2, jitter=(1, 3))
def get_efo_iri(phenotype_name):
    if not phenotype_name:
        return None

    # Trait to store Zooma/OxO results - other attributes not used
    trait = Trait(phenotype_name, None, None)
    processed_trait = process_trait(trait, zooma_filters, HOST, oxo_targets, oxo_distance,
                                    ols_query_fields, ols_field_list, target_ontology, preferred_ontologies)
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
