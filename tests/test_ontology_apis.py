from unittest.mock import patch

from cmat.trait_mapping.trait import Trait, OntologyEntry

from opentargets_pharmgkb.ontology_apis import get_chebi_iri, get_efo_iri


ontology_id_regex = '(^NCIT_|^Orphanet_|^GO_|^HP_|^EFO_|^MONDO_|^DOID_|^MP_|^OTAR_|^PATO_|^OBI_|^OBA_|^OGMS_|^GSSO_|^UBERON_)'


def test_get_chebi_iri():
    # Exactly one exact match
    assert get_chebi_iri('morphine') == 'http://purl.obolibrary.org/obo/CHEBI_17303'

    # Multiple results but none that exactly match
    assert get_chebi_iri('fluorouracil') is None


def test_get_efo_iri(mappings):
    # Exactly one high-confidence match
    assert get_efo_iri('Lymphoma', mappings, ontology_id_regex) == 'http://purl.obolibrary.org/obo/MONDO_0005062'

    # No high-confidence matches
    assert get_efo_iri('neoplasms', mappings, ontology_id_regex) is None

    # Mapping found but doesn't match regex
    assert get_efo_iri('ulcer', mappings, ontology_id_regex) is None


def test_get_efo_iri_cache(mappings):
    with patch('opentargets_pharmgkb.ontology_apis.process_trait') as m_process_trait:
        # Caches result of process_trait
        finished_trait = Trait('test trait name', None, None)
        finished_trait.finished_mapping_set = {OntologyEntry('http://example.com/HP_123', 'test')}
        m_process_trait.return_value = finished_trait

        get_efo_iri('test trait name', mappings, ontology_id_regex)
        assert m_process_trait.call_count == 1
        m_process_trait.reset_mock()
        get_efo_iri('test trait name', mappings, ontology_id_regex)
        assert m_process_trait.call_count == 0

        # Also caches when no mapping found
        finished_trait.finished_mapping_set = set()
        m_process_trait.reset_mock()
        get_efo_iri('something else', mappings, ontology_id_regex)
        assert m_process_trait.call_count == 1
        m_process_trait.reset_mock()
        get_efo_iri('something else', mappings, ontology_id_regex)
        assert m_process_trait.call_count == 0
