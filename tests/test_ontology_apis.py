from opentargets_pharmgkb.ontology_apis import get_chebi_iri, get_efo_iri


def test_get_chebi_iri():
    # Exactly one exact match
    assert get_chebi_iri('morphine') == 'http://purl.obolibrary.org/obo/CHEBI_17303'

    # Multiple results but none that exactly match
    assert get_chebi_iri('fluorouracil') is None


def test_get_efo_iri():
    # Exactly one high-confidence match
    assert get_efo_iri('Lymphoma') == 'http://www.ebi.ac.uk/efo/EFO_0000574'

    # No high-confidence matches
    assert get_efo_iri('neoplasms') is None
