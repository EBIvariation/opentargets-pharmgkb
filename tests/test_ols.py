from opentargets_pharmgkb.ols import get_chebi_iri


def test_get_chebi_iri():
    # Exactly one exact match
    assert get_chebi_iri('morphine') == 'http://purl.obolibrary.org/obo/CHEBI_17303'

    # Multiple results but none that exactly match
    assert get_chebi_iri('fluorouracil') is None
