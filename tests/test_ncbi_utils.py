from opentargets_pharmgkb.ncbi_utils import get_spdi_coords_for_rsid


def test_get_spdi_coords_for_rsid():
    assert get_spdi_coords_for_rsid('rs35068180') == ('NC_000011.10', 102845216, 'AAAAA', ['AAAA', 'AAAAAA', 'AAAAAAA'])
    assert get_spdi_coords_for_rsid('rs3000000000') == (None, None, None, [])
