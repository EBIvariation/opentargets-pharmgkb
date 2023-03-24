from opentargets_pharmgkb.variant_coordinates import get_coordinates_for_clinical_annotation


def test_get_coordinates_single_alt():
    # rs4343: G/A
    # When there's only one alt we don't need the genotypes
    assert get_coordinates_for_clinical_annotation('rs4343', []) == '17_63488670_G_A'


def test_get_coordinates_multiple_alts():
    # rs4659982: T/A/C
    assert get_coordinates_for_clinical_annotation('rs4659982', ['TT', 'TA', 'AA']) == '1_240566955_T_A'
    # Use of different alt allele generates a different identifier
    assert get_coordinates_for_clinical_annotation('rs4659982', ['TT', 'TC', 'CC']) == '1_240566955_T_C'
    # Only two alleles in genotypes, but neither is reference which is always included
    assert get_coordinates_for_clinical_annotation('rs4659982', ['AA', 'CC']) == None


def test_get_coordinates_deletion():
    # rs70991108: -/TCGCGCGTCCCGCCCAGGT/TGGCGCCTCCCGCCCAGGT/TGGCGCGTCCCGCCCAGGT
    assert get_coordinates_for_clinical_annotation(
        'rs70991108',
        ['TGGCGCGTCCCGCCCAGGT/TGGCGCGTCCCGCCCAGGT', 'TGGCGCGTCCCGCCCAGGT/del', 'del/del']
    ) == '5_80654345_-_TGGCGCGTCCCGCCCAGGT'
