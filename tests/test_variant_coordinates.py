from opentargets_pharmgkb.variant_coordinates import Fasta


def test_get_coordinates(fasta: Fasta):
    assert fasta.get_coordinates_for_clinical_annotation(
        'rs1051266',
        'NC_000021.9:33341701',
        ['GG', 'GT', 'TT']) == '21_33341701_G_T'
    # Use of different alt allele in genotypes generates a different identifier
    assert fasta.get_coordinates_for_clinical_annotation(
        'rs1051266',
        'NC_000021.9:33341701',
        ['GG', 'GA', 'AA']) == '21_33341701_G_A'


def test_get_coordinates_multiple_alts(fasta: Fasta):
    # More than two alleles present in genotypes - TODO update once we decide how to treat these
    assert fasta.get_coordinates_for_clinical_annotation(
        'rs1051266',
        'NC_000021.9:33341701',
        ['GG', 'GT', 'GA', 'TA']) == '21_33341701_G_A'


def test_get_coordinates_deletion(fasta: Fasta):
    assert fasta.get_coordinates_for_clinical_annotation(
        'rs1051266',
        'NC_000021.9:33341701',
        ['TGGCGCGTCCCGCCCAGGT/TGGCGCGTCCCGCCCAGGT', 'TGGCGCGTCCCGCCCAGGT/del', 'del/del']
    ) == '21_33341700_A_ATGGCGCGTCCCGCCCAGGT'


def test_get_coordinates_hybrid_genotype_format(fasta: Fasta):
    assert fasta.get_coordinates_for_clinical_annotation(
        'rs35068180',
        'NC_000021.9:45514917',
        ['A/del', 'AA', 'del/del']
    ) == '21_45514916_T_TA'


def test_get_coordinates_range_location(fasta: Fasta):
    assert fasta.get_coordinates_for_clinical_annotation(
        'rs1051266',
        'NC_000021.9:33341701_33341703',
        ['GAC/GAC', 'GAC/del', 'del/del']
    ) == '21_33341700_AGAC_A'


def test_get_coordinates_not_match_reference(fasta: Fasta):
    assert fasta.get_coordinates_for_clinical_annotation(
        'rs1051266',
        'NC_000021.9:33341701_33341703',
        ['TTT/TTT', 'TTT/del', 'del/del']
    ) == None
