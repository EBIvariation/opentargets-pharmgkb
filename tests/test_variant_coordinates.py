from opentargets_pharmgkb.variant_coordinates import Fasta


def test_get_coordinates_single_alt(fasta: Fasta):
    # rs4343: G/A
    assert fasta.get_coordinates_for_clinical_annotation('rs4343', 'NC_000017.11:63488670', ['GG', 'GA', 'AA']) == '17_63488670_G_A'


def test_get_coordinates_multiple_alts(fasta: Fasta):
    # rs4659982: T/A/C
    assert fasta.get_coordinates_for_clinical_annotation('rs4659982', 'NC_000001.11:240566955', ['TT', 'TA', 'AA']) == '1_240566955_T_A'
    # Use of different alt allele in genotypes generates a different identifier
    assert fasta.get_coordinates_for_clinical_annotation('rs4659982', 'NC_000001.11:240566955', ['TT', 'TC', 'CC']) == '1_240566955_T_C'
    # More than two alleles present in genotypes - TODO update once we decide how to treat these
    assert fasta.get_coordinates_for_clinical_annotation('rs4659982', 'NC_000001.11:240566955', ['TT', 'TA', 'TC', 'AC']) == '1_240566955_T_A'


def test_get_coordinates_deletion(fasta: Fasta):
    # rs70991108: -/TCGCGCGTCCCGCCCAGGT/TGGCGCCTCCCGCCCAGGT/TGGCGCGTCCCGCCCAGGT
    assert fasta.get_coordinates_for_clinical_annotation(
        'rs70991108',
        'NC_000005.10:80654345',
        ['TGGCGCGTCCCGCCCAGGT/TGGCGCGTCCCGCCCAGGT', 'TGGCGCGTCCCGCCCAGGT/del', 'del/del']
    ) == '5_80654344_C_CTGGCGCGTCCCGCCCAGGT'
