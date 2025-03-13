from unittest.mock import patch

from opentargets_pharmgkb.variant_coordinates import Fasta, parse_genotype


def test_parse_genotype():
    assert parse_genotype('TC') == ['T', 'C']
    assert parse_genotype('CAG/CAG') == ['CAG', 'CAG']
    assert parse_genotype('A/del') == ['A', 'DEL']
    assert parse_genotype('(CA)2/(CA)3') == ['CACA', 'CACACA']
    assert parse_genotype('(CCCACCCGA)2') == ['CCCACCCGACCCACCCGA']
    assert parse_genotype('something/(CCCACCCGA)2') == []


def test_get_coordinates(fasta: Fasta):
    assert fasta.get_chr_pos_ref(
        'rs1051266',
        'NC_000021.9:33341701',
        [['G', 'G'], ['G', 'T'], ['T', 'T']]) == ('21', 33341701, 'G', {'G': 'G', 'T': 'T'})


def test_get_coordinates_deletion(fasta: Fasta):
    assert fasta.get_chr_pos_ref(
        'rs1051266',
        'NC_000021.9:33341701',
        [['TGGCGCGTCCCGCCCAGGT', 'TGGCGCGTCCCGCCCAGGT'], ['TGGCGCGTCCCGCCCAGGT', 'DEL'], ['DEL', 'DEL']]
    ) == ('21', 33341700, 'A', {'DEL': 'A', 'TGGCGCGTCCCGCCCAGGT': 'ATGGCGCGTCCCGCCCAGGT'})


def test_get_coordinates_range_location(fasta: Fasta):
    assert fasta.get_chr_pos_ref(
        'rs1051266',
        'NC_000021.9:33341701_33341703',
        [['GAC', 'GAC'], ['GAC', 'DEL'], ['DEL', 'DEL']]
    ) == ('21', 33341700, 'AGAC', {'DEL': 'A', 'GAC': 'AGAC'})


def test_get_coordinates_not_match_reference_snp(fasta: Fasta):
    # Multi-allelic SNP, reference not annotated
    assert fasta.get_chr_pos_ref(
        'rs1051266',
        'NC_000021.9:45537880',
        [['C', 'C'], ['C', 'G'], ['G', 'G']]
    ) == ('21', 45537880, 'T', {'C': 'C', 'G': 'G'})


def test_get_coordinates_not_match_reference_del(fasta: Fasta):
    # Deletion
    # Mock NCBI response, so we can simulate an exact scenario on chr21 without searching for the right RS
    with patch('opentargets_pharmgkb.variant_coordinates.get_spdi_coords_for_rsid') as m_spdi_coords:
        m_spdi_coords.return_value = ('NC_000021.9', 45537879, 'TGCTGC', ['TGC'])
        result = fasta.get_chr_pos_ref(
            'rs1051266',
            'NC_000021.9:45537880_45537885',
            [['TGC', 'TGC'], ['TGC', 'DEL'], ['DEL', 'DEL']]
        )
        assert result == ('21', 45537879, 'GTGC', {'TGC': 'GTGC', 'DEL': 'G'})


def test_get_coordinates_not_match_reference_unnormalised_alleles(fasta: Fasta):
    # Test scenario where reference from FASTA is not found and PGKB alleles are not normalised
    # Mock NCBI response, so we can simulate an exact scenario on chr21 without searching for the right RS
    with patch('opentargets_pharmgkb.variant_coordinates.get_spdi_coords_for_rsid') as m_spdi_coords:
        # Definition of RS consists of adding 1 or 2 TTTA repeats at the beginning of the repeating section
        m_spdi_coords.return_value = ('NC_000021.9', 7678480, 'T', ['TTTTA', 'TTTTATTTA'])

        # Case 1: PGKB variant consists of inserting another copy of ATTT in the middle
        result_pos_same = fasta.get_chr_pos_ref(
            'rs1051266',
            'NC_000021.9:7678489',
            [['ATTT', 'ATTT'], ['ATTT', 'ATTTATTT'], ['ATTTATTT', 'ATTTATTT']]
        )
        # PGKB alleles normalise to the same position as NCBI gives us, so we use the normalised alleles
        assert result_pos_same == ('21', 7678481, 'T', {'ATTT': 'T', 'ATTTATTT': 'TTTTA'})

        # Case 2: PGKB variant consists of inserting ATT (not the repeat unit) in the middle
        result_pos_diff = fasta.get_chr_pos_ref(
            'rs1051266',
            'NC_000021.9:7678489',
            [['ATTT', 'ATTT'], ['ATTT', 'ATTTATT'], ['ATTTATT', 'ATTTATT']]
        )
        # PGKB alleles don't normalise to the same position, so we use the original alleles
        assert result_pos_diff == ('21', 7678481, 'T', {'ATTT': 'ATTT', 'ATTTATT': 'ATTTATT'})


def test_normalise(fasta: Fasta):
    # Section of chr21 consisting of ATTT repeats
    assert fasta.get_ref_from_fasta('NC_000021.9', 7678481, 7678480 + (10*4)) \
           == 'TTTTATTTATTTATTTATTTATTTATTTATTTATTTATTT'
    # Variant that deletes one repeat from the middle is normalised to be left-aligned
    assert fasta.normalise('NC_000021.9', 7678489, ['ATTTATTT', 'ATTT']) == ('NC_000021.9', 7678481, ['TTTTA', 'T'])
    # Same variant with different representation (empty allele)
    assert fasta.normalise('NC_000021.9', 7678489, ['ATTT', '']) == ('NC_000021.9', 7678481, ['TTTTA', 'T'])
    # Multiple alternate alleles
    assert fasta.normalise('NC_000021.9', 7678489, ['ATTTATTT', 'ATTT', 'ATTTATTTATTT']) \
           == ('NC_000021.9', 7678481, ['TTTTA', 'T', 'TTTTATTTA'])
