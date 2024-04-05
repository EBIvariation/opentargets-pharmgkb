from opentargets_pharmgkb.variant_coordinates import Fasta, parse_genotype


def test_parse_genotype():
    assert parse_genotype('TC') == ['T', 'C']
    assert parse_genotype('CAG/CAG') == ['CAG', 'CAG']
    assert parse_genotype('A/del') == ['A', 'DEL']


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


def test_get_coordinates_from_ncbi(fasta: Fasta):
    assert fasta.get_chr_pos_ref(
        'rs1051266',
        'NC_000021.9:45537879',
        [['TGA', 'TGA'], ['TGA', 'TGAGA'], ['TGAGA', 'TGAGA']]
    ) == ('21', 45537879, 'T', {'TGA': 'TGA', 'TGAGA': 'TGAGA'})
