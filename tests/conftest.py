import os

import pytest

from opentargets_pharmgkb.variant_coordinates import Fasta


@pytest.fixture
def fasta():
    fasta_path = os.path.join(os.path.dirname(__file__), 'resources', 'chr21.fa')
    return Fasta(fasta_path)
