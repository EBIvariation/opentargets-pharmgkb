import os

import pytest

from opentargets_pharmgkb.variant_coordinates import Fasta


fasta_path = os.path.join(os.path.dirname(__file__), 'resources', 'chr21.fa')


@pytest.fixture
def fasta():
    return Fasta(fasta_path)
