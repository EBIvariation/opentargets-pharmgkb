import logging
from functools import lru_cache
import re

import pandas as pd
from Bio import SeqIO

from opentargets_pharmgkb.ncbi_utils import get_spdi_coords_for_rsid

logger = logging.getLogger(__package__)


def parse_genotype(genotype_string):
    """
    Parse PGKB string representations of genotypes into alleles.

    :param genotype_string: e.g. 'A', 'TA', 'A/del', 'CAG/CAGCAG'
    :return: list of alleles in the genotype, e.g. ['A'], ['T','A'], ['A','DEL'], ['CAG','CAGCAG']
    """
    alleles = []
    # X/Y chrom variants
    if len(genotype_string) == 1:
        alleles.append(genotype_string)

    # SNPs
    if len(genotype_string) == 2:
        alleles.append(genotype_string[0])
        alleles.append(genotype_string[1])

    # short indels
    m = re.match('([ACGT]+|del)/([ACGT]+|del)', genotype_string, re.IGNORECASE)
    if m:
        alleles.append(m.group(1))
        alleles.append(m.group(2))

    if not alleles:
        logger.error(f'Could not parse genotype {genotype_string}')
    # Normalise to uppercase before returning
    return [a.upper() for a in alleles]


class Fasta:

    def __init__(self, path_to_fasta):
        self.record_dict = SeqIO.to_dict(SeqIO.parse(path_to_fasta, 'fasta'))

    def get_chr_pos_ref(self, rsid, location, all_parsed_genotypes):
        """
        Gets chromosome, position, and reference for a variant using location and genotype information.

        :param rsid: rsID of the variant
        :param location: string consisting of RefSeq accession and position separated by a colon,
                         e.g. 'NC_000001.11:46399999'
        :param all_parsed_genotypes: list of genotypes parsed into alleles, e.g. [['T', 'T'], ['T','A'], ['A','A']]
        :return: tuple of (chr, pos, ref, alleles), or Nones if coordinates cannot be determined.
                 alleles is a dict mapping allele string as present in all_parsed_genotypes, to an allele string with
                 context possibly added - e.g. {'AAG': 'CAAG', 'DEL': 'C'}
        """
        if pd.isna(location):
            return None, None, None, None
        chrom, pos = location.strip().split(':')
        if not chrom or not pos:
            return None, None, None, None
        chrom_num = self.get_chrom_num_from_refseq(chrom)

        # Ranges are inclusive of both start and end
        if '_' in pos:
            start, end = pos.split('_')
            start = int(start)
            end = int(end)
        else:
            start = end = int(pos)

        alleles_dict = {alt: alt for genotype in all_parsed_genotypes for alt in genotype}
        # Add context base for deletion alleles
        # This is the only normalisation needed for PGKB alleles, plus bookkeeping for the start/end coordinates
        if 'DEL' in alleles_dict:
            if end == start:
                end -= 1  # keep end == start if they began that way
            start -= 1
            alleles_dict = {allele: self.add_context_base(chrom, start, allele) for allele in alleles_dict}

        if not alleles_dict:
            logger.warning(f'Could not parse any genotypes for {rsid}')
            return chrom_num, start, None, None
        pgkb_alleles = alleles_dict.values()

        # First check for ref based on PGKB's coordinates
        ref = self.get_ref_from_fasta(chrom, start, end)
        if ref in pgkb_alleles:
            return chrom_num, start, ref, alleles_dict

        # If could not determine ref from FASTA, use ref determined from NCBI & normalised
        logger.info(f'Ref from FASTA not in alleles: {rsid}\t{ref}\t{",".join(pgkb_alleles)}')
        chrom, pos, ref = self.get_norm_coords_from_ncbi(rsid)
        logger.info(f'Will use ref from NCBI: {ref}')
        return chrom_num, pos, ref, alleles_dict

    def get_norm_coords_from_ncbi(self, rsid):
        """
        Get normalised coordinates from NCBI for an rsID.

        :param rsid: rsID to query
        :return: normalised coordinates (chrom, pos, ref)
        """
        chrom, pos, ref, alts = get_spdi_coords_for_rsid(rsid)
        chrom, norm_pos, norm_ref, norm_alts = self.normalise_with_ref(chrom, pos, ref, alts)
        # Don't actually need the alts from SPDI
        return chrom, norm_pos, norm_ref

    @lru_cache
    def get_ref_from_fasta(self, chrom, start, end=None):
        """Get reference allele at a given location."""
        if not end:
            end = start
        try:
            return str(self.record_dict[chrom][start-1:end].seq).upper()
        except (KeyError, IndexError):
            logger.warning(f'Could not get reference allele for {chrom}:{start}_{end}')
            return None

    @lru_cache
    def add_context_base(self, chrom, pos, allele):
        """Add left-hand context base to allele at a given location."""
        context_base = self.get_ref_from_fasta(chrom, pos)
        if context_base:
            if allele.lower() == 'del':
                return context_base
            return f'{context_base}{allele}'
        return None

    @lru_cache
    def get_chrom_num_from_refseq(self, chrom_refseq):
        # TODO use the assembly report? API call?
        m = re.search(r'Homo sapiens chromosome (.*?), ', self.record_dict[chrom_refseq].description)
        if m and m.group(1):
            return m.group(1)
        logger.warning(f'Could not get chromosome number for {chrom_refseq}')
        return None

    def normalise(self, chrom, pos, alleles):
        """
        Normalise alleles to be parsimonious and left-aligned.
        See here: https://genome.sph.umich.edu/wiki/Variant_Normalization

        :param chrom: chromosome (RefSeq)
        :param pos: position
        :param alleles: list of alleles to normalise
        :return: chromosome, normalised position, and list of normalised alleles (guaranteed to preserve input order)
        """
        # allow for initially empty alleles
        if any(len(a) == 0 for a in alleles):
            # extend alleles 1 to the left
            pos -= 1
            alleles = [self.add_context_base(chrom, pos, a) for a in alleles]
        # while all alleles end in same nucleotide
        while (len(set(a[-1] for a in alleles)) == 1):
            # truncate rightmost nucleotide
            alleles = [a[:-1] for a in alleles]
            # if exists an empty allele
            if any(len(a) == 0 for a in alleles):
                # extend alleles 1 to the left
                pos -= 1
                alleles = [self.add_context_base(chrom, pos, a) for a in alleles]
        # while all start with same nucleotide and have length 2 or more
        while (len(set(a[0] for a in alleles)) == 1) and all(len(a) >= 2 for a in alleles):
            # truncate leftmost nucleotide
            alleles = [a[1:] for a in alleles]
            pos += 1
        return chrom, pos, alleles

    def normalise_with_ref(self, chrom, pos, ref, alts):
        """Normalisation where reference allele is known."""
        chrom, pos, alleles = self.normalise(chrom, pos, [ref] + list(alts))
        return chrom, pos, alleles[0], alleles[1:]
