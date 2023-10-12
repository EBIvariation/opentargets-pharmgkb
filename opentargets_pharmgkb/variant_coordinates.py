import logging
from functools import lru_cache
import re

import pandas as pd
import requests
from Bio import SeqIO

logger = logging.getLogger(__package__)


@lru_cache
def get_chrom_pos_for_rs_from_ensembl(rsid):
    """Queries Ensembl for chromosome and position of rsid. Returns None if not found."""
    if not rsid.startswith('rs'):
        rsid = f'rs{rsid}'
    ensembl_url = f'https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json'
    resp = requests.get(ensembl_url)
    data = resp.json()
    if 'mappings' in data:
        for mapping in data['mappings']:
            if mapping['assembly_name'] == 'GRCh38':
                chrom = mapping['seq_region_name']
                # Skip things like CHR_HSCHR22_1_CTG7
                if '_' in chrom:
                    continue
                pos = mapping['start']
                return chrom, pos
    return None, None


def parse_genotype(genotype_string):
    """
    Parse PGKB string representations of genotypes into alleles.

    :param genotype_string: e.g. 'A', 'TA', 'A/del', 'CAG/CAGCAG'
    :return: list of alleles in the genotype, e.g. ['A'], ['T','A'], ['A','DEL'], ['CAG','CAGCAG']
    """
    alleles = []
    # X chrom variants
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
        # Correct for deletion alleles
        if 'DEL' in alleles_dict:
            if end == start:
                end -= 1  # keep end == start if they began that way
            start -= 1
            alleles_dict = {allele: self.add_context_base(chrom, start, allele) for allele in alleles_dict}

        if not alleles_dict:
            logger.warning(f'Could not parse any genotypes for {rsid}')
            return chrom_num, start, None, None

        ref = self.get_ref_from_fasta(chrom, start, end)
        # Report & skip if ref is not among the alleles
        # TODO we can now support cases where ref is truly not among the annotated genotypes, but we can't distinguish
        #  this situation from cases where the reference or location is wrong for some reason.
        if ref not in alleles_dict.values():
            logger.warning(f'Ref not in alleles: {rsid}\t{ref}\t{",".join(alleles_dict)}')
            return chrom_num, start, None, None

        return chrom_num, start, ref, alleles_dict

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
