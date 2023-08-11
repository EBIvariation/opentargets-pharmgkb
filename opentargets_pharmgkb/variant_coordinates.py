import logging
from functools import lru_cache
import re

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


class Fasta:

    def __init__(self, path_to_fasta):
        self.record_dict = SeqIO.to_dict(SeqIO.parse(path_to_fasta, 'fasta'))

    def get_coordinates_for_clinical_annotation(self, rsid, location, all_genotypes):
        """
        Gets vcf-style coordinate string (chr_pos_ref_alt) using location and genotype information.

        :param rsid: rsID of the variant
        :param location: string consisting of RefSeq accession and position separated by a colon,
                         e.g. 'NC_000001.11:46399999'
        :param all_genotypes: list of genotype strings, e.g. ['TT', 'TA', 'AA']
        :return: string of form chr_pos_ref_alt, or None if coordinates cannot be determined
        """
        chrom, pos = location.strip().split(':')
        if not chrom or not pos:
            # fall back on Ensembl (should be rare)
            chrom, pos = get_chrom_pos_for_rs_from_ensembl(rsid)
            if not chrom or not pos:
                return None

        # Ranges are inclusive of both start and end
        if '_' in pos:
            start, end = pos.split('_')
            start = int(start)
            end = int(end)
        else:
            start = end = int(pos)
        ref = self.get_ref_from_fasta(chrom, start, end)

        alleles = set()
        for genotype in all_genotypes:
            # X chrom variants
            if len(genotype) == 1:
                alleles.add(genotype)
                continue
            # SNPs
            if len(genotype) == 2:
                alleles.add(genotype[0])
                alleles.add(genotype[1])
                continue
            # short indels
            m = re.match('([ACGT]+|del)/([ACGT]+|del)', genotype, re.IGNORECASE)
            if not m:
                logger.info(f'Could not parse genotype: {genotype}')
                continue
            alleles.add(self.add_context_base(chrom, start, m.group(1)))
            alleles.add(self.add_context_base(chrom, start, m.group(2)))
            # In this case need to subtract 1 from start as we've added a base to the beginning
            start -= 1

        # Remove ref if present among alleles
        alts = alleles - {ref}
        for alt in alts:
            # TODO multiple IDs for multiple alts?
            return f'{chrom}_{start}_{ref}_{alt}'
        return None

    @lru_cache
    def get_ref_from_fasta(self, chrom, start, end=None):
        """Get reference allele at a given location."""
        if not end:
            end = start
        try:
            return self.record_dict[chrom][start-1:end].upper()
        except (KeyError, IndexError):
            logger.warning(f'Could not get reference allele for {chrom}:{start}_{end}')
            return None

    @lru_cache
    def add_context_base(self, chrom, start, allele):
        """Add left-hand context base to allele at a given location."""
        context_base = self.get_ref_from_fasta(chrom, start-1)
        if context_base:
            if allele.lower() == 'del':
                return context_base
            return f'{context_base}{allele}'
        return None
