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
        if pd.isna(location):
            return None
        chrom, pos = location.strip().split(':')
        if not chrom or not pos:
            return None

        # Ranges are inclusive of both start and end
        if '_' in pos:
            start, end = pos.split('_')
            start = int(start)
            end = int(end)
        else:
            start = end = int(pos)

        alleles = set()
        contains_del = False
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
                logger.info(f'Could not parse genotype for {rsid}: {genotype}')
                continue
            if m.group(1) == 'del'.lower() or m.group(2).lower() == 'del':
                contains_del = True
            alleles.add(m.group(1))
            alleles.add(m.group(2))

        # Correct for deletion alleles
        if contains_del:
            if end == start:
                end -= 1  # keep end == start if they began that way
            start -= 1
            alleles = {self.add_context_base(chrom, start, allele) for allele in alleles}

        if not alleles:
            logger.warning(f'Could not parse genotypes: {rsid}\t{",".join(all_genotypes)}')
            return None
        ref = self.get_ref_from_fasta(chrom, start, end)
        # Remove ref if present among alleles; otherwise report & skip
        if ref in alleles:
            alts = alleles - {ref}
        else:
            logger.warning(f'Ref not in alleles: {rsid}\t{ref}\t{"/".join(alleles)}')
            return None
        chrom_num = self.get_chrom_num_from_refseq(chrom)
        for alt in sorted(alts):
            # TODO multiple IDs for multiple alts?
            return f'{chrom_num}_{start}_{ref}_{alt}'
        return None

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
