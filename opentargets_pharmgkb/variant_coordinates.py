import logging
from functools import lru_cache
import re

import requests

logger = logging.getLogger(__package__)


def get_coordinates_for_clinical_annotation(rsid, all_genotypes):
    """
    Gets vcf-style coordinate string (chr_pos_ref_alt) using rsid alone, falling back on genotype information if needed.
    Returns None if coordinates cannot be determined.
    """
    chrom, pos, ref, alts = get_coordinates_for_rs(rsid)
    if not chrom or not pos or not ref or not alts:
        return None
    if len(alts) == 1:
        return f'{chrom}_{pos}_{ref}_{alts[0]}'
    # If multiple alts, check what is referred to in the clinical alleles table
    alleles = {ref}
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
            continue
        alleles.add(convert_allele(m.group(1)))
        alleles.add(convert_allele(m.group(2)))
    if len(alleles) > 2:
        logger.warning(f'Too many alleles for {rsid}, skipping')
        return None
    for a in alleles:
        if a in alts:
            return f'{chrom}_{pos}_{ref}_{a}'
    return None


def convert_allele(allele):
    # TODO "del" alleles are an issue for chr_pos_ref_alt - need context bases
    if allele.lower() == 'del':
        return '-'
    return allele


@lru_cache
def get_coordinates_for_rs(rsid):
    """Queries Ensembl for vcf-style coordinates for rsid. Returns None if not found."""
    # TODO do this in bulk (batches of 200) or replace with reference FASTA check
    if not rsid.startswith('rs'):
        rsid = f'rs{rsid}'
    ensembl_url = f'https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json'
    resp = requests.get(ensembl_url)
    data = resp.json()
    if 'mappings' in data:
        mapping = None
        for mapping in data['mappings']:
            if mapping['assembly_name'] == 'GRCh38':
                break
        if mapping:
            chrom = mapping['seq_region_name']
            pos = mapping['start']
            alleles = mapping['allele_string'].split('/')
            ref = alleles[0]
            alts = alleles[1:]
            return chrom, pos, ref, alts
    return None, None, None, None
