import json
import logging
import multiprocessing
import os
import sys
from collections import defaultdict
from itertools import zip_longest

import numpy as np
import pandas as pd
from cmat.consequence_prediction.common.biomart import query_biomart
from cmat.consequence_prediction.snp_indel_variants.pipeline import process_variants
from cmat.output_generation.consequence_type import get_so_accession_dict

from opentargets_pharmgkb.counts import ClinicalAnnotationCounts
from opentargets_pharmgkb.ontology_apis import get_chebi_iri, get_efo_iri
from opentargets_pharmgkb.pandas_utils import none_to_nan, explode_column, read_tsv_to_df
from opentargets_pharmgkb.validation import validate_evidence_string
from opentargets_pharmgkb.variant_coordinates import Fasta, parse_genotype

logging.basicConfig()
logger = logging.getLogger(__package__)
logger.setLevel(level=logging.DEBUG)

ID_COL_NAME = 'Clinical Annotation ID'


def pipeline(data_dir, fasta_path, created_date, output_path, debug_path=None):
    clinical_annot_path = os.path.join(data_dir, 'clinical_annotations.tsv')
    clinical_alleles_path = os.path.join(data_dir, 'clinical_ann_alleles.tsv')
    clinical_evidence_path = os.path.join(data_dir, 'clinical_ann_evidence.tsv')
    variants_path = os.path.join(data_dir, 'variants.tsv')
    drugs_path = os.path.join(data_dir, 'drugs.tsv')
    relationships_path = os.path.join(data_dir, 'relationships.tsv')
    for p in (clinical_annot_path, clinical_alleles_path, clinical_evidence_path, variants_path, drugs_path):
        if not os.path.exists(p):
            logger.error(f'Missing required data file: {p}')
            raise ValueError(f'Missing required data file: {p}')

    clinical_annot_table = read_tsv_to_df(clinical_annot_path)
    clinical_alleles_table = read_tsv_to_df(clinical_alleles_path)
    clinical_evidence_table = read_tsv_to_df(clinical_evidence_path)
    variants_table = read_tsv_to_df(variants_path)
    drugs_table = read_tsv_to_df(drugs_path)
    relationships_table = read_tsv_to_df(relationships_path)

    # Gather input counts
    counts = ClinicalAnnotationCounts()
    counts.clinical_annotations = len(clinical_annot_table)
    counts.with_rs = len(clinical_annot_table[clinical_annot_table['Variant/Haplotypes'].str.startswith('rs')])

    # Main processing
    merged_with_variants_table = pd.merge(clinical_annot_table, variants_table, left_on='Variant/Haplotypes',
                                          right_on='Variant Name', how='left')
    merged_with_alleles_table = pd.merge(merged_with_variants_table, clinical_alleles_table, on=ID_COL_NAME, how='left')
    counts.exploded_alleles = len(merged_with_alleles_table)

    exploded_pgx_cat = explode_column(merged_with_alleles_table, 'Phenotype Category', 'split_pgx_category')
    counts.exploded_pgx_cat = len(exploded_pgx_cat)

    mapped_drugs = explode_and_map_drugs(exploded_pgx_cat, drugs_table)
    counts.exploded_drugs = len(mapped_drugs)

    mapped_phenotypes = explode_and_map_phenotypes(mapped_drugs)
    counts.exploded_phenotypes = len(mapped_phenotypes)

    with_rs, no_rs = split_df_with_or_without_rs(mapped_phenotypes)

    # BRANCH 1: For rsIds, form genotype IDs + get consequences from VEP
    genotype_ids_table = get_genotype_ids(with_rs, fasta_path, counts)
    rs_consequences_table = get_functional_consequences(genotype_ids_table)

    # BRANCH 2: For named alleles, form haplotype IDs + map PGKB genes with Biomart
    haplotype_ids_table = get_haplotype_ids(no_rs, relationships_table)
    no_rs_genes_table = explode_and_map_genes(haplotype_ids_table)

    # Merge the tables - will fill in NaNs where columns are missing
    consequences_table = pd.concat((rs_consequences_table, no_rs_genes_table))

    # Add clinical evidence with PMIDs
    pmid_evidence = clinical_evidence_table[clinical_evidence_table['PMID'].notna()]
    evidence_table = pd.merge(consequences_table, pmid_evidence.groupby(by=ID_COL_NAME).aggregate(
        publications=('PMID', list)), on=ID_COL_NAME)

    # Gather output counts
    counts.evidence_strings = len(evidence_table)
    counts.with_chebi = evidence_table['chebi'].count()
    counts.with_efo = evidence_table['efo'].count()
    counts.with_consequence = evidence_table['consequence_term'].count()
    counts.with_target_gene = evidence_table['overlapping_gene'].count() + evidence_table['gene_from_pgkb'].count()
    counts.with_haplotype = evidence_table['haplotype_id'].nunique()
    counts.resolved_haplotype_id = evidence_table['pgkb_haplotype_id'].nunique()

    # Generate evidence
    so_accession_dict = get_so_accession_dict()
    so_accession_dict['no_sequence_alteration'] = 'SO_0002073'
    evidence = [
        generate_clinical_annotation_evidence(so_accession_dict, created_date, row)
        for _, row in evidence_table.iterrows()
    ]
    # Validate and write
    invalid_evidence = False
    with open(output_path, 'w+') as output:
        for ev_string in evidence:
            if True:  # TODO change back once schema fixed
            # if validate_evidence_string(ev_string):
                output.write(json.dumps(ev_string)+'\n')
            else:
                invalid_evidence = True

    # Final count report
    # NB. gene comparison conflicts with named allele processing,
    #  restore if needed for https://github.com/EBIvariation/opentargets-pharmgkb/issues/21
    # if not debug_path:
    #     debug_path = f'{output_path.rsplit(".", 1)[0]}_genes.csv'
    # gene_comparison_counts(evidence_table, counts, debug_path=debug_path)
    counts.report()

    # Exit with an error code if any invalid evidence is produced
    # Do this at the very end so we still output counts and any valid evidence strings.
    if invalid_evidence:
        logger.error('Invalid evidence strings occurred, please check the logs for the details')
        sys.exit(1)


def split_df_with_or_without_rs(df):
    m = df['Variant/Haplotypes'].str.startswith('rs')
    return df[m], df[~m]


def genotype_id(chrom, pos, ref, parsed_genotype):
    return f'{chrom}_{pos}_{ref}_{",".join(parsed_genotype)}' if chrom and pos and ref else None


def get_genotype_ids(df, fasta_path, counts=None):
    """
    Get genotype IDs (chr_pos_ref_allele1,allele2) for dataframe.

    :param df: dataframe to annotate (needs 'Genotype/Allele', 'Variant/Haplotypes', 'Location' columns)
    :param fasta_path: path to fasta file to check reference
    :param counts: ClinicalAnnotationCounts; if provided will count multi-allelic variants.
    :return: dataframe with 'genotype_id' column added
    """
    fasta = Fasta(fasta_path)
    # First set a column with all genotypes for a given RS
    df_with_ids = df.assign(parsed_genotype=df['Genotype/Allele'].apply(parse_genotype))
    df_with_ids = pd.merge(df_with_ids, df_with_ids.groupby(by='Variant/Haplotypes').aggregate(
        all_genotypes=('parsed_genotype', list)), on='Variant/Haplotypes')
    # Get coordinates (chromosome, position, reference, and all alternate alleles) for each RS
    rs_to_coords = {}
    for i, row in df_with_ids.drop_duplicates(['Variant/Haplotypes']).iterrows():
        chrom, pos, ref, alleles_dict = fasta.get_chr_pos_ref(row['Variant/Haplotypes'], row['Location'],
                                                              row['all_genotypes'])
        rs_to_coords[row['Variant/Haplotypes']] = (chrom, pos, ref, alleles_dict)
        # Generate per-variant counts, if applicable
        if not counts:
            continue
        counts.total_rs += 1
        if not alleles_dict:
            continue
        counts.rs_with_alleles += 1
        if len(alleles_dict) <= 2:
            continue
        counts.rs_with_more_than_2_alleles += 1
    # Use rs_to_coords to generate genotypeId for each genotype
    for i, row in df_with_ids.iterrows():
        chrom, pos, ref, alleles_dict = rs_to_coords[row['Variant/Haplotypes']]
        if chrom and pos and ref and alleles_dict:
            df_with_ids.at[i, 'genotype_id'] = genotype_id(chrom, pos, ref, sorted([alleles_dict[a]
                                                                                    for a in row['parsed_genotype']]))
        else:
            df_with_ids.at[i, 'genotype_id'] = None
    return df_with_ids


def get_functional_consequences(df):
    """
    Get functional consequences from VEP.

    :param df: dataframe to annotate (needs 'genotype_id' column)
    :return: dataframe with 'overlapping_gene' and 'consequence_term' columns added
    """
    vep_id_to_genotype_ids = defaultdict(list)
    for genotype_id in df['genotype_id'].dropna().drop_duplicates().tolist():
        for vep_id in genotype_id_to_vep_ids(genotype_id):
            vep_id_to_genotype_ids[vep_id].append(genotype_id)
    # Note that variants in a single genotype will have VEP logic applied independently, i.e. most severe consequence
    # for each overlapping gene.
    with multiprocessing.Pool(processes=24) as pool:
        all_consequences = [
            pool.apply(process_to_list, args=(batch,))
            for batch in grouper(vep_id_to_genotype_ids.keys(), 200)
        ]
    mapped_consequences = pd.DataFrame(data=[
        {
            'genotype_id': genotype_id,
            'overlapping_gene': gene_id,
            'consequence_term': consequence_term
        }
        for batch in all_consequences
        for variant_id, gene_id, gene_symbol, consequence_term in batch
        for genotype_id in vep_id_to_genotype_ids[variant_id]
    ]).drop_duplicates()
    # For every VEP id, also record the no_sequence_alteration consequence for the corresponding ref/ref genotype
    ref_ref_consequences = pd.DataFrame(data=[
        {
            'genotype_id': vep_id_to_ref_ref_id(variant_id),
            'overlapping_gene': gene_id,
            'consequence_term': 'no_sequence_alteration'
        }
        for batch in all_consequences
        for variant_id, gene_id, _, _ in batch
    ]).drop_duplicates()
    mapped_consequences = pd.concat((mapped_consequences, ref_ref_consequences))
    return pd.merge(df, mapped_consequences, on='genotype_id', how='left')


def genotype_id_to_vep_ids(coord_id):
    """Converts an underscore-separated genotype identifier (e.g. 15_7237571_C_T,C) to VEP compatible ones."""
    id_fields = coord_id.split('_')
    assert len(id_fields) == 4, 'Invalid identifier supplied (should contain exactly 4 fields)'
    chrom, pos, ref, genotype = id_fields
    genotype = genotype.split(',')
    for alt in genotype:
        # Skip non-variants
        if alt != ref:
            yield f'{chrom} {pos} . {ref} {alt}'


def vep_id_to_ref_ref_id(vep_id):
    """Converts a VEP compatible variant ID to an underscore-separated ref/ref genotype identifier."""
    chrom, pos, iden, ref, alt = vep_id.split(' ')
    return f'{chrom}_{pos}_{ref}_{ref},{ref}'


def grouper(iterable, n):
    args = [iter(iterable)] * n
    return [x for x in zip_longest(*args, fillvalue=None) if x is not None]


def process_to_list(b):
    """Wrapper for process_variants because multiprocessing does not like generators."""
    return list(process_variants(b, False))


def get_haplotype_ids(df, relationships_table):
    """
    Get haplotype IDs (gene symbol + allele name) for dataframe.

    :param df: dataframe to annotate (needs 'Gene' and 'Genotype/Allele' columns)
    :param relationships_table: table of entity relationships in PGKB
    :return: dataframe with 'haplotype_id' and 'pgkb_haplotype_id' columns added
    """
    # Filter out rows where gene column is missing or contains multiple genes.
    # This is just a safeguard, we really shouldn't have these.
    single_gene = df[~df['Gene'].str.contains(';', na=True)]
    if len(single_gene) != len(df):
        logger.warning(
            f'Tried to get haplotype IDs for {len(df)-len(single_gene)} rows with multiple or missing genes!'
        )
    # Construct the haplotype ID
    single_gene['haplotype_id'] = (
            # e.g. "CYP2D6" or "G6PD"
            single_gene['Gene']
            # whether to concatenate with a space or not
            + np.where(single_gene['Genotype/Allele'].str.startswith('*'), '', ' ')
            # e.g. "*3" or "A- 202A_376G"
            + single_gene['Genotype/Allele']
    )
    # Fetch the internal haplotype ID
    single_gene['pgkb_haplotype_id'] = single_gene['haplotype_id'].apply(hap_id_to_pgkb_id, args=(relationships_table,))

    return single_gene


def hap_id_to_pgkb_id(id_, relationships_table):
    if pd.notna(id_):
        vals = set(relationships_table[relationships_table['Entity1_name'] == id_]['Entity1_id'])
        if len(vals) != 1:
            logger.warning(f'Could not determine internal ID for haplotype {id_}')
        else:
            return vals.pop()
    return None


def explode_and_map_genes(df):
    """
    Maps gene symbols to Ensembl gene IDs using BioMart. Explodes multiple genes in single row.

    :param df: dataframe to annotate (should have a 'Gene' column)
    :return: dataframe with 'ensembl_gene_id' column added
    """
    split_genes = explode_column(df, 'Gene', 'split_gene')
    ensembl_ids = query_biomart(
        ('hgnc_symbol', 'split_gene'),
        [('ensembl_gene_id', 'gene_from_pgkb')],
        split_genes['split_gene'].dropna().drop_duplicates().tolist()
    )
    mapped_genes = pd.merge(split_genes, ensembl_ids, on='split_gene', how='left')
    # HGNC could map to more than one ensembl gene id, so must explode again
    mapped_genes = mapped_genes.explode('gene_from_pgkb').reset_index(drop=True)
    return mapped_genes


def explode_and_map_drugs(df, drugs_table):
    """
    Maps drug names to CHEBI IRIs using OLS, falling back on primary drugs data from PharmGKB if needed.
    Explodes multiple drugs in a single row.

    :param df: dataframe to annotate (should have a 'Drug(s)' column)
    :param drugs_table: drugs dataframe
    :return: dataframe with 'chebi' column added
    """
    split_drugs = explode_column(df, 'Drug(s)', 'split_drug')
    # Query OLS in parallel as there are no batch queries currently.
    with multiprocessing.Pool(processes=24) as pool:
        str_to_iri = {
            s: pool.apply(get_chebi_iri, args=(s,))
            for s in split_drugs['split_drug'].drop_duplicates().tolist()
        }
    mapped_drugs = pd.concat(
        split_drugs[split_drugs['split_drug'] == s].assign(chebi=none_to_nan(iri))
        for s, iri in str_to_iri.items()
    )
    # Some drugs we can't unambiguously map using OLS, so we rely on primary data provided by PharmGKB.
    # Using OLS first ensures we get up-to-date IDs wherever possible.
    drugs_table['chebi_id'] = drugs_table['Cross-references'].str.extract(r'CHEBI:(?P<chebi_id>\d+)', expand=False)
    mapped_drugs = pd.merge(mapped_drugs, drugs_table, left_on='split_drug', right_on='Name', how='left')
    mapped_drugs['chebi'].fillna(mapped_drugs['chebi_id'].apply(chebi_id_to_iri), inplace=True)
    return mapped_drugs


def chebi_id_to_iri(id_):
    if pd.notna(id_):
        return f'http://purl.obolibrary.org/obo/CHEBI_{id_}'
    return None


def explode_and_map_phenotypes(df):
    """
    Maps phenotype text to EFO IRIs using Zooma. Explodes multiple phenotypes in single row.

    :param df: dataframe to annotate (should have a 'Phenotype(s)' column)
    :return: dataframe with 'efo' column added
    """
    df['Phenotype(s)'].fillna('', inplace=True)
    split_phenotypes = explode_column(df, 'Phenotype(s)', 'split_phenotype')
    with multiprocessing.Pool(processes=24) as pool:
        str_to_iri = {
            s: pool.apply(get_efo_iri, args=(s,))
            for s in split_phenotypes['split_phenotype'].drop_duplicates().tolist()
        }
    mapped_phenotypes = pd.concat(
        split_phenotypes[split_phenotypes['split_phenotype'] == s].assign(efo=none_to_nan(iri))
        for s, iri in str_to_iri.items()
    )
    return mapped_phenotypes


def iri_to_code(iri):
    """Convert iri (e.g. http://purl.obolibrary.org/obo/CHEBI_4792) to code, per Open Targets request."""
    return iri.split('/')[-1] if iri and pd.notna(iri) else None


def generate_clinical_annotation_evidence(so_accession_dict, created_date, row):
    """Generates an evidence string for a PharmGKB clinical annotation."""
    partial_evidence_string = {
        # DATA SOURCE ATTRIBUTES
        'datasourceId': 'pharmgkb',
        'datasourceVersion': created_date,

        # RECORD ATTRIBUTES
        'datatypeId': 'clinical_annotation',
        'studyId': row[ID_COL_NAME],
        'evidenceLevel': row['Level of Evidence'],
        'literature': [str(x) for x in row['publications']],

        # GENOTYPE/ALLELE ATTRIBUTES
        'genotype': row['Genotype/Allele'],
        'genotypeAnnotationText': row['Annotation Text'],
        'alleleFunction': row['Allele Function'],

        # PHENOTYPE ATTRIBUTES
        'drugFromSource': row['split_drug'],
        'drugId': iri_to_code(row['chebi']),
        'pgxCategory': row['split_pgx_category'].lower(),
        'phenotypeText': row['split_phenotype'],
        'phenotypeFromSourceId': iri_to_code(row['efo'])
    }
    evidence_string = add_variant_haplotype_attributes(so_accession_dict, row, partial_evidence_string)
    # Remove the attributes with empty values (either None or empty lists).
    evidence_string = {key: value for key, value in evidence_string.items()
                       if value and (isinstance(value, list) or pd.notna(value))}
    return evidence_string


def add_variant_haplotype_attributes(so_accession_dict, row, evidence_string):
    """Adds attributes to the evidence string depending on whether the record is for an rsID variant or not."""
    if row['Variant/Haplotypes'].startswith('rs'):
        evidence_string.update({
            'genotypeId': row['genotype_id'],
            'variantRsId': row['Variant/Haplotypes'],
            'variantFunctionalConsequenceId': so_accession_dict.get(row['consequence_term'], None),
            'targetFromSourceId': row['overlapping_gene'],
        })
    else:
        evidence_string.update({
            'haplotypeId': row['haplotype_id'],
            'internalHaplotypeId': row['pgkb_haplotype_id'],
            'targetFromSourceId': row['gene_from_pgkb']
        })
    return evidence_string


def gene_comparison_counts(df, counts, debug_path=None):
    # Map PGKB genes
    mapped_genes = explode_and_map_genes(df)
    # Re-group by ID column
    genes_table = mapped_genes.groupby(by=ID_COL_NAME).aggregate(
        all_pgkb_genes=('gene_from_pgkb', lambda x: set(x.dropna())),
        all_vep_genes=('overlapping_gene', lambda x: set(x.dropna()))
    )
    # Compare sets of genes
    counts.annot_with_pgkb_genes = len(genes_table[genes_table['all_pgkb_genes'] != set()])
    counts.annot_with_vep_genes = len(genes_table[genes_table['all_vep_genes'] != set()])
    neq_genes_table = genes_table[genes_table['all_pgkb_genes'] != genes_table['all_vep_genes']]
    counts.pgkb_vep_gene_diff = len(neq_genes_table)
    # Debug dump genes table
    if debug_path:
        neq_genes_table.to_csv(debug_path)
