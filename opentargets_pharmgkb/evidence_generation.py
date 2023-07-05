import json
import multiprocessing
import os
from itertools import zip_longest

import pandas as pd
from cmat.consequence_prediction.common.biomart import query_biomart
from cmat.consequence_prediction.snp_indel_variants.pipeline import process_variants

from opentargets_pharmgkb.counts import ClinicalAnnotationCounts
from opentargets_pharmgkb.ontology_apis import get_chebi_iri, get_efo_iri
from opentargets_pharmgkb.pandas_utils import none_to_nan, explode_column
from opentargets_pharmgkb.variant_coordinates import get_coordinates_for_clinical_annotation

ID_COL_NAME = 'Clinical Annotation ID'


def pipeline(data_dir, created_date, output_path):
    # TODO test file existence
    clinical_annot_path = os.path.join(data_dir, 'clinical_annotations.tsv')
    clinical_alleles_path = os.path.join(data_dir, 'clinical_ann_alleles.tsv')
    clinical_evidence_path = os.path.join(data_dir, 'clinical_ann_evidence.tsv')
    drugs_path = os.path.join(data_dir, 'drugs.tsv')

    clinical_annot_table = read_tsv_to_df(clinical_annot_path)
    clinical_alleles_table = read_tsv_to_df(clinical_alleles_path)
    clinical_evidence_table = read_tsv_to_df(clinical_evidence_path)
    drugs_table = read_tsv_to_df(drugs_path)

    # Restrict to variants with rsIDs
    rs_only_table = clinical_annot_table[clinical_annot_table['Variant/Haplotypes'].str.contains('rs')]

    # Gather input counts
    counts = ClinicalAnnotationCounts()
    counts.clinical_annotations = len(clinical_annot_table)
    counts.with_rs = len(rs_only_table)

    # Main processing
    merged_with_alleles_table = pd.merge(rs_only_table, clinical_alleles_table, on=ID_COL_NAME, how='left')
    coordinates_table = get_vcf_coordinates(merged_with_alleles_table)
    consequences_table = get_functional_consequences(coordinates_table)
    # mapped_genes = explode_and_map_genes(consequences_table)
    mapped_drugs = explode_and_map_drugs(consequences_table, drugs_table)
    mapped_phenotypes = explode_and_map_phenotypes(mapped_drugs)

    # Add clinical evidence with PMIDs
    pmid_evidence = clinical_evidence_table[clinical_evidence_table['PMID'].notna()]
    evidence_table = pd.merge(mapped_phenotypes, pmid_evidence.groupby(by=ID_COL_NAME).aggregate(
        publications=('PMID', list)), on=ID_COL_NAME)

    # Gather output counts
    counts.evidence_strings = len(evidence_table)
    counts.with_chebi = evidence_table['chebi'].count()
    counts.with_efo = evidence_table['efo'].count()
    counts.with_consequence = evidence_table['consequence_term'].count()
    # counts.with_pgkb_gene = evidence_table['gene_from_pgkb'].count()
    counts.with_vep_gene = evidence_table['overlapping_gene'].count()

    # Generate evidence
    evidence = [generate_clinical_annotation_evidence(created_date, row) for _, row in evidence_table.iterrows()]
    with open(output_path, 'w+') as output:
        output.write('\n'.join(json.dumps(ev) for ev in evidence))

    # Final count report
    gene_comparison_counts(evidence_table, counts, debug_path=f'{output_path}_genes.csv')
    counts.report()


def read_tsv_to_df(path):
    return pd.read_csv(path, sep='\t', dtype=str)


def get_vcf_coordinates(df):
    """
    Get VCF-style coordinates (chr_pos_ref_alt) for dataframe.

    :param df: dataframe to annotate (needs 'Genotype/Allele' and 'Variant/Haplotypes' columns)
    :return: dataframe with 'vcf_coords' column added
    """
    # First set a column with all genotypes for a given rs
    df_with_coords = pd.merge(df, df.groupby(by=ID_COL_NAME).aggregate(
        all_genotypes=('Genotype/Allele', list)), on=ID_COL_NAME)
    # Then get coordinates for each row
    for i, row in df_with_coords.iterrows():
        df_with_coords.at[i, 'vcf_coords'] = get_coordinates_for_clinical_annotation(
            row['Variant/Haplotypes'], row['all_genotypes'])
    return df_with_coords


def get_functional_consequences(df):
    """
    Get functional consequences from VEP.

    :param df: dataframe to annotate (needs 'vcf_coords' column)
    :return: dataframe with 'overlapping_gene' and 'consequence_term' columns added
    """
    vep_id_to_coords = {
        coord_id_to_vep_id(x): x for x in df['vcf_coords'].dropna().drop_duplicates().tolist()
    }
    with multiprocessing.Pool(processes=24) as pool:
        all_consequences = [
            pool.apply(process_to_list, args=(batch,))
            for batch in grouper(vep_id_to_coords.keys(), 200)
        ]
    mapped_consequences = pd.DataFrame(data=[
        {
            'vcf_coords': vep_id_to_coords[variant_id],
            'overlapping_gene': gene_id,
            'consequence_term': consequence_term
        }
        for batch in all_consequences
        for variant_id, gene_id, gene_symbol, consequence_term in batch
    ])
    return pd.merge(df, mapped_consequences, on='vcf_coords', how='left')


def coord_id_to_vep_id(coord_id):
    """Converts an underscore-separated coordinate identifier (e.g. 15_7237571_C_T) to VEP compatible one."""
    id_fields = coord_id.split('_')
    assert len(id_fields) == 4, 'Invalid identifier supplied (should contain exactly 4 fields)'
    return '{} {} . {} {}'.format(*id_fields)


def grouper(iterable, n):
    args = [iter(iterable)] * n
    return [x for x in zip_longest(*args, fillvalue=None) if x is not None]


def process_to_list(b):
    """Wrapper for process_variants because multiprocessing does not like generators."""
    return list(process_variants(b))


def explode_and_map_genes(df):
    """
    Maps gene symbols to Ensembl gene IDs using BioMart. Explodes multiple genes in single row.

    :param df: dataframe to annotate (should have a 'Gene' column)
    :return: dataframe with 'ensembl_gene_id' column added
    """
    split_genes = explode_column(df, 'Gene', 'split_gene')
    ensembl_ids = query_biomart(
        ('hgnc_symbol', 'split_gene'),
        ('ensembl_gene_id', 'gene_from_pgkb'),
        split_genes['split_gene'].drop_duplicates().tolist()
    )
    mapped_genes = pd.merge(split_genes, ensembl_ids, on='split_gene')
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


def generate_clinical_annotation_evidence(created_date, row):
    """Generates an evidence string for a PharmGKB clinical annotation."""
    evidence_string = {
        # DATA SOURCE ATTRIBUTES
        'datasourceId': 'pharmgkb',
        'datasourceVersion': created_date,

        # RECORD ATTRIBUTES
        'datatypeId': 'clinical_annotation',
        'studyId': row[ID_COL_NAME],
        'evidenceLevel': row['Level of Evidence'],
        'literature': [str(x) for x in row['publications']],

        # VARIANT ATTRIBUTES
        'variantId': row['vcf_coords'],
        'variantRsId': row['Variant/Haplotypes'],
        # 'originalSourceGeneId': row['gene_from_pgkb'],
        'variantFunctionalConsequenceId': row['consequence_term'],
        'targetFromSourceId': row['overlapping_gene'],

        # GENOTYPE ATTRIBUTES
        'genotype': row['Genotype/Allele'],
        'genotypeAnnotationText': row['Annotation Text'],

        # PHENOTYPE ATTRIBUTES
        'drugFromSource': row['split_drug'],
        'drugId': row['chebi'],
        'pgxCategory': row['Phenotype Category'],
        'phenotypeText': row['split_phenotype'],
        'phenotypeFromSourceId': row['efo']
    }
    # Remove the attributes with empty values (either None or empty lists).
    evidence_string = {key: value for key, value in evidence_string.items()
                       if value and (isinstance(value, list) or pd.notna(value))}
    return evidence_string


def gene_comparison_counts(df, counts, debug_path=None):
    # Map PGKB genes
    mapped_genes = explode_and_map_genes(df)
    # Re-group by ID column
    genes_table = mapped_genes.groupby(by=ID_COL_NAME).aggregate(
        all_pgkb_genes=('gene_from_pgkb', set),
        all_vep_genes=('overlapping_gene', set)
    )
    # Compare sets of genes
    counts.pgkb_vep_gene_diff = len(genes_table[genes_table['all_pgkb_genes'] != genes_table['all_vep_genes']])
    # Debug dump genes table
    if debug_path:
        genes_table.to_csv(debug_path)
