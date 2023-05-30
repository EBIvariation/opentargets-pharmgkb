import json
import multiprocessing

import pandas as pd
from cmat.consequence_prediction.common.biomart import query_biomart

from opentargets_pharmgkb.ols import get_chebi_iri
from opentargets_pharmgkb.pandas_utils import none_to_nan, explode_column
from opentargets_pharmgkb.variant_coordinates import get_coordinates_for_clinical_annotation

ID_COL_NAME = 'Clinical Annotation ID'


def pipeline(clinical_annot_path, clinical_alleles_path, clinical_evidence_path, drugs_path, created_date, output_path):
    clinical_annot_table = pd.read_csv(clinical_annot_path, sep='\t')
    clinical_alleles_table = pd.read_csv(clinical_alleles_path, sep='\t')
    clinical_evidence_table = pd.read_csv(clinical_evidence_path, sep='\t')
    drugs_table = pd.read_csv(drugs_path, sep='\t')

    merged_table = pd.merge(clinical_annot_table, clinical_alleles_table, on=ID_COL_NAME, how='left')
    # Restrict to variants with rsIDs
    rs_only_table = merged_table[merged_table['Variant/Haplotypes'].str.contains('rs')]
    # Also provide a column with all genotypes for a given rs
    rs_only_table = pd.merge(rs_only_table, rs_only_table.groupby(by=ID_COL_NAME).aggregate(
        all_genotypes=('Genotype/Allele', list)), on=ID_COL_NAME)

    mapped_genes = explode_and_map_genes(rs_only_table)
    mapped_drugs = explode_and_map_drugs(mapped_genes, drugs_table)

    # Add clinical evidence with PMIDs
    pmid_evidence = clinical_evidence_table[clinical_evidence_table['PMID'].notna()].astype({'PMID': 'int'})
    evidence_table = pd.merge(mapped_drugs, pmid_evidence.groupby(by=ID_COL_NAME).aggregate(
        publications=('PMID', list)), on=ID_COL_NAME)

    # Generate evidence
    evidence = [generate_clinical_annotation_evidence(created_date, row) for _, row in evidence_table.iterrows()]
    with open(output_path, 'w+') as output:
        output.write('\n'.join(json.dumps(ev) for ev in evidence))


def explode_and_map_genes(df):
    """
    Maps gene symbols to Ensembl gene IDs using BioMart. Explodes multiple genes in single row.

    :param df: dataframe to annotate (should have a 'Gene' column)
    :return: dataframe with 'ensembl_gene_id' column added
    """
    split_genes = explode_column(df, 'Gene', 'split_gene')
    ensembl_ids = query_biomart(
        ('hgnc_symbol', 'split_gene'),
        ('ensembl_gene_id', 'ensembl_gene_id'),
        split_genes['split_gene'].drop_duplicates().tolist()
    )
    mapped_genes = pd.merge(split_genes, ensembl_ids, on='split_gene')
    # HGNC could map to more than one ensembl gene id, so must explode again
    mapped_genes = mapped_genes.explode('ensembl_gene_id').reset_index(drop=True)
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
    return f'http://purl.obolibrary.org/obo/CHEBI_{id_}'


def generate_clinical_annotation_evidence(created_date, row):
    """Generates an evidence string for a PharmGKB clinical annotation."""
    vcf_full_coords = get_coordinates_for_clinical_annotation(row['Variant/Haplotypes'], row['all_genotypes'])
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
        'variantId': vcf_full_coords,
        'variantRsId': row['Variant/Haplotypes'],
        'targetFromSourceId': row['ensembl_gene_id'],
        # TODO need to use consequence prediction from clinvar repo
        'variantFunctionalConsequenceId': None,
        'variantOverlappingGeneId': None,

        # GENOTYPE ATTRIBUTES
        'genotype': row['Genotype/Allele'],
        'genotypeAnnotationText': row['Annotation Text'],

        # PHENOTYPE ATTRIBUTES
        'drugText': row['split_drug'],
        'drugId': row['chebi'],
        'pgxCategory': row['Phenotype Category'],
        'phenotypeFromSourceId': row['Phenotype(s)']  # TODO EFO - needs to be exploded & mapped
    }
    # Remove the attributes with empty values (either None or empty lists).
    evidence_string = {key: value for key, value in evidence_string.items()
                       if value and (isinstance(value, list) or pd.notna(value))}
    return evidence_string
