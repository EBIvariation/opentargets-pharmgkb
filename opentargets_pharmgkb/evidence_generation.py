import json
import multiprocessing

import pandas as pd
from cmat.consequence_prediction.common.biomart import query_biomart

from opentargets_pharmgkb.ols import get_chebi_iri
from opentargets_pharmgkb.variant_coordinates import get_coordinates_for_clinical_annotation


def pipeline(clinical_annot_path, clinical_alleles_path, created_date, output_path):
    clinical_annot_table = pd.read_csv(clinical_annot_path, sep='\t')
    clinical_alleles_table = pd.read_csv(clinical_alleles_path, sep='\t')

    merged_table = pd.merge(clinical_annot_table, clinical_alleles_table, on='Clinical Annotation ID', how='left')
    # Restrict to variants with rsIDs
    rs_only_table = merged_table[merged_table['Variant/Haplotypes'].str.contains('rs')]
    # Also provide a column with all genotypes for a given rs
    rs_only_table = pd.merge(rs_only_table, rs_only_table.groupby(by='Clinical Annotation ID').aggregate(
        all_genotypes=('Genotype/Allele', list)), on='Clinical Annotation ID')

    mapped_genes = explode_and_map_genes(rs_only_table)
    mapped_drugs = explode_and_map_drugs(mapped_genes)

    # Generate evidence
    evidence = [generate_clinical_annotation_evidence(created_date, row) for _, row in mapped_drugs.iterrows()]
    with open(output_path, 'w+') as output:
        output.write('\n'.join(json.dumps(ev) for ev in evidence))


def explode_and_map_genes(df):
    # Explode on the 'Gene' column (which consists of gene symbols) and convert each to Ensembl gene ID
    split_genes = df.assign(split_gene=df['Gene'].str.split(';')).explode('split_gene')
    ensembl_ids = query_biomart(
        ('hgnc_symbol', 'split_gene'),
        ('ensembl_gene_id', 'ensembl_gene_id'),
        split_genes['split_gene'].drop_duplicates().tolist()
    )
    mapped_genes = pd.merge(split_genes, ensembl_ids, on='split_gene')
    # HGNC could map to more than one ensembl gene id, so must explode again
    mapped_genes = mapped_genes.explode('ensembl_gene_id')
    return mapped_genes


def explode_and_map_drugs(df):
    # TODO compare this with using drugs.tsv directly
    # Explode on the 'Drug(s)' column and map to CHEBI iri in parallel
    split_drugs = df.assign(split_drug=df['Drug(s)'].str.split(';')).explode('split_drug')
    with multiprocessing.Pool(processes=24) as pool:
        str_to_iri = {
            s: pool.apply(get_chebi_iri, args=(s,))
            for s in split_drugs['split_drug'].drop_duplicates().tolist()
        }
    mapped_drugs = pd.concat(
        split_drugs[split_drugs['split_drug'] == s].assign(chebi=iri)
        for s, iri in str_to_iri.items()
    )
    return mapped_drugs


def generate_clinical_annotation_evidence(created_date, row):
    """Generates an evidence string for a PharmGKB clinical annotation."""
    vcf_full_coords = get_coordinates_for_clinical_annotation(row['Variant/Haplotypes'], row['all_genotypes'])
    evidence_string = {
        # DATA SOURCE ATTRIBUTES
        'datasourceId': 'pharmgkb',
        'datasourceVersion': created_date,

        # RECORD ATTRIBUTES
        'datatypeId': 'clinical_annotation',
        'studyId': row['Clinical Annotation ID'],
        'evidenceLevel': row['Level of Evidence'],

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
        'drugId': row['chebi'],
        'pgxCategory': row['Phenotype Category'],
        'phenotypeFromSourceId': row['Phenotype(s)']  # TODO EFO - needs to be exploded & mapped
    }
    # Remove the attributes with empty values (either None or empty lists).
    evidence_string = {key: value for key, value in evidence_string.items() if value and pd.notna(value)}
    return evidence_string
