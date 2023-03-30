import json

import pandas as pd

from opentargets_pharmgkb.variant_coordinates import get_coordinates_for_clinical_annotation

pgkb_column_names = [
    'Clinical Annotation ID',
    'Variant/Haplotypes',
    'Gene',
    'Level of Evidence',
    'Phenotype Category',
    'Drug(s)',
    'Phenotype(s)',
    'Genotype/Allele',
    'Annotation Text'
]


def pipeline(clinical_annot_path, clinical_alleles_path, created_date, output_path):
    clinical_annot_table = pd.read_csv(clinical_annot_path, sep='\t')
    clinical_alleles_table = pd.read_csv(clinical_alleles_path, sep='\t')

    merged_table = pd.merge(clinical_annot_table, clinical_alleles_table, on='Clinical Annotation ID', how='left')
    rs_only_table = merged_table[merged_table['Variant/Haplotypes'].str.contains('rs')][pgkb_column_names]

    # Also provide a column with all genotypes for a given rs
    rs_only_table = pd.merge(rs_only_table, rs_only_table.groupby(by='Clinical Annotation ID').aggregate(
        all_genotypes=('Genotype/Allele', list)), on='Clinical Annotation ID')

    # Generate evidence
    evidence = [generate_clinical_annotation_evidence(created_date, row) for _, row in rs_only_table.iterrows()]
    with open(output_path, 'w+') as output:
        output.write('\n'.join(json.dumps(ev) for ev in evidence))


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
        'targetFromSourceId': row['Gene'],  # TODO needs to be exploded & mapped
        # TODO need to use consequence prediction from clinvar repo
        'variantFunctionalConsequenceId': None,
        'variantOverlappingGeneId': None,

        # GENOTYPE ATTRIBUTES
        'genotype': row['Genotype/Allele'],
        'genotypeAnnotationText': row['Annotation Text'],

        # PHENOTYPE ATTRIBUTES
        'drugId': row['Drug(s)'],  # TODO CHEBI - needs to be exploded & mapped
        'pgxCategory': row['Phenotype Category'],
        'phenotypeFromSourceId': row['Phenotype(s)']  # TODO EFO - needs to be exploded & mapped
    }
    # Remove the attributes with empty values (either None or empty lists).
    evidence_string = {key: value for key, value in evidence_string.items() if value and pd.notna(value)}
    return evidence_string
