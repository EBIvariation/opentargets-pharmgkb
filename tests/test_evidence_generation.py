import os

import pandas as pd

from opentargets_pharmgkb import evidence_generation
from opentargets_pharmgkb.evidence_generation import get_functional_consequences, explode_and_map_drugs, read_tsv_to_df

resources_dir = os.path.join(os.path.dirname(__file__), 'resources')


def test_get_functional_consequences():
    df = pd.DataFrame(columns=['vcf_coords'], data=[['10_100980986_C_T']])
    annotated_df = get_functional_consequences(df)
    assert annotated_df.shape == (2, 3)
    assert 'ENSG00000095539' in annotated_df['overlapping_gene'].values
    assert 'intron_variant' in annotated_df['consequence_term'].values


def test_explode_and_map_drugs():
    drugs_table = read_tsv_to_df(os.path.join(resources_dir, 'drugs.tsv'))
    df = pd.DataFrame(columns=['Drug(s)'], data=[['tamoxifen; fluorouracil'], ['peginterferon alfa-2a']])
    annotated_df = explode_and_map_drugs(df, drugs_table)
    assert annotated_df.shape[0] == 3
    annotated_df = annotated_df.set_index('split_drug')
    assert annotated_df.loc['tamoxifen']['chebi'] == 'http://purl.obolibrary.org/obo/CHEBI_41774'
    assert annotated_df.loc['fluorouracil']['chebi'] == 'http://purl.obolibrary.org/obo/CHEBI_46345'
    assert annotated_df.loc['peginterferon alfa-2a']['chebi'] is None


def test_pipeline():
    output_path = os.path.join(resources_dir, 'test_output.json')
    expected_path = os.path.join(resources_dir, 'expected_output.json')
    evidence_generation.pipeline(
        clinical_annot_path=os.path.join(resources_dir, 'clinical_annotations.tsv'),
        clinical_alleles_path=os.path.join(resources_dir, 'clinical_ann_alleles.tsv'),
        clinical_evidence_path=os.path.join(resources_dir, 'clinical_ann_evidence.tsv'),
        drugs_path=os.path.join(resources_dir, 'drugs.tsv'),
        created_date='2023-03-23',
        output_path=output_path
    )

    with open(output_path) as test_output, open(expected_path) as expected_output:
        assert test_output.readlines() == expected_output.readlines()

    if os.path.exists(output_path):
        os.remove(output_path)
