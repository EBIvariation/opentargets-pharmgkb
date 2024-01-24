import os

import numpy as np
import pandas as pd
import pytest

from opentargets_pharmgkb import evidence_generation
from opentargets_pharmgkb.evidence_generation import get_functional_consequences, explode_and_map_drugs, \
    read_tsv_to_df, explode_and_map_genes, get_genotype_ids, get_haplotype_ids
from tests.conftest import fasta_path

resources_dir = os.path.join(os.path.dirname(__file__), 'resources')


def test_get_genotype_ids():
    df = pd.DataFrame(columns=['Variant/Haplotypes', 'Location', 'Genotype/Allele'],
                      data=[['rs1051266', 'NC_000021.9:33341701', 'GG'],
                            ['rs1051266', 'NC_000021.9:33341701', 'GT'],
                            ['rs1051266', 'NC_000021.9:33341701', 'GA'],
                            ['rs1051266', 'NC_000021.9:33341701', 'TA']])
    annotated_df = get_genotype_ids(df, fasta_path)
    assert set(annotated_df['genotype_id'].values) == {
        '21_33341701_G_G,G',
        '21_33341701_G_G,T',
        '21_33341701_G_A,G',
        '21_33341701_G_A,T'}


def test_get_haplotype_ids():
    df = pd.DataFrame(columns=['Gene', 'Genotype/Allele'],
                      data=[['CYP2D6', '*1'],
                            ['CYP2D6', '*1xN'],
                            ['G6PD', 'A- 202A_376G'],
                            ['G6PD', 'Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham'],
                            ['GSTT1', 'null/non-null'],
                            [None, '*1'],
                            ['G6PD;GSTT1', '*1']])
    annotated_df = get_haplotype_ids(df)
    assert set(annotated_df['haplotype_id'].values) == {
        'CYP2D6*1',
        'CYP2D6*1xN',
        'G6PD A- 202A_376G',
        'G6PD Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham',
        'GSTT1 null/non-null'}


def test_get_functional_consequences():
    df = pd.DataFrame(columns=['genotype_id'], data=[['10_100980986_C_C,T']])
    annotated_df = get_functional_consequences(df)
    assert annotated_df.shape == (2, 3)
    assert 'ENSG00000095539' in annotated_df['overlapping_gene'].values
    assert 'intron_variant' in annotated_df['consequence_term'].values


def test_get_functional_consequences_ref_ref():
    df = pd.DataFrame(columns=['genotype_id'], data=[['10_100980986_C_C,C']])
    # If there is no variant-containing genotype we will get an error
    with pytest.raises(KeyError):
        get_functional_consequences(df)
    # Adding a variant-containing genotype yields 1 no_sequence_alteration row for each overlapping gene
    df = pd.concat((df, pd.DataFrame(columns=['genotype_id'], data=[['10_100980986_C_C,T']])))
    annotated_df = get_functional_consequences(df)
    assert annotated_df.shape == (4, 3)
    assert 'ENSG00000095539' in annotated_df['overlapping_gene'].values
    assert 'no_sequence_alteration' in annotated_df['consequence_term'].values


def test_explode_and_map_drugs():
    drugs_table = read_tsv_to_df(os.path.join(resources_dir, 'drugs.tsv'))
    df = pd.DataFrame(columns=['Drug(s)'], data=[['tamoxifen; fluorouracil'], ['peginterferon alfa-2a']])
    annotated_df = explode_and_map_drugs(df, drugs_table)
    assert annotated_df.shape[0] == 3
    annotated_df = annotated_df.set_index('split_drug')
    assert annotated_df.loc['tamoxifen']['chebi'] == 'http://purl.obolibrary.org/obo/CHEBI_41774'
    assert annotated_df.loc['fluorouracil']['chebi'] == 'http://purl.obolibrary.org/obo/CHEBI_46345'
    assert annotated_df.loc['peginterferon alfa-2a']['chebi'] is None


def test_explode_and_map_genes():
    df = pd.DataFrame(columns=['Gene'], data=[['IFNL3;IFNL4'], ['HLA-G'], [np.nan]])
    annotated_df = explode_and_map_genes(df)
    assert annotated_df.shape == (13, 3)
    assert 'ENSG00000272395' in annotated_df['gene_from_pgkb'].values
    assert 'ENSG00000235346' in annotated_df['gene_from_pgkb'].values
    assert pd.isna(annotated_df['gene_from_pgkb']).any()


def test_pipeline():
    output_path = os.path.join(resources_dir, 'test_output.json')
    expected_path = os.path.join(resources_dir, 'expected_output.json')
    evidence_generation.pipeline(
        data_dir=resources_dir,
        fasta_path=fasta_path,
        created_date='2023-03-23',
        output_path=output_path,
        debug_path=f'{output_path}.csv'
    )

    with open(output_path) as test_output, open(expected_path) as expected_output:
        assert sorted(test_output.readlines()) == sorted(expected_output.readlines())

    if os.path.exists(output_path):
        os.remove(output_path)
    if os.path.exists(f'{output_path}.csv'):
        os.remove(f'{output_path}.csv')


def test_pipeline_missing_file():
    output_path = os.path.join(resources_dir, 'test_output.json')
    with pytest.raises(ValueError):
        evidence_generation.pipeline(
            data_dir=os.path.join(resources_dir, 'nonexistent'),
            fasta_path=fasta_path,
            created_date='2023-03-23',
            output_path=output_path
        )
