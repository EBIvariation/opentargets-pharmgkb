import numpy as np
import pandas as pd

from opentargets_pharmgkb.variant_annotations import get_variant_annotations, associate_annotations_with_alleles


def test_associate_annotations_with_alleles_snp():
    # Tests SNP genotype parsing and associating variant annotations with clinical annotation genotypes
    clinical_alleles_df = pd.DataFrame([
        ['1', 'AA'],
        ['1', 'AG'],
        ['1', 'GG']
    ], columns=['Clinical Annotation ID', 'Genotype/Allele'])
    var_annotations_df = pd.DataFrame([
        ['1.1', 'AA + AG'],
        ['1.2', 'AA'],
        ['1.3', 'A'],
        ['1.4', 'TT']
    ], columns=['Variant Annotation ID', 'Alleles'])
    expected_df = pd.DataFrame([
        ['1', 'AA', '1.1', 'AA + AG'],
        ['1', 'AA', '1.2', 'AA'],
        ['1', 'AA', '1.3', 'A'],
        ['1', 'AG', '1.1', 'AA + AG'],
        ['1', 'AG', '1.3', 'A'],
        ['1', 'GG', np.nan, np.nan]
    ], columns=['Clinical Annotation ID', 'Genotype/Allele', 'Variant Annotation ID', 'Alleles'])
    result = associate_annotations_with_alleles(var_annotations_df, clinical_alleles_df)
    assert result[[
        'Clinical Annotation ID', 'Genotype/Allele', 'Variant Annotation ID', 'Alleles'
    ]].reset_index(drop=True).equals(expected_df)


def test_associate_annotations_with_alleles_star():
    # Tests star allele parsing and associating variant annotations with clinical annotation alleles
    clinical_alleles_df = pd.DataFrame([
        ['2', '*1'],
        ['2', '*2'],
        ['2', '*3'],
        ['2', '*4']
    ], columns=['Clinical Annotation ID', 'Genotype/Allele'])
    var_annotations_df = pd.DataFrame([
        ['2.1', '*1'],
        ['2.2', '*1/*2'],
        ['2.3', '*1/*2 + *1/*3'],
        ['2.4', '*5']
    ], columns=['Variant Annotation ID', 'Alleles'])
    expected_df = pd.DataFrame([
        ['2', '*1', '2.1', '*1'],
        ['2', '*1', '2.2', '*1/*2'],
        ['2', '*1', '2.3', '*1/*2 + *1/*3'],
        ['2', '*2', '2.2', '*1/*2'],
        ['2', '*2', '2.3', '*1/*2 + *1/*3'],
        ['2', '*3', '2.3', '*1/*2 + *1/*3'],
        ['2', '*4', np.nan, np.nan]
    ], columns=['Clinical Annotation ID', 'Genotype/Allele', 'Variant Annotation ID', 'Alleles'])
    result = associate_annotations_with_alleles(var_annotations_df, clinical_alleles_df)
    assert result[[
        'Clinical Annotation ID', 'Genotype/Allele', 'Variant Annotation ID', 'Alleles'
    ]].reset_index(drop=True).equals(expected_df)


def test_get_variant_annotations():
    clinical_alleles_df = pd.DataFrame([
        ['1', 'AA'],
        ['1', 'AG'],
        ['1', 'GG'],
        ['2', '*1'],
        ['2', '*2']
    ], columns=['Clinical Annotation ID', 'Genotype/Allele'])
    var_annotations_df = pd.DataFrame([
        ['1.1', 'AA + AG', '123', 'sentence', 'Associated with', 'increased', 'metabolism of', 'nicotine', np.nan],
        ['1.2', 'AA', '456', 'sentence', 'Associated with', 'decreased', 'response to', 'morphine', 'AG + GG'],
        ['2.1', '*1', '789', 'sentence', 'Associated with', 'increased', 'risk of', 'vomiting', '*2'],
        ['2.2', '*1/*2', '000', 'sentence', 'Not associated with', 'increased', 'risk of', 'vomiting', '*2/*2']
    ], columns=['Variant Annotation ID', 'Alleles', 'PMID', 'Sentence', 'Is/Is Not associated', 'Direction of effect',
                'effect_term', 'object_term', 'Comparison Allele(s) or Genotype(s)'])
    clinical_evidence_df = pd.DataFrame([
        ['1', '1.1'],
        ['1', '1.2'],
        ['2', '2.1'],
        ['2', '2.2']
    ], columns=['Clinical Annotation ID', 'Evidence ID'])
    expected_df = pd.DataFrame([
        ['AA', ['123', '456'], ['sentence', 'sentence'], ['AA + AG', 'AA'], ['increased', 'decreased'],
         ['metabolism of', 'response to'], ['nicotine', 'morphine'], [np.nan, 'AG + GG'], '1'],
        ['AG', ['123'], ['sentence'], ['AA + AG'], ['increased'], ['metabolism of'], ['nicotine'], [np.nan], '1'],
        ['GG', [], [], [], [], [], [], [], '1'],
        ['*1', ['789'], ['sentence'], ['*1'], ['increased'], ['risk of'], ['vomiting'], ['*2'], '2'],
        ['*2', [], [], [], [], [], [], [], '2']
    ], columns=['Genotype/Allele', 'PMID', 'Sentence', 'Alleles', 'Direction of effect', 'effect_term', 'object_term',
                'Comparison Allele(s) or Genotype(s)', 'Clinical Annotation ID'])
    result = get_variant_annotations(clinical_alleles_df, clinical_evidence_df, var_annotations_df)
    assert result.reset_index(drop=True).equals(expected_df)
