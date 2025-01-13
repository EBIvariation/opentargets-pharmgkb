import re

import pandas as pd

from opentargets_pharmgkb.pandas_utils import split_and_explode_column

ID_COL_NAME = 'Clinical Annotation ID'
VAR_ID_COL_NAME = 'Variant Annotation ID'
EFFECT_COL_NAME = 'effect_term'
OBJECT_COL_NAME = 'object_term'
ASSOC_COL_NAME = 'Is/Is Not associated'
DOE_COL_NAME = 'Direction of effect'
COMPARISON_COL_NAME = 'Comparison Allele(s) or Genotype(s)'


def merge_variant_annotation_tables(var_drug_table, var_pheno_table):
    """
    Return a single dataframe with selected columns from each variant annotation table.

    :param var_drug_table: variant drug annotation table
    :param var_pheno_table: variant phenotype annotation table
    :return: unified dataframe
    """
    # Select relevant columns
    # TODO confirm which columns we want to include
    drug_df = var_drug_table[[
        'Variant Annotation ID', 'PMID', 'Sentence', 'Alleles', 'Is/Is Not associated',
        'Direction of effect', 'PD/PK terms', 'Drug(s)',
        'Comparison Allele(s) or Genotype(s)'
    ]]
    phenotype_df = var_pheno_table[[
        'Variant Annotation ID', 'PMID', 'Sentence', 'Alleles', 'Is/Is Not associated',
        'Direction of effect', 'Side effect/efficacy/other', 'Phenotype',
        'Comparison Allele(s) or Genotype(s)'
    ]]
    # Rename differing columns so we can concat
    drug_df = drug_df.rename(columns={'PD/PK terms': EFFECT_COL_NAME, 'Drug(s)': OBJECT_COL_NAME})
    phenotype_df = phenotype_df.rename(columns={'Side effect/efficacy/other': EFFECT_COL_NAME,
                                                'Phenotype': OBJECT_COL_NAME})

    # TODO determine if we want to include functional evidence or not
    # functional_df = var_fa_ann[[
    #     'Variant Annotation ID', 'PMID', 'Sentence', 'Alleles', 'Is/Is Not associated',
    #     'Direction of effect', 'Functional terms', 'Gene/gene product',
    #     'Comparison Allele(s) or Genotype(s)'
    # ]]
    # functional_df = functional_df.rename(columns={'Functional terms': EFFECT_COL_NAME,
    #                                               'Gene/gene product': OBJECT_COL_NAME})

    return pd.concat((drug_df, phenotype_df))


# TODO modify to use dataframes directly rather than dict
def get_variant_annotations(clinical_alleles_df, clinical_evidence_df, variant_annotations):
    """Main method for getting resulting associations"""
    caid_to_vaid = {
        caid: clinical_evidence_df[clinical_evidence_df[ID_COL_NAME] == caid]['Evidence ID'].to_list()
        for caid in clinical_evidence_df[ID_COL_NAME]
    }
    results = {}
    for caid, vaids in caid_to_vaid.items():
        clinical_alleles_df = clinical_alleles_df[clinical_alleles_df[ID_COL_NAME] == caid][[ID_COL_NAME, 'Genotype/Allele', 'Annotation Text']]
        variant_ann_for_caid = variant_annotations[variant_annotations[VAR_ID_COL_NAME].isin(vaids)]
        results[caid] = get_associations(variant_ann_for_caid, clinical_alleles_df)
    return results


def get_associations(annotation_df, clinical_alleles_df):
    """

    :param annotation_df: variant annotation dataframe, already filtered to only contain one clinical annotation ID
    :param clinical_alleles_df: clinical annotation alleles dataframe, already filtered to contain one CAID
    :return:
    """
    # Split on +
    split_ann_df = split_and_explode_column(annotation_df, 'Alleles', 'split_alleles_1', sep='\+')
    # Split on /
    split_ann_df = split_and_explode_column(split_ann_df, 'split_alleles_1', 'split_alleles_2', sep='/')
    # Get alleles from clinical annotations - same logic as for getting ids
    split_clin_df = clinical_alleles_df.assign(
        parsed_genotype=clinical_alleles_df['Genotype/Allele'].apply(extended_parse_genotype))
    split_clin_df = split_clin_df.explode('parsed_genotype').reset_index(drop=True)

    # Match by +-split and /-split
    merged_df = pd.merge(split_clin_df, split_ann_df, how='outer', left_on='Genotype/Allele',
                         right_on='split_alleles_1')
    merged_df_2 = pd.merge(split_clin_df, split_ann_df, how='outer', left_on='parsed_genotype',
                           right_on='split_alleles_2')
    # TODO match also on comparison genotype/allele

    # If a genotype in a clinical annotation doesn't have evidence, want this listed with nan's
    all_results = []
    for _, genotype, _, parsed_genotype in split_clin_df.itertuples(index=False):
        # Rows that matched on genotype
        rows_first_match = merged_df[
            (merged_df['Genotype/Allele'] == genotype) & (merged_df['parsed_genotype'] == parsed_genotype) & (
                ~merged_df[VAR_ID_COL_NAME].isna())]
        # Rows that matched on parsed genotype
        rows_second_match = merged_df_2[
            (merged_df_2['Genotype/Allele'] == genotype) & (merged_df_2['parsed_genotype'] == parsed_genotype) & (
                ~merged_df_2[VAR_ID_COL_NAME].isna())]

        # If neither matches, add with nan's
        if rows_first_match.empty and rows_second_match.empty:
            all_results.append(merged_df[(merged_df['Genotype/Allele'] == genotype) & (
                        merged_df['parsed_genotype'] == parsed_genotype)])
        else:
            all_results.extend([rows_first_match, rows_second_match])

    final_result = pd.concat(all_results).drop_duplicates()

    # If _no_ part of a variant annotation is associated with any clinical annotation, want this listed with nan's
    for idx, row in split_ann_df.iterrows():
        vaid = row[VAR_ID_COL_NAME]
        alleles = row['Alleles']
        split_1 = row['split_alleles_1']
        split_2 = row['split_alleles_2']
        results_with_vaid = final_result[final_result[VAR_ID_COL_NAME] == vaid]
        if results_with_vaid.empty:
            final_result = pd.concat((final_result,
                                      merged_df[(merged_df[VAR_ID_COL_NAME] == vaid) & (
                                                  merged_df['Alleles'] == alleles) &
                                                (merged_df['split_alleles_1'] == split_1) & (
                                                            merged_df['split_alleles_2'] == split_2)]
                                      ))
    return final_result


def extended_parse_genotype(genotype_string):
    """
    Parse PGKB string representations of genotypes into alleles. Extended to include star alleles.
    TODO can we get rid of this method?
    """
    alleles = [genotype_string]

    # SNPs
    if len(genotype_string) == 2 and '*' not in genotype_string:
        alleles = [genotype_string[0], genotype_string[1]]

    # others
    m = re.match('([^/]+)/([^/]+)', genotype_string, re.IGNORECASE)
    if m:
        alleles = [m.group(1), m.group(2)]

    return alleles
