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
from opentargets_pharmgkb.ontology_apis import get_efo_iri
from opentargets_pharmgkb.pandas_utils import none_to_nan, split_and_explode_column, read_tsv_to_df, nan_to_empty
from opentargets_pharmgkb.validation import validate_evidence_string
from opentargets_pharmgkb.variant_annotations import merge_variant_annotation_tables, get_variant_annotations, \
    DOE_COL_NAME, EFFECT_COL_NAME, OBJECT_COL_NAME, COMPARISON_COL_NAME, BASE_ALLELE_COL_NAME
from opentargets_pharmgkb.variant_coordinates import Fasta, parse_genotype

logging.basicConfig()
logger = logging.getLogger(__package__)
logger.setLevel(level=logging.DEBUG)

ID_COL_NAME = 'Clinical Annotation ID'
GENOTYPE_ALLELE_COL_NAME = 'Genotype/Allele'
VARIANT_HAPLOTYPE_COL_NAME = 'Variant/Haplotypes'


def pipeline(data_dir, fasta_path, created_date, output_path, with_doe=False):
    clinical_annot_path = os.path.join(data_dir, 'clinical_annotations.tsv')
    clinical_alleles_path = os.path.join(data_dir, 'clinical_ann_alleles.tsv')
    clinical_evidence_path = os.path.join(data_dir, 'clinical_ann_evidence.tsv')
    variants_path = os.path.join(data_dir, 'variants.tsv')
    relationships_path = os.path.join(data_dir, 'relationships.tsv')
    var_drug_path = os.path.join(data_dir, 'var_drug_ann.tsv')
    var_pheno_path = os.path.join(data_dir, 'var_pheno_ann.tsv')
    check_data_files_present((clinical_annot_path, clinical_alleles_path, clinical_evidence_path, variants_path))
    if with_doe:
        check_data_files_present((var_drug_path, var_pheno_path))

    clinical_annot_table = read_tsv_to_df(clinical_annot_path)
    clinical_alleles_table = read_tsv_to_df(clinical_alleles_path)
    clinical_evidence_table = read_tsv_to_df(clinical_evidence_path)
    variants_table = read_tsv_to_df(variants_path)
    relationships_table = read_tsv_to_df(relationships_path)
    if with_doe:
        unified_var_ann_table = merge_variant_annotation_tables(read_tsv_to_df(var_drug_path),
                                                                read_tsv_to_df(var_pheno_path))

    # Gather input counts
    # TODO think how to update counts
    counts = ClinicalAnnotationCounts()
    counts.clinical_annotations = len(clinical_annot_table)
    counts.with_rs = len(clinical_annot_table[clinical_annot_table[VARIANT_HAPLOTYPE_COL_NAME].str.startswith('rs')])

    # Main processing
    merged_with_variants_table = pd.merge(clinical_annot_table, variants_table, left_on=VARIANT_HAPLOTYPE_COL_NAME,
                                          right_on='Variant Name', how='left')
    merged_with_alleles_table = pd.merge(merged_with_variants_table, clinical_alleles_table, on=ID_COL_NAME, how='left')
    counts.exploded_alleles = len(merged_with_alleles_table)

    exploded_pgx_cat = split_and_explode_column(merged_with_alleles_table, 'Phenotype Category', 'split_pgx_category')
    counts.exploded_pgx_cat = len(exploded_pgx_cat)

    mapped_drugs = explode_drugs(exploded_pgx_cat)
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
        all_publications=('PMID', list)), on=ID_COL_NAME)
    if with_doe:
        parsed_var_ann_df = get_variant_annotations(evidence_table, pmid_evidence, unified_var_ann_table)
        evidence_table = pd.merge(evidence_table, parsed_var_ann_df, on=(ID_COL_NAME, GENOTYPE_ALLELE_COL_NAME))

    # Gather output counts
    counts.evidence_strings = len(evidence_table)
    counts.with_efo = evidence_table['efo'].count()
    counts.with_consequence = evidence_table['consequence_term'].count()
    counts.with_target_gene = evidence_table['overlapping_gene'].count() + evidence_table['gene_from_pgkb'].count()
    counts.with_haplotype = evidence_table['haplotype_id'].nunique()
    counts.resolved_haplotype_id = evidence_table['pgkb_haplotype_id'].nunique()

    # Generate evidence
    so_accession_dict = get_so_accession_dict()
    so_accession_dict['no_sequence_alteration'] = 'SO_0002073'
    evidence = [
        generate_clinical_annotation_evidence(so_accession_dict, created_date, row, with_doe)
        for _, row in evidence_table.iterrows()
    ]
    # Validate and write
    invalid_evidence = False
    with open(output_path, 'w+') as output:
        for ev_string in evidence:
            # DoE evidence will not validate for now
            if with_doe or validate_evidence_string(ev_string):
                output.write(json.dumps(ev_string)+'\n')
            else:
                invalid_evidence = True

    # Final count report
    counts.report()

    # Exit with an error code if any invalid evidence is produced
    # Do this at the very end so we still output counts and any valid evidence strings.
    if invalid_evidence:
        logger.error('Invalid evidence strings occurred, please check the logs for the details')
        sys.exit(1)


def check_data_files_present(paths):
    for p in paths:
        if not os.path.exists(p):
            logger.error(f'Missing required data file: {p}')
            raise ValueError(f'Missing required data file: {p}')


def split_df_with_or_without_rs(df):
    m = df[VARIANT_HAPLOTYPE_COL_NAME].str.startswith('rs')
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
    df_with_ids = df.assign(parsed_genotype=df[GENOTYPE_ALLELE_COL_NAME].apply(parse_genotype))
    df_with_ids = pd.merge(df_with_ids, df_with_ids.groupby(by=VARIANT_HAPLOTYPE_COL_NAME).aggregate(
        all_genotypes=('parsed_genotype', list)), on=VARIANT_HAPLOTYPE_COL_NAME)
    # Get coordinates (chromosome, position, reference, and all alternate alleles) for each RS
    rs_to_coords = {}
    for i, row in df_with_ids.drop_duplicates([VARIANT_HAPLOTYPE_COL_NAME]).iterrows():
        chrom, pos, ref, alleles_dict = fasta.get_chr_pos_ref(row[VARIANT_HAPLOTYPE_COL_NAME], row['Location'],
                                                              row['all_genotypes'])
        rs_to_coords[row[VARIANT_HAPLOTYPE_COL_NAME]] = (chrom, pos, ref, alleles_dict)
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
        chrom, pos, ref, alleles_dict = rs_to_coords[row[VARIANT_HAPLOTYPE_COL_NAME]]
        if chrom and pos and ref and alleles_dict and row['parsed_genotype']:
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
            + np.where(single_gene[GENOTYPE_ALLELE_COL_NAME].str.startswith('*'), '', ' ')
            # e.g. "*3" or "A- 202A_376G"
            + single_gene[GENOTYPE_ALLELE_COL_NAME]
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
    split_genes = split_and_explode_column(df, 'Gene', 'split_gene')
    ensembl_ids = query_biomart(
        ('hgnc_symbol', 'split_gene'),
        [('ensembl_gene_id', 'gene_from_pgkb')],
        split_genes['split_gene'].dropna().drop_duplicates().tolist()
    )
    mapped_genes = pd.merge(split_genes, ensembl_ids, on='split_gene', how='left')
    # HGNC could map to more than one ensembl gene id, so must explode again
    mapped_genes = mapped_genes.explode('gene_from_pgkb').reset_index(drop=True)
    return mapped_genes


def explode_drugs(df):
    """
    Explodes multiple drugs in a single row, unless they are known to be a drug combination.

    :param df: dataframe to annotate (should have a 'Drug(s)' column)
    :return: dataframe with 'split_drug' column added
    """
    # Drugs on same row but not explicitly annotated as combinations
    split_drugs = split_and_explode_column(df, 'Drug(s)', 'split_drug')
    # Drugs explicitly annotated as combinations are kept as a list of drug names
    split_drugs = split_and_explode_column(split_drugs, 'split_drug', 'split_drug', sep='/', split_only=True)
    return split_drugs


def explode_and_map_phenotypes(df):
    """
    Maps phenotype text to EFO IRIs using Zooma. Explodes multiple phenotypes in single row.

    :param df: dataframe to annotate (should have a 'Phenotype(s)' column)
    :return: dataframe with 'efo' column added
    """
    df['Phenotype(s)'].fillna('', inplace=True)
    split_phenotypes = split_and_explode_column(df, 'Phenotype(s)', 'split_phenotype')
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
    # Temporary workaround for MPATH terms until https://github.com/EBIvariation/CMAT/issues/417 is fixed
    return iri.split('/')[-1] if iri and pd.notna(iri) and 'MPATH' not in iri else None


def generate_clinical_annotation_evidence(so_accession_dict, created_date, row, with_doe):
    """Generates an evidence string for a PharmGKB clinical annotation."""
    partial_evidence_string = {
        # DATA SOURCE ATTRIBUTES
        'datasourceId': 'pharmgkb',
        'datasourceVersion': created_date,

        # RECORD ATTRIBUTES
        'datatypeId': 'clinical_annotation',
        'studyId': row[ID_COL_NAME],
        'evidenceLevel': row['Level of Evidence'],
        'literature': [str(x) for x in row['all_publications']],

        # GENOTYPE/ALLELE ATTRIBUTES
        'genotype': row[GENOTYPE_ALLELE_COL_NAME],
        'genotypeAnnotationText': row['Annotation Text'],
        'directionality': row['Allele Function'],

        # PHENOTYPE ATTRIBUTES
        'drugs': [{'drugFromSource': d} for d in row['split_drug']],
        'pgxCategory': row['split_pgx_category'].lower(),
        'phenotypeText': row['split_phenotype'],
        'phenotypeFromSourceId': iri_to_code(row['efo'])
    }
    evidence_string = add_variant_haplotype_attributes(so_accession_dict, row, partial_evidence_string)
    if with_doe:
        evidence_string = add_direction_of_effect_attributes(row, evidence_string)
    # Remove the attributes with empty values (either None or empty lists).
    evidence_string = {key: value for key, value in evidence_string.items()
                       if value and (isinstance(value, list) or pd.notna(value))}
    return evidence_string


def add_variant_haplotype_attributes(so_accession_dict, row, evidence_string):
    """Adds attributes to the evidence string depending on whether the record is for an rsID variant or not."""
    if row[VARIANT_HAPLOTYPE_COL_NAME].startswith('rs'):
        evidence_string.update({
            'genotypeId': row['genotype_id'],
            'variantRsId': row[VARIANT_HAPLOTYPE_COL_NAME],
            'variantFunctionalConsequenceId': so_accession_dict.get(row['consequence_term'], None),
            'targetFromSourceId': row['overlapping_gene'],
        })
    else:
        evidence_string.update({
            'haplotypeId': row['haplotype_id'],
            'haplotypeFromSourceId': row['pgkb_haplotype_id'],
            'targetFromSourceId': row['gene_from_pgkb']
        })
    return evidence_string


def add_direction_of_effect_attributes(row, evidence_string):
    evidence_string['evidenceFromSource'] = [
        {
            # Convert columns to a short summary statement, e.g. "increased metabolism of nicotine"
            'directionOfEffect': ' '.join((nan_to_empty(doe), nan_to_empty(effect), nan_to_empty(obj))).strip(),
            'baseAlleleOrGenotype': base_allele,
            'comparisonAlleleOrGenotype': comp_allele,
            'PMID': pmid,
            'annotationText': sentence
        }
        # Note pandas groupby().aggregate(list) will preserve order, so the use of zip is safe
        for pmid, doe, effect, obj, base_allele, comp_allele, sentence in zip(
            row['PMID'], row[DOE_COL_NAME], row[EFFECT_COL_NAME], row[OBJECT_COL_NAME],
            row[BASE_ALLELE_COL_NAME], row[COMPARISON_COL_NAME], row['Sentence']
        )
    ]
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
