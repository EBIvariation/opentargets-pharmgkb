#!/usr/bin/env python3

import argparse
from opentargets_pharmgkb import evidence_generation

parser = argparse.ArgumentParser('Generates Open Targets evidence strings from PharmGKB data')
parser.add_argument('--clinical-annot-path', help='Path to clinical_annotations.tsv', required=True)
parser.add_argument('--clinical-alleles-path', help='Path to clinical_ann_alleles.tsv', required=True)
parser.add_argument('--drugs-path', help='Path to drugs.tsv', required=True)
parser.add_argument('--created-date', help='Created date of downloaded files (provided by PharmGKB)', required=True)
parser.add_argument('--output-path', help='Path to output evidence strings', required=True)


if __name__ == '__main__':
    args = parser.parse_args()
    evidence_generation.pipeline(
        clinical_annot_path=args.clinical_annot_path,
        clinical_alleles_path=args.clinical_alleles_path,
        drugs_path=args.drugs_path,
        created_date=args.created_date,
        output_path=args.output_path
    )
