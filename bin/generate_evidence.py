#!/usr/bin/env python3

import argparse
from opentargets_pharmgkb import evidence_generation

parser = argparse.ArgumentParser('Generates Open Targets evidence strings from PharmGKB data')
parser.add_argument('--data-dir', help='Directory containing necessary .tsv files from PharmGKB', required=True)
parser.add_argument('--fasta', help='Path to FASTA file for GRCh38 (should use RefSeq contigs)')
parser.add_argument('--created-date', help='Created date of downloaded files (provided by PharmGKB)', required=True)
parser.add_argument('--output-path', help='Path to output evidence strings', required=True)
parser.add_argument('--with-doe', help='Generate evidence including direction of effect from variant annotations',
                    action='store_true', default=False)


if __name__ == '__main__':
    args = parser.parse_args()
    evidence_generation.pipeline(
        data_dir=args.data_dir,
        fasta_path=args.fasta,
        created_date=args.created_date,
        output_path=args.output_path,
        with_doe=args.with_doe
    )
