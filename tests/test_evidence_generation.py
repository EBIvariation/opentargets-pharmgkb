import os

from opentargets_pharmgkb import evidence_generation

resources_dir = os.path.join(os.path.dirname(__file__), 'resources')


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
