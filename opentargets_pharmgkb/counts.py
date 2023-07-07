class ClinicalAnnotationCounts:
    """Simple class to hold counts for generating clinical annotation evidence."""

    def __init__(self):
        # Input counts (before exploding by allele, etc.)
        self.clinical_annotations = 0
        self.with_rs = 0
        # Output counts (after exploding)
        self.evidence_strings = 0
        self.with_chebi = 0
        self.with_efo = 0
        self.with_consequence = 0
        # self.with_pgkb_gene = 0
        self.with_vep_gene = 0
        # Evaluation counts - after annotation but without exploding
        self.annot_with_pgkb_genes = 0
        self.annot_with_vep_genes = 0
        self.pgkb_vep_gene_diff = 0

    def report(self):
        report_str = f'\nTotal clinical annotations: {self.clinical_annotations}\n'
        report_str += f'\tWith RS: {self.with_rs}\n'
        report_str += f'Total evidence strings: {self.evidence_strings}\n'
        report_str += f'\tWith CHEBI: {self.with_chebi}\n'
        report_str += f'\tWith EFO phenotype: {self.with_efo}\n'
        report_str += f'\tWith functional consequence: {self.with_consequence}\n'
        # report_str += f'\tWith PGKB gene: {self.with_pgkb_gene}\n'
        report_str += f'\tWith VEP gene: {self.with_vep_gene}\n'
        report_str += f'Gene comparisons per annotation\n'
        report_str += f'\tWith PGKB genes: {self.annot_with_pgkb_genes}\n'
        report_str += f'\tWith VEP genes: {self.annot_with_vep_genes}\n'
        report_str += f'\tPGKB genes != VEP genes: {self.pgkb_vep_gene_diff}\n'
        print(report_str)
        return report_str
