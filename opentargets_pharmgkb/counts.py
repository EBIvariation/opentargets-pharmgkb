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

    def report(self):
        report_str = f'\nTotal clinical annotations: {self.clinical_annotations}\n'
        report_str += f'\tWith RS: {self.with_rs}\n'
        report_str += f'Total evidence strings: {self.evidence_strings}\n'
        report_str += f'\tWith CHEBI: {self.with_chebi}\n'
        report_str += f'\tWith EFO phenotype: {self.with_efo}\n'
        print(report_str)
        return report_str
