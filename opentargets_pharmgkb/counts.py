def format_percent(part, total):
    return f'{part / total:.2%}'


def format_decimal(part, total):
    return f'{part / total:.2}'


class ClinicalAnnotationCounts:
    """Simple class to hold counts for generating clinical annotation evidence."""

    def __init__(self):
        # Input counts (before annotation and exploding by drug, etc.)
        self.clinical_annotations = 0
        self.with_rs = 0
        # Counts after exploding by each attribute
        self.exploded_alleles = 0
        self.exploded_drugs = 0
        self.exploded_phenotypes = 0
        # Output counts (after annotation and exploding)
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
        # Variant counts
        self.total_rs = 0
        self.rs_with_alleles = 0
        self.rs_with_multiple_alleles = 0

    def report(self):
        report_str = f'\nTotal clinical annotations: {self.clinical_annotations}\n'
        report_str += f'\tWith RS: {self.with_rs} ({format_percent(self.with_rs, self.clinical_annotations)})\n'
        report_str += (f'\t\t1. Exploded by allele: {self.exploded_alleles} '
                       f'({format_decimal(self.exploded_alleles, self.with_rs)}x)\n')
        report_str += (f'\t\t2. Exploded by drug: {self.exploded_drugs} '
                       f'({format_decimal(self.exploded_drugs, self.exploded_alleles)}x)\n')
        report_str += (f'\t\t3. Exploded by phenotype: {self.exploded_phenotypes}'
                       f' ({format_decimal(self.exploded_phenotypes, self.exploded_drugs)}x)\n')
        report_str += f'Total evidence strings: {self.evidence_strings}\n'
        report_str += f'\tWith CHEBI: {self.with_chebi} ({format_percent(self.with_chebi, self.evidence_strings)})\n'
        report_str += (f'\tWith EFO phenotype: {self.with_efo}'
                       f' ({format_percent(self.with_efo, self.evidence_strings)})\n')
        report_str += (f'\tWith functional consequence: {self.with_consequence} '
                       f'({format_percent(self.with_consequence, self.evidence_strings)})\n')
        # report_str += f'\tWith PGKB gene: {self.with_pgkb_gene}\n'
        report_str += (f'\tWith VEP gene: {self.with_vep_gene} '
                       f'({format_percent(self.with_vep_gene, self.evidence_strings)})\n')
        report_str += f'Gene comparisons per annotation\n'
        report_str += (f'\tWith PGKB genes: {self.annot_with_pgkb_genes} '
                       f'({format_percent(self.annot_with_pgkb_genes, self.clinical_annotations)})\n')
        report_str += (f'\tWith VEP genes: {self.annot_with_vep_genes} '
                       f'({format_percent(self.annot_with_vep_genes, self.clinical_annotations)})\n')
        report_str += (f'\tPGKB genes != VEP genes: {self.pgkb_vep_gene_diff} '
                       f'({format_percent(self.pgkb_vep_gene_diff, self.clinical_annotations)})\n')
        report_str += f'Total RS: {self.total_rs}\n'
        report_str += (f'\tWith parsed alleles: {self.rs_with_alleles} '
                       f'({format_percent(self.rs_with_alleles, self.total_rs)})\n')
        report_str += (f'\t\tWith >2 alleles: {self.rs_with_multiple_alleles} '
                       f'({format_percent(self.rs_with_multiple_alleles, self.rs_with_alleles)})\n')
        print(report_str)
        return report_str
