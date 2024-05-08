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
        self.exploded_pgx_cat = 0
        self.exploded_drugs = 0
        self.exploded_phenotypes = 0
        # Output counts (after annotation and exploding)
        self.evidence_strings = 0
        self.with_efo = 0
        self.with_consequence = 0
        self.with_target_gene = 0
        self.with_haplotype = 0
        self.resolved_haplotype_id = 0  # indicates we were able to resolve the haplotype to a PGKB internal ID
        # Variant counts
        self.total_rs = 0
        self.rs_with_alleles = 0
        self.rs_with_more_than_2_alleles = 0

    def report(self):
        report_str = f'\nTotal clinical annotations: {self.clinical_annotations}\n'
        report_str += f'\tWith RS: {self.with_rs} ({format_percent(self.with_rs, self.clinical_annotations)})\n'
        report_str += (f'\t\t1. Exploded by allele: {self.exploded_alleles} '
                       f'({format_decimal(self.exploded_alleles, self.with_rs)}x)\n')
        report_str += (f'\t\t2. Exploded by PGx category: {self.exploded_pgx_cat} '
                       f'({format_decimal(self.exploded_pgx_cat, self.exploded_alleles)}x)\n')
        report_str += (f'\t\t3. Exploded by drug: {self.exploded_drugs} '
                       f'({format_decimal(self.exploded_drugs, self.exploded_pgx_cat)}x)\n')
        report_str += (f'\t\t4. Exploded by phenotype: {self.exploded_phenotypes}'
                       f' ({format_decimal(self.exploded_phenotypes, self.exploded_drugs)}x)\n')
        report_str += f'Total evidence strings: {self.evidence_strings}\n'
        report_str += (f'\tWith EFO phenotype: {self.with_efo}'
                       f' ({format_percent(self.with_efo, self.evidence_strings)})\n')
        report_str += (f'\tWith functional consequence: {self.with_consequence} '
                       f'({format_percent(self.with_consequence, self.evidence_strings)})\n')
        report_str += (f'\tWith target gene: {self.with_target_gene} '
                       f'({format_percent(self.with_target_gene, self.evidence_strings)})\n')
        report_str += f'Total RS: {self.total_rs}\n'
        report_str += (f'\tWith parsed alleles: {self.rs_with_alleles} '
                       f'({format_percent(self.rs_with_alleles, self.total_rs)})\n')
        report_str += (f'\t\tWith >2 alleles: {self.rs_with_more_than_2_alleles} '
                       f'({format_percent(self.rs_with_more_than_2_alleles, self.rs_with_alleles)})\n')
        report_str += f'Total haplotype IDs: {self.with_haplotype}\n'
        report_str += (f'\tResolved to PGKB IDs: {self.resolved_haplotype_id} '
                       f'({format_percent(self.resolved_haplotype_id, self.with_haplotype)})\n')
        print(report_str)
        return report_str
