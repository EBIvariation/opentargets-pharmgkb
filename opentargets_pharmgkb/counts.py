def format_percent(part, total):
    return f'{part} ({part / total:.2%})'


def format_decimal(part, total):
    return f'{part} ({part / total:.2}x)'


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
        self.with_doe = 0
        self.mean_num_doe = 0
        self.median_num_doe = 0
        self.max_num_doe = 0
        # Variant annotation counts
        self.variant_annotations = 0
        self.variant_anns_with_clinical_anns = 0
        self.matchable_variant_anns = 0
        self.matched_variant_anns = 0
        # Haplotype counts
        self.with_haplotype = 0
        self.resolved_haplotype_id = 0  # indicates we were able to resolve the haplotype to a PGKB internal ID
        # Variant counts
        self.total_rs = 0
        self.rs_with_alleles = 0
        self.rs_with_more_than_2_alleles = 0

    def report(self):
        report_str = f'\nTotal clinical annotations: {self.clinical_annotations}\n'
        report_str += f'\tWith RS: {format_percent(self.with_rs, self.clinical_annotations)}\n'
        report_str += f'\t\t1. Exploded by allele: {format_decimal(self.exploded_alleles, self.with_rs)}\n'
        report_str += (f'\t\t2. Exploded by PGx category: '
                       f'{format_decimal(self.exploded_pgx_cat, self.exploded_alleles)}\n')
        report_str += (f'\t\t3. Exploded by drug: '
                       f'{format_decimal(self.exploded_drugs, self.exploded_pgx_cat)}\n')
        report_str += (f'\t\t4. Exploded by phenotype: '
                       f'{format_decimal(self.exploded_phenotypes, self.exploded_drugs)}\n')
        report_str += f'Total evidence strings: {self.evidence_strings}\n'
        report_str += f'\tWith EFO phenotype: {format_percent(self.with_efo, self.evidence_strings)}\n'
        report_str += (f'\tWith functional consequence: '
                       f'{format_percent(self.with_consequence, self.evidence_strings)}\n')
        report_str += f'\tWith target gene: {format_percent(self.with_target_gene, self.evidence_strings)}\n'
        if self.with_doe > 0:
            report_str += (f'\tWith direction of effect annotation: '
                           f'{format_percent(self.with_doe, self.evidence_strings)}\n')
            report_str += f'\t\tMean per evidence string (when present): {self.mean_num_doe:.2f}\n'
            report_str += f'\t\tMedian per evidence string (when present): {int(self.median_num_doe)}\n'
            report_str += f'\t\tMax per evidence string: {self.max_num_doe}\n'
            report_str += f'Total variant annotations: {self.variant_annotations}\n'
            report_str += (f'\tAssociated to a clinical annotation: '
                           f'{format_percent(self.variant_anns_with_clinical_anns, self.variant_annotations)}\n')
            report_str += (f'\t\tMatchable by genotype: '
                           f'{format_percent(self.matchable_variant_anns, self.variant_anns_with_clinical_anns)}\n')
            report_str += (f'\t\t\tMatched by genotype: '
                           f'{format_percent(self.matched_variant_anns, self.matchable_variant_anns)}\n')
        report_str += f'Total RS: {self.total_rs}\n'
        report_str += f'\tWith parsed alleles: {format_percent(self.rs_with_alleles, self.total_rs)}\n'
        report_str += (f'\t\tWith >2 alleles: '
                       f'{format_percent(self.rs_with_more_than_2_alleles, self.rs_with_alleles)}\n')
        report_str += f'Total haplotype IDs: {self.with_haplotype}\n'
        report_str += f'\tResolved to PGKB IDs: {format_percent(self.resolved_haplotype_id, self.with_haplotype)}\n'
        print(report_str)
        return report_str
