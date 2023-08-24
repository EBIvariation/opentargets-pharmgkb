# opentargets-pharmgkb
Pipeline to provide evidence strings for Open Targets from PharmGKB

## How to run
```
# Download data
DATA_DIR=<directory for data>
wget https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip
wget https://api.pharmgkb.org/v1/download/file/data/drugs.zip
wget https://api.pharmgkb.org/v1/download/file/data/variants.zip

unzip -j clinicalAnnotations.zip "*.tsv" -d $DATA_DIR
unzip -j clinicalAnnotations.zip "CREATED*.txt" -d $DATA_DIR
unzip -j drugs.zip "*.tsv" -d $DATA_DIR
unzip -j variants.zip "*.tsv" -d $DATA_DIR
rm clinicalAnnotations.zip drugs.zip variants.zip

# Run pipeline
generate_evidence.py --data-dir $DATA_DIR --fasta <path to fasta> --created-date <created date> --output-path evidence.json
```

## Schema documentation

Unless otherwise mentioned, data is taken directly from PharmGKB.

Field | Description | Example
--|--|--
datasourceId | Identifier for data source | `"pharmGKB"`
datasourceVersion | Date when data dump was generated, formatted YYYY-MM-DD | `"2023-03-23"`
datatypeId | Type of data corresponding to this evidence string (currently only clinical annotation) | `"clinical_annotation"`
studyId | Clinical Annotation ID | `"1449309937"`
evidenceLevel |  Level of evidence (see [here](https://www.pharmgkb.org/page/clinAnnLevels)) | `"1A"`
literature | List of PMIDs associated with this clinical annotation | `["7857962", "1389482"]`
variantId | VCF-style (`chr_pos_ref_alt`) identifier of variant; computed as described [below](#variant-coordinate-computation) | `"21_45514946_CAAG_C"`
variantRsId | RS ID of variant | `"rs121918596"`
variantFunctionalConsequenceId | Sequence Ontology term, currently from VEP only | `"SO_0001624"`
targetFromSourceId | Ensembl stable gene ID, currently from VEP only | `"ENSG00000173638"`
genotype | Genotype string | SNP `"TA"`, indel `"del/AAG"`, repeat `"(CA)16/(CA)17"`
genotypeAnnotationText | Full annotation string for genotype | `"Patients with the rs121918596 del/AAG genotype may develop malignant hyperthermia when treated with volatile anesthetics [...]"`
drugFromSource | Drug name | `"isoflurane"`
drugId | CHEBI ID of drug, mapped through OLS | `"CHEBI_6015"`
pgxCategory | Pharmacogenomics phenotype category | `"Toxicity"`
phenotypeText | Phenotype name | `"Malignant Hyperthermia"`
phenotypeFromSourceId | EFO ID of phenotype, mapped through ZOOMA / OXO | `"Orphanet_423"`

### Variant coordinate computation

TODO: describe this in words

````mermaid
graph TD
    J[PharmGKB]
    H[FASTA files]
    E[Clinical alleles table]
    A[Variant table]    
    D[Generate 'chr_pos_ref_alt' identifier]
    S[NCBI Genome Assembly]
    J --> A
    J --> E
    S --> H
    A --> |locations 'Chr+Pos'| D
    H --> |Reference + context| D
    E --> |Alternate alleles| D
````
