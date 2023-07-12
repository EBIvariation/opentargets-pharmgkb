# opentargets-pharmgkb
Pipeline to provide evidence strings for Open Targets from PharmGKB

```
# Download data
DATA_DIR=<directory for data>
wget https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip
wget https://api.pharmgkb.org/v1/download/file/data/drugs.zip

unzip -j clinicalAnnotations.zip "*.tsv" -d $DATA_DIR
unzip -j clinicalAnnotations.zip "CREATED*.txt" -d $DATA_DIR
unzip -j drugs.zip "*.tsv" -d $DATA_DIR
rm clinicalAnnotations.zip drugs.zip

# Run pipeline
generate_evidence.py --data-dir $DATA_DIR --created-date <created date> --output-path evidence.json
```