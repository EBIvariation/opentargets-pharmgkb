# opentargets-pharmgkb
Pipeline to provide evidence strings for Open Targets from PharmGKB

Download data:
```
wget https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip
wget https://api.pharmgkb.org/v1/download/file/data/drugs.zip

unzip -j clinicalAnnotations.zip "*.tsv" -d .
unzip -j clinicalAnnotations.zip "CREATED*.txt" -d .
unzip -j drugs.zip "*.tsv" -d .
rm clinicalAnnotations.zip drugs.zip
```