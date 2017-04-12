# Use oncotator to annotate somatic mutations

[Oncotator](http://gatkforums.broadinstitute.org/gatk/discussion/4154/howto-install-and-run-oncotator-for-the-first-time#latest) has been widely used in annotating somatic mutations in TCGA projects. It provides comprehensive annotation result and many downstream bioinformatics tools receive its `maf` output file as standard input such as MutSigCV and MEMo etc.


## Download prerequisite file

1. Download datasource corpus from [GATK forum](http://gatkforums.broadinstitute.org/gatk/discussion/4154/howto-install-and-run-oncotator-for-the-first-time#latest) and uncompress it. The file will consume 40G of your hard disk.
2. Oncotator source code [broad institute github repository](https://github.com/broadinstitute/oncotator/releases) and install it accordingly.
3. (Optional) you may also need to download the other two files from [GATK forum](http://gatkforums.broadinstitute.org/gatk/discussion/4154/howto-install-and-run-oncotator-for-the-first-time#latest) `tx_exact_uniprot_matches.txt` and `tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt` to replace the default version in the datasource corpus.
3. (Optional) Since oncotator is written in python, I may recommend you to install [anaconda python](https://www.continuum.io/downloads), which has many packages pre-installed.


## An example to running oncotator

```bash
#!/bin/bash

infile=icgc_release_18_liver_cancer_and_ACRG_HCC.simple_somatic_mutation.open2.txt
outfile=icgc_release_18_liver_cancer_and_ACRG_HCC.simple_somatic_mutation.open.maf
/ifshk7/BC_PS/yangchao/bin/Oncotator -c tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt --input_format=MAFLITE -v --db-dir oncotator_v1_ds_Jan262014/ $infile $outfile hg19

```

