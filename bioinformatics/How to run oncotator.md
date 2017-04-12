# Use oncotator to annotate somatic mutations

[Oncotator](http://gatkforums.broadinstitute.org/gatk/discussion/4154/howto-install-and-run-oncotator-for-the-first-time#latest) has been widely used in annotating somatic mutations in TCGA projects. It provides comprehensive annotation result and many downstream bioinformatics tools receive its `maf` output file as standard input such as MutSigCV and MEMo etc.


## Download prerequisite file

1. Download datasource corpus from [GATK forum](http://gatkforums.broadinstitute.org/gatk/discussion/4154/howto-install-and-run-oncotator-for-the-first-time#latest) and uncompress it. The file will consume 40G of your hard disk.

```bash
# download from Broad Institute FTP
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/oncotator/oncotator_v1_ds_April052016.tar.gz

# uncompress database
tar zxvf oncotator_v1_ds_April052016.tar.gz
```

2. Oncotator source code [broad institute github repository](https://github.com/broadinstitute/oncotator/releases) and install it accordingly.

```bash
# Download oncotator source code
wget https://github.com/broadinstitute/oncotator/archive/v1.9.2.0.tar.gz 

# Install
pip install v1.9.2.0.tar.gz
```

3. (Optional) you may also need to download the other two files from [GATK forum](http://gatkforums.broadinstitute.org/gatk/discussion/4154/howto-install-and-run-oncotator-for-the-first-time#latest) `tx_exact_uniprot_matches.txt` and `tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt` to replace the default version in the datasource corpus.

4. (Optional) Since oncotator is written in python, I may recommend you to install [anaconda python](https://www.continuum.io/downloads), which has many packages pre-installed.


## An example to running oncotator

```bash
#!/bin/bash
infile=../release_23/icgc.pancreatic_cancer.txt 
outfile=icgc.pancreatic_cancer.maf
Oncotator -c /expt/lixc/software/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt --input_format=MAFLITE -v --db-dir /expt/lixc/software/oncotator_v1_ds_April052016 $infile $outfile hg19
```

