

# Follow instructions in http://patchwork.r-forge.r-project.org to install required packages:
# R packages: DNAcopy, patchworkData, patchwork; Python package: Pysam

# Note: samtools, bcftools and vcfutils.pl must be in PATH variable!!

pileup_ops <- function(reference, bamfile, outfile) {
  cmd <- sprintf("samtools mpileup -f %s %s >%s", reference, bamfile, outfile)
  flag <- system(cmd) # execute
  return(flag == 0) # if successfully, flag is 0
}

vcf_ops <- function(reference, bamfile, outfile) {
  cmd <- sprintf("samtools mpileup -uf %s %s | bcftools view -bvcg - > raw.bcf && bcftools view raw.bcf | vcfutils.pl varFilter -D1000 > %s", reference, bamfile, outfile)
  flag <- system(cmd)
  return(flag == 0)
}

plp_vcf_ops <- function(reference, bamfile) {
  
  output.prefix <- basename(bamfile)
  plp.fl = paste(output.prefix, '.plp', sep="")
  vcf.fl = paste(output.prefix, '.vcf', sep="")
  
  cmd <- sprintf("samtools mpileup -f %s %s >%s", reference, bamfile, plp.fl)
  flag1 <- system(cmd) # execute
  
  cmd <- sprintf("samtools mpileup -uf %s %s | bcftools view -bvcg - > raw.bcf && bcftools view raw.bcf | vcfutils.pl varFilter -D1000 > %s", reference, bamfile, vcf.fl)
  flag2 <- system(cmd)
  
  flag <- all(c(flag1, flag2) == 0) # if successfully, flag is 0
  
  return(flag)
}

library(parallel)
library(patchworkData)
library(patchwork)

tumor.bam.file = 'tumor.bam'
normal.bam.file = 'normal.bam'
reference = 'reference.fasta'

tumor.prefx = basename(tumor.bam.file)
normal.prefix = basename(normal.bam.file)

bamfiles = c(tumor.bam.file, normal.bam.file)

r <- mclapply(bamfiles, function(bamfile) {
  plp_vcf_ops(reference = reference, bamfile = bamfile)
})

if (r != TRUE) {
  warning("Maybe there is error in plp_vcf_ops(...), check by hand to proceed!")
}

tumor.plp.file = paste(tumor.prefx, '.plp', sep="")
normal.plp.file = paste(normal.prefix, '.plp', sep="")

tumor.vcf.file = paste(tumor.prefx, '.vcf', sep="")
normal.vcf.file = paste(normal.prefix, '.vcf', sep="")

# Note: to get patchwork to read compressed vcf and pileup files, replace patchwork internal perl scripts
#+in patchwork subfolder `perl` with mpile2alleles.pl and pile2alleles.pl.

patchwork.plot(
  tumor.bam.file,
  tumor.plp.file,
  tumor.vcf.file,
  normal.bam.file,
  normal.plp.file,
  normal.vcf.file,
  Alpha = 0.0001,
  SD = 1
)


cn2 <- 1
delta <- 0.28
het <- 0.21
hom <- 0.79

myCNFile = list.files(".", pattern = "_copynumbers.Rdata", full.names = TRUE)[1]
patchwork.copynumbers(CNfile=myCNFile, cn2=cn2, delta=delta, het=het, hom=hom)

#save.image()




