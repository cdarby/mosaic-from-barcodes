# mosaic-from-barcodes
Detect mosaic variants using 10x Genomics linked-read barcode information.

### Overview

Mosaic variants can be confidently identified from sequence data using sequence read phasing (see Figure 1 of [this article][plos-genetics]). `mosaic-from-barcodes` contains software to identify mosaic variants from [VCF files][vcf-format] with barcode information incorporated in `BX` tags, such as those produced by [Long Ranger][long-ranger].

The identification of mosaic variants occurs according to the following steps:
* Variants are read from the input file and stored in a list. The barcodes of possible mosaic and germline variants are stored in a Python dictionary (hash table).
* Possible mosaic variants that occur some distance (default is 500,000bp) upstream of the current variant are phased to nearby germline variants.
  - A sparse matrix of allele information from nearby germline variants (rows) and barcodes (columns) is constructed.
  - The haplotype of each barcode is determined by majority vote. Low-quality barcodes are filtered.
  - The candidate mosaic variant is phased to the identified haplotypes using the barcode information.
  - If phasing indicates high-confidence mosaic status of the variant, a MOSAIC tag is added to the variant's INFO field.
* Variants twice as far upstream as the phasing position from the current variant are output and their barcode information is removed from the dictionary.

### Usage

Clone the repository.
```
git clone git@github.com:DonFreed/mosaic-from-barcodes.git
```

Make a directory for storing test data and download input and annotation data.
```
mkdir mosaic-from-barcodes/test
cd mosaic-from-barcodes/test
# Download the input variants #
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/genome/NA12878_WGS_210/NA12878_WGS_210_phased_variants.vcf.gz
# Download a dbSNP and 1000 Genomes VCF #
wget --user=gsapubftp-anonymous \
ftp://ftp.broadinstitute.org/bundle/2.8/b37/dbsnp_138.b37.vcf.gz \
ftp://ftp.broadinstitute.org/bundle/2.8/b37/1000G_phase3_v4_20130502.sites.vcf.gz
cd ..
```

Annotate the 10x VCF with information on variants in dbSNP and 1000 Genomes.
```
(zcat test/NA12878_WGS_210_phased_variants.vcf.gz | head -n 91; \
zcat  test/NA12878_WGS_210_phased_variants.vcf.gz | tail -n +92 | sort -k1,1V -k2,2n) | \
python3 annotate_variants.py --ref_var_name ONEKG DBSNP \
--ref_var_file test/1000G_phase3_v4_20130502.sites.vcf \
test/dbsnp_138.b37.vcf > test/NA12878_10x_annotated.vcf
```

Identify mosaic variants from the annotated VCF.
```
python3 main.py --infile test/NA12878_10x_annotated.vcf --marks_germline ONEKG \
--absent_marks_possible_mosaic DBSNP --outfile test/NA12878_mosaic5.vcf \
--ignore_filter 10X_PHASING_INCONSISTENT > test/NA12878_10x_mosaic.vcf
```

[plos-genetics]: http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006245
[vcf-format]: https://samtools.github.io/hts-specs/VCFv4.2.pdf
[long-ranger]: http://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger