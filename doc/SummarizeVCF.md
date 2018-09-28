# SummarizeVCF
**SummarizeVCF.py** is a command line tool for producing per-sample summary statistics for variants in a VCF file.
Conceptually, **SummarizeVCF.py** is intended to be for VCF files what [FastQC] is for FASTQ files.

## Usage

``` sh
cd Pipeline-Tools
python ./SummarizeVCF.py <summary_type> --vcf <vcf_file> [options]

summary types:
 Mutect, Multisample

```

Detailed description of additional options available in help menu.

``` sh

python ./SummarizeVCF.py --help


usage: SummarizeVCF.py <summary_type> [options]

positional arguments:
  {Mutect,Multisample}  VCF Summary type.

optional arguments:
  -h, --help            show this help message and exit
  --vcf VCF_FILE        Path to vcf file to summarize.
  --max-records MAX_RECORDS
                        Maximum number of records to process. Default: ALL.
  --max-indel-len MAX_INDEL_LEN
                        Upper bound of indel length summary.
  --max-depth MAX_DEPTH
                        Upper bound of variant depth summary.
  --max-qual MAX_QUAL   Upper bound of variant quality summary.
  --afs-bins NUM_AFS_BINS
                        Number of bins to use for Allele Frequency Spectrum.
  -v                    Increase verbosity of the program.Multiple -v's
                        increase the verbosity level: 0 = Errors 1 = Errors +
                        Warnings 2 = Errors + Warnings + Info 3 = Errors +
                        Warnings + Info + Debug
```

## Why are there multiple summary modes?

Although [VCF] is supposed to be a standardized format, variant calling programs differ in how they utilize specific fields.
Currently, **VCFSummary.py** is able to handle both standard multisample VCF files (produced by GATK [GenotypeGVCFs] or Samtools [mpileup]):

``` sh
python ./SummarizeVCF.py Multisample --vcf <vcf_file> [options]
```

And VCF files produced by somatic callers like [Mutect2].

``` sh
python ./SummarizeVCF.py Mutect --vcf <vcf_file> [options]
```
**Disclaimer**: Choosing the correct summary type is important as a mismatch will cause the summary program to crash.
Make sure you know what kind of VCF you're working with. 

## What kind of summary statistics does VCFSummary.py compute?
For each sample in VCF, **SummarizeVCF.py** reports the number of:
* Missing GT 
* Called GT 
* Variant GT 
* Heterozygous loci
* Homozygous-Alt 
* Deletions 
* Insertions
* Monomorphs 
* SNPs 
* Ts 
* Tv
* 12 SNP transition types (e.g. C -> A, G -> T)
* dbSNP variants
* Structural variants

If annotation information is provided, the following can also reported if present:
* Number of variants by impact (HIGH, MODERATE, LOW, MODIFIER)
* Number of variants by type
    * intergenic
    * intronic 
    * synonymous_SNV 
    * nonsynonymous_SNV 
    * stoploss,stopgain,
    * onframeshift_deletion
    * nonframeshift_insertion
    * frameshift_deletion
    * frameshift_insertion
    
 **SummarizeVCF.py** also produces per-sample distribution summaries of the following:
* Genotype Quality
* Read depth
* Allele frequency (AF=percentage of samples with that variant)
* Insertion length
* Deletion length

If run in *Mutect* mode, the following are also computed for each sample:
* Number variants passed Mutect filter


## Input Assumptions

1. VCF format v4.0+
2. Annotated with either [Annovar], [SnpEff], both, or no annotation.
3. No multi-allelic variants
    * Obviously real loci are multi-alleleic but VCF files should be normalized using [BCFtools] Norm to break them into separate lines in VCF
    * Allows each allele to have separate annotations
4. If in Mutect mode, VCF file **must** be Mutect output or merged VCF from multiple Mutect runs
5. If in Multisample mode, VCF **must not** be Mutect output

## Additional options
*--max-records* option can be used to subsample the number of VCF records processed for faster runtimes. Default is to process all records.

*--max-indel-len* sets the upper bound for summarizing the indel size distribution

*--max-depth* sets the upper bound for summarizing the variant read depth distribution

*--afs-bins* specifies the number of bins for summarizing the allele frequency spectrum of alternate alleles

## Output format

**SummarizeVCF.py** produces a tab-delimited output file with 6 sections:
1. COUNTS: counts described above for each sample.
2. DEPTH: variant depth distribution for each sample
3. QUAL:  variant quality score distribution for each sample
4. INSERT_LEN: insertion length distribution for each sample
5. DELETE_LEN: deletion length distribution for each sample
6. AAFS: Allele frequency spectrum for each sample

All sections have a first header row followed by subsequent rows for each sample. 

Check out the [Example VCFSummary](./vcf_summary_example.txt) for a better idea.

## Parallelization with CatVCFSummary.py
The helper program CatVCFSummary.py is designed to merge VCFSummaries to facilitate parallelized processing.

### Usage

``` sh

cd Pipeline-Tools
python CatVCFSummary.py --help

usage: CatVCFSummary [-h] -i INPUT_FILES [INPUT_FILES ...] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILES [INPUT_FILES ...]
                        Space-delimited list of VCFSummary files to combine
  -v                    Increase verbosity of the program.Multiple -v's
                        increase the verbosity level: 0 = Errors 1 = Errors +
                        Warnings 2 = Errors + Warnings + Info 3 = Errors +
                        Warnings + Info + Debug
```

### When is CatVCFSummary.py useful?
To drastically decrease processing time, VCF files can be split by chromosome using [SnpEff] and summarized in parallel. 
**CatVCFSummary.py** is designed to merge these splits back into a single VCFSummary.

Example:

Summarize variants in VCF splits using **SummarizeVCF.py**
``` sh
cd ./Pipeline-Tools

python ./SummarizeVCF.py Multisample -vcf chr1.vcf > chr1.sum.txt
python ./SummarizeVCF.py Multisample -vcf chr2.vcf > chr2.sum.txt
python ./SummarizeVCF.py Multisample -vcf chr3.vcf > chr3.sum.txt
```
Merge using **CatVCFSummary.py**

``` sh
cd ./Pipeline-Tools

python ./CatVCFSummary.py -i chr1.sum.txt chr2.sum.txt chr3.sum.txt \
    -o full.sum.txt
```

[VCF]:http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/
[Annovar]:http://annovar.openbioinformatics.org/en/latest/
[SnpEff]:http://snpeff.sourceforge.net/
[BCFTools]:https://samtools.github.io/bcftools/
[FastQC]:https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[Mutect2]:https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php
[GenotypeGVCFs]:https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php
[mpileup]:http://samtools.sourceforge.net/mpileup.shtml

