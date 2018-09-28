# RecodeVCF
**RecodeVCF.py** is a command line tool for recoding genotype calls 
contained in a VCF file and outputting variants in tabular format with annotation fields expanded.

## Usage

``` sh
cd Pipeline-Tools
python ./RecodeVCF.py --vcf <vcf_file> --output <output_file>
```

Detailed description of additional options available in help menu.

``` sh

python ./RecodeVCF.py --help

usage: RecodeVCF [-h] --vcf VCF_FILE --output OUT_FILE
                 [--info-columns INFO_COLUMNS]
                 [--min-call-depth MIN_CALL_DEPTH]
                 [--missing-data-char MISSING_DATA_CHAR]
                 [--missing-gt-char MISSING_GT_CHAR] [--multiallelic] [-v]

optional arguments:
  -h, --help            show this help message and exit
  --vcf VCF_FILE        Path to vcf file to recode.
  --output OUT_FILE     Path to recoded output file.
  --info-columns INFO_COLUMNS
                        Column-delimited list of INFO columns to include in
                        output. NO SPACES ALLOWED or list will not be parsed!
  --min-call-depth MIN_CALL_DEPTH
                        Minimum read depth required to call variant genotype.
  --missing-data-char MISSING_DATA_CHAR
                        Character used as placeholder for missing VCF info.
  --missing-gt-char MISSING_GT_CHAR
                        Character used as placeholder for missing genotypes.
  --multiallelic        Flag allowing variant records to contain more than one
                        alternate allele. This flag shouldn't really be used.
  -v                    Increase verbosity of the program.Multiple -v's
                        increase the verbosity level: 0 = Errors 1 = Errors +
                        Warnings 2 = Errors + Warnings + Info 3 = Errors +
                        Warnings + Info + Debug
```

## Conceptual overview

[VCF] is the standard format for storing genotype information for a set of individuals.
When doing variant analysis, we're typically interested in two types of information:

1. Sample genotypes at each position
2. Annotations associated with variant alleles

**Problem:** sample genotypes and annotations are not searchable/sortable in VCF format

*Example VCF genotypes*

    Sample_1                    Sample_2
    0/1:1,12:13:11:425,0,11	    ./.:0,0:0:.:0,0,0

*Example VCF annotations*

    Func.refGene=splicing;Gene.refGene=NADK;GeneDetail.refGene=NM_001198995:exon2:c.168-1G>T\x3bNM_001198993

Further complicating the issue, there are several different programs (e.g. Annovar, SnpEff) that each add variant annotations with their own standard format.


**RecodeVCF.py** solves this problem by converting a VCF file into a RecodedVCF format.
**RecodedVCF.py** makes variant information searchable/sortable by:

1. Recoding each genotype into a single number between -1 and 1 
2. Parsing variant annotation data into separate columns

### Recoding algorithm
The goal of variant recoding is to summarize an 
individual's genotype and the confidence associated with that genotype call in a single number.
 
    num_reads = number of reads needed to make a genotype call
    
    For each genotype:
        1. If genotype is ./. (i.e. no call was made), recode as 0
        2. If genotype is homozygous ref:
            If read_depth for sample is > num_reads:
                recoded_gt = -1
            Else:
                recoded_gt = -1 * (reads supporting ref allele / num_reads)
        3. If genotype is alternate allele:
            If read depth for sample > num_reads:
                recoded_gt = 1
            Else:
                recoded_gt = reads supporting alt allele / num_reads
            
Examples:

    let num_reads = 10 for all examples
    
    
    1. GT=Ref, Depth=140 -> Recoded GT = -1
    2. GT=Ref, Depth=7   -> Recoded GT = -0.7
    3. GT=Alt, Depth=40  -> Recoded GT = 1
    4. GT=Alt, Depth=2   -> Recoded GT = 0.2
    5. GT=None, Depth=1  -> Recoded GT = 0

## Input Assumptions

1. VCF format v4.0+
2. Annotated with either [Annovar], [SnpEff], both, or no annotation.
3. No multi-allelic variants
    * Obviously real loci are multi-alleleic but VCF files should be normalized using [BCFtools] Norm to break them into separate lines in VCF
    * Allows each allele to have separate annotations
    
## Output file
RecodeVCF.py outputs a tab-delimited file which has the following properties:
1. First 7 columns are always CHROM, POS, ID, REF, ALT, QUAL, FILTER as they are in VCF without modification
2. Each field in the INFO column becomes a separate column in output (Column headers included)
3. Each sample gets its own column for storing its recoded genotype at each locus
### Advanced parameters

By default the *--min-call-depth* parameter is set to 10. This means that at least 10 reads are required to support an allele in order to make a genotype call. 
This can be adjusted as necessary depending on the specifics of your analysis
    
By default **RecodeVCF** creates columns for ALL info info fields. Specific columns to return can be specified by providing a comma-delimited list using the *--info-columns* parameter.
This could be helpful if you're VCF has lots of annotation fields and your looking to cut down on the output file size.

## Parallelization with CatRecodedVCF.py
The helper program CatRecodedVCF.py is designed to merge RecodedVCFs to facilitate parallelized processing.

### Usage

``` sh

cd Pipeline-Tools
python CatRecodedVCF.py --help

usage: CatRecodeVCF [-h] -i INPUT_FILES [INPUT_FILES ...] --output OUT_FILE
                    [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILES [INPUT_FILES ...]
                        Space-delimited list of RecodedVCF files to combine
  --output OUT_FILE     Path to recoded output file.
  -v                    Increase verbosity of the program.Multiple -v's
                        increase the verbosity level: 0 = Errors 1 = Errors +
                        Warnings 2 = Errors + Warnings + Info 3 = Errors +
                        Warnings + Info + Debug
```

### When is CatRecodedVCF.py useful?
To drastically decrease processing time, VCF files can be split by chromosome using [SnpEff] and recoded in parallel. 
**CatRecodedVCF** is designed to merge these splits back into a single RecodedVCF.

Example:

Recode VCF splits using **RecodeVCF.py**
``` sh
cd ./Pipeline-Tools

python ./RecodeVCF.py -i gt.chr1.vcf -o gt.chr1.rec.vcf
python ./RecodeVCF.py -i gt.chr2.vcf -o gt.chr2.rec.vcf
python ./RecodeVCF.py -i gt.chr3.vcf -o gt.chr3.rec.vcf
```
Merge using **CatRecodedVCF.py**

``` sh
cd ./Pipeline-Tools

python ./CatRecodedVCF.py -i gt.chr1.rec.vcf gt.chr2.rec.vcf gt.chr3.rec.vcf \
    -o genotypes.recoded.vcf
```

[VCF]:http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/
[Annovar]:http://annovar.openbioinformatics.org/en/latest/
[SnpEff]:http://snpeff.sourceforge.net/
[BCFTools]:https://samtools.github.io/bcftools/