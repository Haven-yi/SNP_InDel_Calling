# SNP and InDel Calling
This is about detecting Single Nucleotide Polymorphisms (SNPs) and Insertions/Deletions (INDELs) from high-throughput sequencing data. The workflow involves individual-level and population-level variant calling, followed by filtering to retain high-confidence variants. The tools used include GATK, Bcftools, Vcftools, and Bedtools.  

Required Software  
•	GATK (Genome Analysis Toolkit) – for variant calling  
•	Bcftools – for variant annotation and filtering  
•	Vcftools – for additional variant filtering  
•	Bedtools – for variant region manipulation  
Required Files   
•	Reference Genome: reference.fasta (FASTA format)  
•	Indexed BAM Files: Sorted and duplicate-removed BAM files (*.bam and *.bai)
# Step 1: Perform HaplotypeCaller Variant Calling Per Sample
```
gatk --java-options "-Xmx36g -Djava.io.tmpdir=./tmp" HaplotypeCaller \
     --native-pair-hmm-threads 6
     -R reference.fasta \  # Reference genome
     -I $species.rmdup.sort.bam \  # Input BAM file
     -ERC GVCF \  # Generate a GVCF file
     --heterozygosity 0.01 \  # Expected heterozygosity (adjustable based on species)
     -O $species.g.vcf.gz  # Output GVCF file
```
# Step 2: Combine GVCFs from Multiple Samples
1 Combine all individual GVCF files into a single file using CombineGVCFs:
```
gatk --java-options "-Xmx36g -Djava.io.tmpdir=./tmp"  CombineGVCFs \
     -R reference.fasta \  # Reference genome
     --variant sample1.g.vcf.gz \  # Input GVCF file 1
     --variant sample2.g.vcf.gz \  # Input GVCF file 2
     -O combined.g.vcf.gz  # Output combined GVCF
```
Convert the combined GVCF into a final VCF file containing genotype information:
```
gatk --java-options "-Xmx36g -Djava.io.tmpdir=./tmp" GenotypeGVCFs \
     -R reference.fasta \  # Reference genome
     -V combined.g.vcf.gz \  # Input combined GVCF
     -O combined.genotyped.vcf.gz  # Output final VCF file
 ```    
2 Combine all individual GVCF files into a single file using GenomicsDBImport
```
gatk --java-options "-Xmx36g -Djava.io.tmpdir=./tmp -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenomicsDBImport \
    --sample-name-map gvcf.map \
    --genomicsdb-workspace-path genomeDB.chr \
    --reader-threads 1 \
    --batch-size 50 \
    --tmp-dir ./tmp 
```
Perform GenotypeGVCFs using GenomicsDB
```
gatk --java-options "-Xmx36g -Djava.io.tmpdir=./tmp" \
    GenotypeGVCFs \
    -R reference.fasta \
    -V gendb://genomeDB.chr \
    -O genome.combined.genotyped.vcf.gz
```
# Step 3: Variant Filtering
To ensure each variant has a unique identifier
```
bcftools annotate --set-id +'%CHROM\:%POS' -O z -o variant.annote.vcf.gz variant.vcf.gz
```
Remove SNPs Near INDELs (5bp Window)
```
#Separate SNPs and INDELs
vcftools --gzvcf variant.annote.vcf.gz --remove-indels --maf 0.0000001 --recode --recode-INFO-all --out variant.annote.snp
bgzip variant.annote.snp.recode.vcf

vcftools --gzvcf variant.annote.vcf.gz --keep-only-indels --maf 0.0000001 --recode --recode-INFO-all --out variant.annote.indel
bgzip variant.annote.indel.recode.vcf
```
Filter SNPs Near INDELs
```
bcftools view -h variant.annote.vcf.gz > header.txt
bedtools window -a variant.annote.snp.recode.vcf -b variant.annote.indel.recode.vcf -w 5 > variant.rm_indel_mark.vcf
cat header.txt variant.rm_indel_mark.vcf | bgzip > variant.rm_indel_mark.vcf.gz
```
Retain Only Biallelic SNPs
```
vcftools --gzvcf variant.rm_indel_mark.vcf.gz --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all --out variant.biallelic
bgzip variant.biallelic.recode.vcf
```
Filter by Read Depth (DP > 5)
```
#remove low-depth variants
vcftools --gzvcf variant.biallelic.recode.vcf.gz --minDP 5 --recode --recode-INFO-all --out variant.DP5
bgzip variant.DP5.recode.vcf
```
Filter by Genotyping Quality (GQ > 10)
```
#retain high-confidence genotypes
vcftools --gzvcf variant.DP5.recode.vcf.gz --minGQ 10 --recode --recode-INFO-all --out variant.GQ10
bgzip variant.GQ10.recode.vcf
```
Filter by Missing Data and Minor Allele Frequency
```
vcftools --gzvcf variant.GQ10.recode.vcf.gz --max-missing 0.8 --maf 0.01 --recode --recode-INFO-all --out variant.filtered
bgzip variant.filtered.recode.vcf
```


