# ConGen 2025 - Deleterious variants tutorial

## Introduction

[SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/) and [SnpSift](https://pcingola.github.io/SnpEff/snpsift/introduction/) are powerful tools designed for annotating and analyzing genetic variants, with a focus on deleterious variants. This tutorial will guide you through the process of using these tools to work with non-model organisms. 
SnpEff is a tool for annotating and analyzing genetic variants, with a focus on deleterious variants. There are prebuild databases for thousands of species, but if you're working with non-model organisms, you might have to create a custom SnpEff database using the organism's reference genome and annotation files. The annotation files usually include gene models in GFF3, GTF, or Gencode format. Another option is to map your reads to a closely related species that already has its genome in the SnpEff database.

### 1) Mapping the genome against a reference with SnpEff database

To map the genome of your non-model organism against a reference with a SnpEff database, you will first need to align the reads to the reference genome using an aligner like BWA or Bowtie2. Once you have the alignments in BAM or SAM format, you can call variants using a tool like GATK, FreeBayes, or Samtools.

a. Prepare the *vcf* file
If we want we could just use the vcf file as it is, but we can also do some more filtering. Here I decided to remove indels, and only keep bi-allelic sites (SNPs with two alleles, not more) with no missing data:

```javascript
vcftools --vcf /data/genomics/workshops/smsc_2024/Deleterious/GuamRails_ptg000009l.vcf --remove-indels  \
 --max-missing 1.0 --max-alleles 2 --recode --out GuamRails_ptg000009l_SNPs
```

### 2) Preparing the annotation database - Creating a custom SnpEff database

To create a custom SnpEff database, follow these steps:

a. Download the SnpEff software and set up the environment:

```javascript
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff
```

b. Create a new directory for your non-model organism in the data folder:

```javascript
mkdir data/your_organism
```

c. Copy the reference genome (FASTA) and annotation (GFF3, GTF, or Gencode) files into the new directory:

```javascript
cp /path/to/reference_genome.fasta data/your_organism/
cp /path/to/annotation_file.gff3 data/your_organism/
```

d. Modify the snpEff.config file to include the new genome:

```javascript
echo "your_organism.genome : Your_Organism" >> snpEff.config
```

e. Build the SnpEff database:

```javascript
java -jar snpEff.jar build -gff3 -v your_organism
```

### 3) Run the SnpEff to annotating variants 

With the custom SnpEff database created and the variants called, you can now annotate the variants using SnpEff:

```javascript
java -jar snpEff.jar Elephant input.vcf -stats ${OUT}/${prefix}_summary.html -csvStats ${OUT}/${prefix}_summary.csv > annotated_output.vcf
```

Familiarize yourself with the outputs. SnpEff will produce:
1) An annotated vcf file. Check the info field (column X do the vcf) containing the annotation.

![Screenshot 2024-12-16 160012](https://hackmd.io/_uploads/BJyjgp641e.png)

2) a text file summarizing the number of variant types per gene
3) an HTML file containing summary statistics about the variants and their annotations. Copy it to your local computer and open in a web browser.
 
```
# Replace YOUR_USER and YOUR_FOLDER to your own username and folder name
# Open a terminal on your ON computer (not connected to hydra) and type: 

scp YOUR_USER@160.111.215.42://scratch/genomics/YOUR_FOLDER/Deleterious/vep_rail.txt_summary.html .
```


### 4) Analyzing deleterious variants with SnpSift

SnpSift is a collection of tools that can be used to filter, sort, and analyze the annotated VCF files. 


#### 4.1) Impact factors

Let's start filtering the VCF annotation file by impact.
SnpEff classify the following impact factors for all our variants:

Impact factor|	Description
---|---
LOW | synonymous variants
MODERATE | non-synonymous variants
HIGH | non-sense variants (affect start, stop or reading frame)
MODIFIER | all other variants (for example intronic, UTRs, intergenic)

A more detailed description can be found [here](https://pcingola.github.io/SnpEff/snpeff/inputoutput)
We will focus on the classes 'MODERATE' and 'HIGH', assuming that they represent potentially deleterious variants (This will not always be true, but we will never know the exact effect of all mutations, not even in model organisms, and that's just something we have to live with!). However, to have something potentially neutral to compare with, we will keep the LOW category.

Create subsets of annotated SNP with LOW, HIGH and MODERATE impact variants:

```javascript
java -jar SnpSift.jar filter "ANN[0].IMPACT has 'LOW'" input.vcf > low_impact_output.vcf

java -jar SnpSift.jar filter "ANN[0].IMPACT has 'HIGH'" input.vcf > high_impact_output.vcf

java -jar SnpSift.jar filter "ANN[0].IMPACT has 'MODERATE'" input.vcf > moderate_impact_output.vcf

```

#### 4.2) Masked and realized load

From the lectures, we recall that masked load comes from heterozygous deleterious mutations, and realized load comes from homozygous deleterious mutations. First we can just look at different genotype counts from one individual:

```javascript
java -jar SnpSift.jar filter "ANN[0].IMPACT has 'LOW'" input.vcf > low_impact_output.vcf

java -jar SnpSift.jar filter "ANN[0].IMPACT has 'HIGH'" input.vcf > high_impact_output.vcf

java -jar SnpSift.jar filter "ANN[0].IMPACT has 'MODERATE'" input.vcf > moderate_impact_output.vcf
```

Do you see a difference between the individuals? Or a difference between 'HIGH' impact mutations and 'MODERATE' impact mutations?

Remenber, this is just a very rough estimation of genetic load! Apart from finding the correct deleterious allele (which we ignored above), it might be necessary to account for differences in sequencing (if some individuals have more missing data, for example). One way to do this is to calculate load as the number of deleterious alleles per genotyped sites.


#### 4.3) Analyze the alleles
You may wonder: which of the two alleles in a site is the deleterious one?? SnpEff doesn't provide us with this information. We know that a mutation in a certain position causes a non-synonymous variation, and that this could be harmful. But for all we know, it could be our reference individual who is carrying the deleterious allele, and the other individuals carrying the 'normal' variant.

There are a few different ways to figure this out (see box below), but for now we will assume that the REFERENCE allele is 'normal', and that the ALTERNATIVE allele is deleterious.

:::success
**Finding out which allele is deleterious**
In a large population, deleterious variants will most likely be selected against, and will never reach high frequencies. Therefore it is often safe to assume that the minor allele is the deleterious variant. But in a small population, we know that also deleterious variants can reach high frequencies just due to drift. 
And what if we only have a couple of samples in the first place? It is hard to tell which is the minor! 
Another option is to polarize the variants, which involves determining the ancestral allele and designating it as the 'normal' variant. This strategy has been widely used in conservation genetics. Polarization is particularly valuable when working with multiple populations/species and using a reference genome that is more closely related to one of the populations. In such cases, the population chosen as the reference may appear to have less variations simply because it is being treated as the reference. Thus, polarization helps mitigate any bias introduced by the reference genome.
:::


A good method to see if our potentially deleterious sites are under more selective constraints than for example synonymous mutations, is to compare their allele frequency spectra.

a) First we'll just look at the genotypes (remove everything else)

```javascript
grep -v "##" GuamRail_moderate_impact.recode.vcf | awk '{out=$1"\t"$2; for(i=10; i<=11; i++){split($i,s,":"); out=out"\t"s[1]}; print out}' |less
```

This is a good start to just get a feeling for the data.
Now we will use vcftools to calculate allele frequency for each site. This will create .frq output files that we can summarize into a little table.

b) Create a table header:

```javascript
echo "Type Number Ref_freq Alt_freq" |sed 's/ /\t/g' >SFS.txt
```

c) Loop over the types, create a frequency table and summarize:

```javascript
for type in "low" "moderate" "high"
do
  vcftools --vcf GuamRail_${type}_impact.recode.vcf --freq --out GuamRail_${type}_impact
  tail -n+2 GuamRail_${type}_impact.frq  |cut -f5,6 |sed 's/:/\t/g' | cut -f2,4 |sort |uniq -c |awk -v t=$type -v OFS="\t" '{print t,$0}' >>SFS.txt
done
```

The file SFS.txt contains the site frequency spectra for all the three types of mutations. 

d) Below is some R code you can use for plotting:


```r
require(tidyverse)

file<-"SFS.txt"
SFS<-file %>% read.table(header=TRUE) %>% as_tibble()

# Order the types
SFS$Type<-factor(SFS$Type, levels=c("high","moderate","low"))

ggplot(SFS, aes(x=Alt_freq, y=Number, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge())
```

What do we see here? The sizes are so different between the types so they are hard to compare! We can try making the bars using relative sizes instead:

```r
#With relative sizes
SFS_rel <- SFS %>% group_by(Type) %>% mutate(Rel=Number/sum(Number))

ggplot(SFS_rel, aes(x=Alt_freq, y=Rel, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge())
```

Now it looks better! We see that the 'HIGH' category is shifted to the left, and the 'LOW' category has relatively more sites with higher alternative frequency. Can you explain why?



