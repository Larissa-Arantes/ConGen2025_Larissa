# ConGen 2025 - Genetic diversity tutorial

## Introduction

Heterozygosity is a simple and informative statistic that can be obtained by analyzing whole-genome data. You can calculate the average heterozygosity of an individual or assess the local heterozygosity throughout the genome to answer more specific questions.

There are many methods available to achieve this. In this tutorial, we will use two strategies. To obtain the average heterozygosity, we will manipulate the VCF file using bcftools. To explore heterozygosity variation across the genome, we will use ANGSD and other associated tools.

### Average heterozygosity

#### a) Use the following script to calculate the average heterozygosity per individual:

```bash
#!/bin/bash
# FILENAME:  avg_het
#SBATCH -A bio240351  # Allocation name
#SBATCH --nodes=1         # Total # of nodes (must be 1 for serial job)
#SBATCH --ntasks=1        # Total # of MPI tasks (should be 1 for serial job)
#SBATCH --time=1:30:00    # Total run time limit (hh:mm:ss)
#SBATCH -J avg_het     # Job name
#SBATCH -o /home/YOUR_USERNAME/logs/avg_het.o%j      # Name of stdout output file
#SBATCH -e /home/YOUR_USERNAME/logs/avg_het.e%j      # Name of stderr error file
#SBATCH -p wholenode  # Queue (partition) name
#SBATCH --mail-user=x-YOUR_USERNAME@anvil.rcac.purdue.edu
#SBATCH --mail-type=all   # Send email to above address at begin and end of job

#Load bcftools
module load biocontainers/default
module load bcftools/1.17

#Define path to input files - do not change it
VCF_FILE="/home/x-larantes/00_raw_data/elephants.chr.11.12.subset.vcf.gz"
POPFILE="/home/x-larantes/x_scripts/popfile.txt"
GENOME_LENGTH="152468515" #Sum of the chromosome lengths included in this analysis

#Make a directory for saving outputs. Replace it with your actual username
OUTPUT_PATH="/home/x-YOUR_USERNAME/05_heterozygosity"
mkdir $OUTPUT_PATH

#Define the output name. Replace it with your actual username
OUTPUT_FILE="/home/x-YOUR_USERNAME/05_heterozygosity/elephants_heterozygosity.tsv"

#Get a list of sample names from the VCF file
SAMPLES=$(bcftools query -l $VCF_FILE)

#Write a header line to the output file
echo -e "Sample\tHeterozygous_sites\tHeterozygosity" > $OUTPUT_FILE

#Loop through each sample and calculate the heterozygosity
for SAMPLE in $SAMPLES; do
  HETEROZYGOUS=$(bcftools view -s $SAMPLE $VCF_FILE | grep -v "#" | grep -o "0/1" | wc -l)
  HETEROZYGOSITY=$(echo "scale=7; $HETEROZYGOUS / $GENOME_LENGTH" | bc)
  echo -e "$SAMPLE\t$HETEROZYGOUS\t$HETEROZYGOSITY" >> $OUTPUT_FILE
done

#Add population origin of each sample to the heterozygosity file
awk 'NR==FNR {pop[$1]=substr($0, index($0, $2)); next} FNR==1 {print $0"\tPopulation"} FNR>1 {print $0"\t"pop[$1]}' $POPFILE $OUTPUT_FILE > temp && mv temp $OUTPUT_FILE
```

This is a bash script that calculates the heterozygosity of each sample in a VCF file. It uses the bcftools command line tool to extract the list of sample names and calculate the number of heterozygous variants for each sample. The results are written to a tab-separated file with the name "heterozygosity.tsv" and two columns: "Sample" and "Heterozygosity". Then, the last command line add a last column with population origin of each individual. 
Copy this script to a file, give a name that ends with `.sh`, and run it on Anvil using `sbatch command.sh`.

#### b) Plot the results using the R

You can plot it with R in your local computer or using the interactive mode in the server.
If you decide to use your own computer, copy the file `elephants_heterozygosity.tsv` to your local computer using `scp`.
If you decide to use Anvil interact mode (login with `scp -X`), load the R module using `module load r/4.1.0` and then open R typing `R` in the terminal.

```r
# Call ggplot library
library(ggplot2)

# Choose the working directory
setwd("/home/x-YOUR_USERNAME/heterozygosity/")

# Define the input file
INPUT_FILE <- "/home/x-YOUR_USERNAME/heterozygosity/elephants_heterozygosity.tsv"

# Read in the data from the input file
data <- read.table(INPUT_FILE, header=TRUE, sep="\t")

# Create a boxplot for heterozygosity by population
plot <- ggplot(data, aes(x=Population, y=Heterozygosity, color=Population)) +
  geom_boxplot(outlier.shape = NA, width=0.6, alpha=0.6, color="black") +  # Boxplot without outliers, adjust box width and transparency
  geom_jitter(size=2, width=0.1, alpha=0.7) +  # Jittered dots for each individual with some horizontal spacing
  theme(axis.text.x=element_text(angle=45, hjust=1)) +  # Angle the x-axis labels for better readability
  labs(x="Population", y="Heterozygosity", color="Population") +
  scale_color_discrete(name="Population")  # Color legend for populations

# Save the plot to a file
ggsave("heterozygosity_boxplot_by_population.png", plot, width=10, height=6, dpi=300)

# Show the plot
print(plot)

# Open with interactive mode
svg("heterozygosity_boxplot_by_population.svg", width=10, height=6)
print(plot)
dev.off()
```

This R code block generates a boxplot of heterozygosity values, grouping individuals by population, using data from a tab-separated value (TSV) file. The ggplot2 library is loaded, and the name of the input file is specified. The data is read in using the read.table function, and a boxplot is created using ggplot with "Populations" on the x-axis and "Heterozygosity" on the y-axis. The resulting plot is saved as a png file named "heterozygosity_boxplot_by_population.png".

> [!IMPORTANT]
> :elephant::grey_question: How does the genetic diversity of savanna elephants compare to that of forest elephants, and what factors might contribute to these differences? What could explain the observed levels of heterozygosity in the savanna elephants of Queen Elizabeth National Park (QENP), Uganda?
  
### Genome-wide heterozygosity

Another way to display heterozygosity is by using a window-based approach. This approach allows for the visualization of variation throughout the genome and the ability to focus on regions of particularly high or low diversity. It also permits comparison of results with other methods, such as Runs of Homozygosity (RoH). There are several software programs that can perform this task in a similar way, including ANGSD.

ANGSD (Analysis of Next Generation Sequencing Data) is an open-source software designed to analyze low-depth Next Generation Sequencing (NGS) data. Developed by researchers at the University of Copenhagen, ANGSD was specifically designed to handle large-scale sequencing data, particularly from non-model organisms. ANGSD provides a flexible and efficient framework for analyzing complex genomic data, enabling researchers to perform various types of population genetics and genomics analyses.

Key capabilities of ANGSD include:

1.  Genotype calling: ANGSD can estimate genotype likelihoods from sequencing data without actually calling genotypes, which helps reduce biases and errors introduced by hard genotype calls.
2.  SNP discovery: The software can identify Single Nucleotide Polymorphisms (SNPs) and estimate their allele frequencies while accounting for sequencing errors and varying levels of sequencing depth.
3.  Genetic association studies: ANGSD can perform Genome-Wide Association Studies (GWAS) and estimate genotype-phenotype associations using mixed linear models or generalized linear models.
4.  Population structure: ANGSD can estimate population structure and admixture proportions using principal component analysis (PCA) or model-based approaches.
5.  Selection scans: The software can detect signals of positive or balancing selection, using various statistics like FST, Tajima's D, and nucleotide diversity.
6.  Handling of various file formats: ANGSD can work with various input file formats, including BAM, CRAM, and VCF files, and output results in multiple formats suitable for downstream analyses.

Overall, ANGSD is a versatile and powerful tool for analyzing NGS data, particularly for non-model organisms or low-coverage sequencing projects. It provides a comprehensive suite of analysis tools that cater to various research objectives in population genetics and genomics.

#### a) Run heterozygosity and realSFS with ANGSD using the bash script below:

```bash
#!/bin/bash
# FILENAME:  genome_wide_het
#SBATCH -A bio240351  # Allocation name
#SBATCH --nodes=1         # Total # of nodes (must be 1 for serial job)
#SBATCH --ntasks=1        # Total # of MPI tasks (should be 1 for serial job)
#SBATCH --time=1:30:00    # Total run time limit (hh:mm:ss)
#SBATCH -J genome_wide_het     # Job name
#SBATCH -o /home/x-YOUR_USERNAME/logs/genome_wide_het.o%j      # Name of stdout output file
#SBATCH -e /home/x-YOUR_USERNAME/logs/genome_wide_het.e%j      # Name of stderr error file
#SBATCH -p wholenode  # Queue (partition) name
#SBATCH --mail-user=x-YOUR_USERNAME@anvil.rcac.purdue.edu
#SBATCH --mail-type=all   # Send email to above address at begin and end of job

#Load angsd
module load biocontainers/default
module load angsd/0.940

#Choose individual
SAMPLE="SE2109" #CH0893 and SE2109

#Choose the output directory
output_directory="/home/x-YOUR_USERNAME/heterozygosity"

#Set variables. Do not change this part.
input_bam_file="/home/x-larantes/00_raw_data/${SAMPLE}.loxAfr4_NC000934"
ancestral_fasta_file="/home/x-larantes/00_raw_data/loxAfr4.fasta"
reference_fasta_file="/home/x-larantes/00_raw_data/loxAfr4.fasta"

#Loop through scaffolds 11 to 14

for chr in {11..14}; do

    #Run ANGSD
    angsd -P 1 -i ${input_bam_file}_${chr}.bam -anc ${ancestral_fasta_file} -dosaf 1 -gl 1 -C 50 -minQ 20 -minmapq 30 -out ${output_directory}/$SAMPLE.chr_${chr} -ref ${reference_fasta_file}

    #Run realSFS
    realSFS -nsites 200000 ${output_directory}/$SAMPLE.chr_${chr}.saf.idx -fold 1 > ${output_directory}/$SAMPLE.chr_${chr}.est.ml

done
```

Replace the `<YOUR_USERNAME>` with your Anvil username, define the path where the output files will be saved, and choose a `$SAMPLE` variable to start running it.
Because of our time and resource limitations, we are running this analysis only for 4 short chromosomes, from chr 11 to chr 14.
This script will loop through chromosomes 11 to 14, running the ANGSD command and then the realSFS command for each chromosome.
ANGSD will calculate the site allele frequency (SAF) likelihood based on individual genotype likelihoods assuming HWE. It generates 3 output binary files (angsdput.saf.idx, angsdput.saf.pos.gz and angsdput.saf.gz).
realSFS will run an optimization of the .saf file, which will estimate the Site Frequency Spectrum (SFS). We are estimating it in windows of a fixed number of sites (200,000 bp). It generates an output file .est.ml, which we will further manipulate to plot the results.

Make sure you understand each option in the ANGSD command line:

-   `P <threads>` \- Sets the number of threads to be used in parallel.
-   `i <input_bam_file>` \- Specifies the input file as a BAM file.
-   `anc <ancestral_fasta_file>` \- Specifies the ancestral fasta reference file.
-   `dosaf <dosaf_value>` \- Computes the Site Frequency Spectrum (SFS) based on the genotype likelihoods.
-   `gl <genotype_likelihood_method>` \- Specifies the method used for calculating genotype likelihoods.
-   `C <base_quality_adjustment>` \- Adjusts the base quality score by a specified value before using it.
-   `minQ <min_base_quality>` \- Sets the minimum base quality score required.
-   `minmapq <min_mapping_quality>` \- Sets the minimum mapping quality score required.
-   `fold <fold_value>` \- Indicates whether you are analyzing folded SFS or unfolded SFS.
-   `out <output_file>` \- Specifies the output file path and name.
-   `ref <reference_fasta_file>` \- Specifies the reference fasta file.
-   `r <region_of_interest>` \- Specifies the region of interest. ANGSD does not differentiate the chromosomes if you run the whole genome at once, that is why we need to use the ‘region’ variable when running ANGSD to specify the chromosome/scaffold. This will be useful to look the heterozygosity throughout the genome.

Also, make sure you understand each option in the realSFS command line:

-   `nsites <number_of_sites>` \- Specifies the number of sites to be considered for the estimation. 
-   `<input_saf_idx_file>` \- Specifies the input .saf.idx file, which is the output from the ANGSD program that contains information about the SFS.
-   `> <output_est_ml_file>` \- Specifies the output file path and name for the maximum likelihood estimate of the SFS.


#### b) Prepare files and plot heterozygosity along the scaffold

Now we can add the sample name and scaffold number for each line in our output. This will make our work easier when we want to plot our results.

```bash
#!/bin/bash

output_directory="path/to/your/output_directory"
SAMPLE="your_sample_name"

# Loop through scaffolds 1 to 18
for i in {1..18}; do
  # Calculate the number of lines in the output file
  num_lines=$(wc -l < ${output_directory}/$SAMPLE.SUPER_${i}.est.ml)

  # Add the number of lines, sample name, and scaffold number to each line of the output file
  awk -v lines="$num_lines" -v sample="$SAMPLE" -v scaffold="$i" '{print lines, sample, "SUPER_" scaffold, $0}' ${output_directory}/$SAMPLE.SUPER_${i}.est.ml > ${output_directory}/$SAMPLE.SUPER_${i}.est.ml.annotated

  # Optional: Move the annotated file to the original file
  mv ${output_directory}/$SAMPLE.SUPER_${i}.est.ml.annotated ${output_directory}/$SAMPLE.SUPER_${i}.est.ml

done
```

Concatenate files

```bash
#!/bin/bash

Set the input and output directory

input_directory="path/to/your/input_directory"
output_directory="path/to/your/output_directory"

Create a new output file

output_file="${output_directory}/all_est_ml_concatenated.txt"

Remove the output file if it already exists

if [ -f "${output_file}" ]; then
rm "${output_file}"
fi

Loop through all "[est.ml](http://est.ml/)" files and concatenate them
for file in ${input_directory}/*.est.ml; do
cat "${file}" >> "${output_file}"
done
```

Plot the results for one chromosome using R

```r
# Load required libraries
library(tidyverse) # Collection of packages for data manipulation and visualization
library(viridis)   # Package for generating color palettes
library(scales)    # Package for scaling and formatting axis labels

# Read data file and store it in the variable 'het_master'
het_master <- read.table("/path/to/file/NN_SUPER.est.ml")

# Manipulate the data
het_master %>% 
  rename(sample=V2) %>%          # Rename V2 as 'sample'
  rename(chromosome = V3) %>%    # Rename V3 as 'chromosome'
  mutate(heterozygosity = V5/(V4 + V5)) %>% # Calculate heterozygosity as V5 / (V4 + V5)
  mutate(position = ((V1*200000)-200000))   %>% # Calculate position as (V1 * 200000) - 200000
  filter(chromosome == "SUPER_2") %>% # Filter data to include only rows where 'chromosome' is 'SUPER_2'

  # Create a ggplot2 plot
  ggplot(aes(x=position, y=heterozygosity)) + # Set x-axis as 'position' and y-axis as 'heterozygosity'
  geom_line(colour="grey",alpha=0.5) + # Add a line plot with grey color and 0.5 alpha (transparency)
  geom_point(aes(colour=factor(chromosome))) + # Add points, color them based on 'chromosome' factor variable
  scale_color_viridis(discrete = TRUE) + # Use viridis color palette for the points
  facet_grid(sample ~ chromosome,scales = "free_x") + # Create a facet grid with 'sample' on the y-axis and 'chromosome' on the x-axis, set x-axis scales to be free
  labs(x = NULL, y = "Heterozygosity\n") + # Remove x-axis labels and set y-axis label to "Heterozygosity\n"
  scale_y_continuous(labels = comma) + # Format y-axis labels with commas
  scale_x_continuous(labels = comma) + # Format x-axis labels with commas
  theme_minimal() + # Apply a minimal theme to the plot
  theme(legend.position = "none", # Remove legend
        strip.text.x = element_text(face = "bold"), # Set strip text for x-axis to bold
        strip.text.y = element_text(face = "bold"), # Set strip text for y-axis to bold
        panel.grid.major.x = element_blank(), # Remove major x-axis gridlines
        panel.grid.minor.x = element_blank(), # Remove minor x-axis gridlines
        panel.spacing.x = unit(0, "line"), # Set panel spacing to zero
        panel.border = element_rect(color = "black", fill = NA, size = 0.25)) # Add a black border around the panels
```

> [!IMPORTANT]
> :elephant::grey_question: How is heterozygosity distributed along the scaffold? Is the distribution uniform, or are there noticeable regions with higher or lower heterozygosity? Are there gaps or regions with no detectable heterozygosity?



