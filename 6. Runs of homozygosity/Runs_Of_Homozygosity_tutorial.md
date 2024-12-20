# ConGen 2025 - Runs of homozygosity tutorial

Runs of Homozygosity (ROH) are defined as the uninterrupted stretches of homozygous genotypes within an individual's genome. These regions can provide valuable insights into the demographic history, inbreeding levels, and disease susceptibility of a population. By analyzing the length and distribution of ROH, we can infer the population structure, migration patterns, and effective population size of a group.

Moreover, ROH can also be used to identify deleterious mutations and genomic regions under positive selection. Inbreeding depression, which is caused by the accumulation of deleterious alleles, can be estimated by measuring the frequency and length of ROH. Longer ROH segments are associated with increased homozygosity and reduced genetic diversity, which can lead to reduced fitness and increased risk of disease.

In this tutorial, we will demonstrate how to estimate ROH using the 'bcftools roh' plugin, which is a widely used tool for detecting ROH from VCF files. We will also discuss the interpretation and application of ROH results in different research contexts, such as conservation genetics, human population genetics, and animal breeding. By the end of this tutorial, you will have a better understanding of the biological significance and practical utility of ROH analysis.

Please see below the location of the VCF file that we will use for the tutorial.

```r
/pool/genomics/figueiroh/SMSC_2023_class/vcf/NN_6samples_HD_PASS_DP5.vcf.gz
```

For the RoHs plots, we are only interested in the “RG” ???? portion of the files, where it contains the homozygous blocks in the genome. These blocks are important because they are indicative of long stretches of DNA that are identical in the two chromosomes, which can occur when the parents are related.

There are several ways of showing the results, and it will mostly depend on your main question. For example, if you are interested in the frequency of RoHs in different populations, you can create histograms that show the distribution of the length of these blocks. On the other hand, if you want to study the relationship between RoHs and disease, you may want to compare the number and length of RoHs between cases and controls, and perform statistical tests to determine if there is an association. In either case, it is important to consider the study design and the underlying biological mechanisms that could affect the results.

```bash
#!/bin/bash
# FILENAME:  roh
#SBATCH -A bio240351  # Allocation name
#SBATCH --nodes=1         # Total # of nodes (must be 1 for serial job)
#SBATCH --ntasks=1        # Total # of MPI tasks (should be 1 for serial job)
#SBATCH --time=1:30:00    # Total run time limit (hh:mm:ss)
#SBATCH -J roh     # Job name
#SBATCH -o /home/x-YOUR_USERNAME/logs/roh.o%j      # Name of stdout output file
#SBATCH -e /home/x-YOUR_USERNAME/logs/roh.e%j      # Name of stderr error file
#SBATCH -p wholenode  # Queue (partition) name
#SBATCH --mail-user=x-YOUR_USERNAME@anvil.rcac.purdue.edu
#SBATCH --mail-type=all   # Send email to above address at begin and end of job

#Load bcftools
module load biocontainers/default
module load bcftools/1.17

#Define path to input files - do not change it
VCF_FILE="/home/x-larantes/00_raw_data/elephants.chr.11.12.subset.vcf.gz"

#Define path to save output files
OUTPUT_PATH="/home/x-YOUR_USERNAME/roh"
mkdir ${OUTPUT_PATH}

bcftools roh --AF-dflt 0.4 -I -G30 --rec-rate 1.4e-9 ${VCF_FILE} > ${OUTPUT_PATH}/elephant.roh.txt
```

-   `bcftools roh`: This command runs the `roh` plugin from `bcftools` to detect runs of homozygosity.
-   `-AF-dflt 0.4`: This option sets the default allele frequency to 0.4 when the frequency is missing in the input file.
-   `I`: This option tells the plugin to perform the imputation of missing genotypes.
-   `G30`: This option sets the phred-scaled genotype quality threshold to 30. Genotypes below this quality threshold will be treated as missing.
-   `-rec-rate 1.4e-9`: This option sets the recombination rate to 1.4e-9.

After running the 'bcftools roh' plugin, you may want to filter and process the output file to retain specific information. For example, you can extract the chromosome, start, and end positions of ROH using the following command:

```bash
grep "RG" elephant.roh.txt | cut -f 2,3,6 > elephant.roh.edited.txt
```

You can run this command on interactive mode.

Two of the most basic statistics we can obtain from this analysis are the number of runs of homozygosity blocks (NROH) and the total length of genome that showed runs of homozygosity (SROH). These two values provide highly informative data, and the ratio between them is known as the inbreeding coefficient (FROH).

You can use the following R script to estimate these values and plot the results.

```r
# Set the working directory
setwd("C:/Users/laris/Documentos/Work/Disciplinas_Courses/2024_ConGen_Namibia/Tutorials")

# Load libraries and read data
library(tidyverse)
library(ggrepel)

# Read data with read_delim() for better control over input file parsing
roh <- read_delim("elephant.roh.with_population.txt", delim = "\t", skip = 1, col_names = c("Sample", "Chromosome", "RoH_length", "Population"))

# Compute NROH and SROH
nroh <- roh %>% 
  group_by(Sample) %>% 
  summarize(NROH = n())

sroh <- roh %>% 
  group_by(Sample) %>% 
  summarize(SROH = sum(RoH_length))

data <- inner_join(nroh, sroh, by = "Sample") %>% 
  inner_join(roh %>% select(Sample, Population) %>% distinct(), by = "Sample")

# Inspect the resulting dataset to ensure it is correct
head(data)

# Create the plot using the preprocessed dataset
snroh_plot <- data %>% 
  ggplot(aes(x = SROH, y = NROH, color = Population)) + # Map color to Population
  geom_point(size = 3) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "NROH vs. SROH", color = "Population") # Add a color legend title

# Print the plot
print(snroh_plot)
ggsave("sroh_nroh.png", snroh_plot, width = 8, height = 6, dpi = 300)

# Create RoH length categories
roh_categories <- roh %>% 
  mutate(
    RoH_category = case_when(
      RoH_length >= 500000 & RoH_length < 1000000 ~ "500 kb - 1 Mb",
      RoH_length >= 1000000 & RoH_length < 2000000 ~ "1 Mb - 2 Mb",
      RoH_length >= 2000000 ~ "> 2 Mb",
      TRUE ~ NA_character_ # Assign NA to lengths < 500k
    )
  ) %>% 
  filter(!is.na(RoH_category)) # Remove rows where RoH_category is NA

# Summing RoH_length by Sample, RoH_category, and Population
summed_roh <- roh_categories %>%
  group_by(Sample, RoH_category, Population) %>%
  summarise(total_RoH_length = sum(RoH_length), .groups = "drop")  # Sum the RoH_length for each group

# Create the stacked bar plot
stacked_barplot <- summed_roh %>%
  ggplot(aes(x = Sample, y = total_RoH_length/1000000, fill = RoH_category)) +  # Stack by RoH_category
  geom_bar(stat = "identity") +
  facet_wrap(~ Population, nrow =1, scales = "free_x") +  # Facet by Population
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    strip.text = element_text(angle = 90, size = 10)) +  # Emphasize facet labels
  scale_y_continuous(labels = scales::comma) +  # Format y-axis labels as regular numbers with commas
  labs(x = "Sample", y = "Total RoH Length (Mb)", fill = "RoH Category")

# Print the plot
print(stacked_barplot)
ggsave("stacked_barplot.png", stacked_barplot, width = 8, height = 6, dpi = 300)

#Estimate FROH

# Define the total chromosome length in the analysis
total_chromosome_length <- 152468515

# Calculate FROH for each sample
frohs <- roh_categories %>%
  filter(RoH_length > 1000000) %>%  # Only include RoH_length > 1 Mb
  group_by(Sample) %>%              # Group by Sample
  summarise(total_RoH_above_1Mb = sum(RoH_length)) %>%  # Sum the RoH_length for each sample
  mutate(FROH = total_RoH_above_1Mb / total_chromosome_length)  # Calculate FROH

# Join it with the original roh_categories to retain Population information
frohs_with_population <- frohs %>%
  left_join(roh_categories %>% select(Sample, Population) %>% distinct(), by = "Sample")

# Reorder the Sample factor by Population
frohs_with_population$Sample <- factor(frohs_with_population$Sample, 
                                       levels = unique(frohs_with_population$Sample[order(frohs_with_population$Population)]))

# Create a boxplot per population with individual dots per sample
boxplot_froh <- frohs_with_population %>%
  ggplot(aes(x = Population, y = FROH, color = Population)) +  # Map Population to x-axis and color
  geom_boxplot(outlier.shape = NA, fill = "transparent", color = "black") +  # Boxplot without outliers shown
  geom_jitter(aes(color = Population), size = 3, width = 0.2) +  # Add individual dots per sample
  theme_minimal() +  # Minimal theme for clean visualization
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    axis.title.x = element_blank(),  # Optional: Remove x-axis title
    axis.title.y = element_text(size = 12),  # Set y-axis title size
    legend.position = "none" ) + # Remove legend, as colors are already mapped to Population
  labs(y = "FROH", x = "Population")

# Print the plot
print(boxplot_froh)
ggsave("boxplot_froh.png", boxplot_froh, width = 8, height = 6, dpi = 300)
```

> [!IMPORTANT]
>:elephant::grey_question: Compare inbreeding level between savanna and forest elephants. When did the inbreeding happen?


