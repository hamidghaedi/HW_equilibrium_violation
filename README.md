# HW equilibrium violation in human exomes
Identifying regions in exomes with excess of heterozygotes 

#### Note : This is an AI-assisted content


Within the human exome, specific regions exhibit a notable abundance of variants characterized by an excess of heterozygotes, ultimately resulting in a violation of Hardy-Weinberg equilibrium. GnomAD employs the inbreeding coefficient as a metric to identify instances where Hardy-Weinberg equilibrium is breached, traditionally associating this violation with inbreeding.

While the inbreeding coefficient effectively detects regions influenced by inbreeding, it also serves as an indicator for areas characterized by an elevated level of heterozygotes and a marked depletion of homozygotes. Intriguingly, these regions are recognized for their suboptimal sequencing coverage, particularly in exome sequencing. The consequence is that variants within these regions, despite potentially harboring somatic mutations, may evade conventional filtration steps applied in analytical pipelines. These variants often pass through different filtering stages designed to exclude spurious variant calls.

The objective of this study is to systematically identify and catalog these stretches of DNA characterized by an excess of heterozygotes and a deficiency of homozygotes. By doing so, we aim to shed light on the unique characteristics of these regions, uncover their potential impact on variant calling, and contribute to a comprehensive understanding of genomic variations within the human exome.


## 1. Data download and vcf filteration

```bash
#!/bin/bash
#SBATCH --job-name=gnomAD
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=2-00:00:00
#SBATCH --mem=50gb
#SBATCH --output=gnomad.%J.out
#SBATCH --error=gnomad.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=*

module load StdEnv/2020 gcc/9.3.0
module load vcflib/1.0.3  bcftools/1.16 bedtools/2.31.0


cd /home/ghaedi/projects/def-gooding-ab/ghaedi/gnomad

# Define the BED file
# the gtf2bed.ipynb shows the process of making this file
# Then used bed tools recommendation to sort the bed file : sort -k1,1 -k2,2n exome_coordinates.bed > in.sorted.bed
BED_FILE="in.sorted.bed"


for chromosome in {1..22} X; do
    # Check if the exome VCF file already exists
    EXOME_FILE="gnomad.exome.v4.0.sites.chr${chromosome}.vcf.bgz"
    if [ ! -f "$EXOME_FILE" ]; then
        echo "Downloading exome VCF file for chromosome $chromosome..."
        wget "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/$EXOME_FILE" -O "$EXOME_FILE"
    fi

    # Step 2: Filter exome VCF file based on inbreeding coefficient
    echo "Filtering exome VCF file for chromosome $chromosome based on inbreeding coefficient..."
    bcftools view -i 'INFO/inbreeding_coeff <= -0.3' "$EXOME_FILE" -Oz -o "chr${chromosome}_exome_filtered_by_InbreedCoef.vcf.gz"

    # Step 3: Intersect exome VCF file with exons
    echo "Intersecting exome VCF file for chromosome $chromosome with exons..."
    bedtools intersect -a "chr${chromosome}_exome_filtered_by_InbreedCoef.vcf.gz" -b "$BED_FILE" -sorted -header > "chr${chromosome}_exome_intersected.vcf"

    # Step 4: Convert exome VCF file to TSV
    echo "Converting exome intersected VCF file for chromosome $chromosome to TSV..."
    vcf2tsv "chr${chromosome}_exome_intersected.vcf" > "chr${chromosome}_exome_intersected.tsv"
done

    # Check if the genome VCF file already exists
    GENOME_FILE="gnomad.genome.v4.0.sites.chr${chromosome}.vcf.bgz"
    if [ ! -f "$GENOME_FILE" ]; then
        echo "Downloading genome VCF file for chromosome $chromosome..."
        wget "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes/$GENOME_FILE" -O "$GENOME_FILE"
    fi

    # Step 2: Filter genome VCF file based on inbreeding coefficient
    echo "Filtering genome VCF file for chromosome $chromosome based on inbreeding coefficient..."
    bcftools view -i 'INFO/inbreeding_coeff <= -0.3' "$GENOME_FILE" -Oz -o "chr${chromosome}_genome_filtered_by_InbreedCoef.vcf.gz"

    # Step 3: Intersect genome VCF file with exons
    echo "Intersecting genome VCF file for chromosome $chromosome with exons..."
    bedtools intersect -a "chr${chromosome}_genome_filtered_by_InbreedCoef.vcf.gz" -b "$BED_FILE" -sorted -header > "chr${chromosome}_genome_intersected.vcf"

    # Step 4: Convert genome VCF file to TSV
    echo "Converting genome intersected VCF file for chromosome $chromosome to TSV..."
    vcf2tsv "chr${chromosome}_genome_intersected.vcf" > "chr${chromosome}_genome_intersected.tsv"

    echo "Processing for chromosome $chromosome completed."
done

```
## 2. Processing the filtered vcf files


Then the * _genome_intersected.vcf were moved to another directory and used the following chunk to process the files :

```bash 
# cd to a directory where all *_exome_intersected.vcf have been moved
############################ EXOME #################################################
# List of chromosome numbers
chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

# Output directory for results
output_dir="exome_result"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Output file for all intervals with length
all_intervals_file="$output_dir/all_intervals_with_length.bed"

# Loop through each chromosome
for chromosome in "${chromosomes[@]}"; do
    # Exome file
    exome_file="chr${chromosome}_exome_intersected.vcf"
    
    # Merge adjacent intervals with a 300 bp gap
    bedtools merge -i "$exome_file" -d 300 -c 1 -o count > "$output_dir/testchr${chromosome}.bed"
    
    # Remove intervals with 1 variant
    awk '$4 != 1' "$output_dir/testchr${chromosome}.bed" > "$output_dir/filtered_intervals_chr${chromosome}.bed"
    
    # Add a column for interval length
    awk '{print $0, $3 - $2}' "$output_dir/filtered_intervals_chr${chromosome}.bed" > "$output_dir/intervals_with_length_chr${chromosome}.bed"
    
    # Append to the output file
    cat "$output_dir/intervals_with_length_chr${chromosome}.bed" >> "$all_intervals_file"
    
    echo "Processing for exome chromosome $chromosome completed."
done

echo "All intervals with length have been appended to $all_intervals_file."


############################ GENOME #################################################
# List of chromosome numbers
chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

# Output directory for results
output_dir="genome_result"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Output file for all intervals with length
all_intervals_file="$output_dir/all_intervals_with_length.bed"

# Loop through each chromosome
for chromosome in "${chromosomes[@]}"; do
    # Genome file
    genome_file="chr${chromosome}_genome_intersected.vcf"
    
    # Merge adjacent intervals with a 300 bp gap
    bedtools merge -i "$genome_file" -d 300 -c 1 -o count > "$output_dir/testchr${chromosome}.bed"
    
    # Remove intervals with 1 variant
    awk '$4 != 1' "$output_dir/testchr${chromosome}.bed" > "$output_dir/filtered_intervals_chr${chromosome}.bed"
    
    # Add a column for interval length
    awk '{print $0, $3 - $2}' "$output_dir/filtered_intervals_chr${chromosome}.bed" > "$output_dir/intervals_with_length_chr${chromosome}.bed"
    
    # Append to the output file
    cat "$output_dir/intervals_with_length_chr${chromosome}.bed" >> "$all_intervals_file"
    
    echo "Processing for genome chromosome $chromosome completed."
done

echo "All intervals with length have been appended to $all_intervals_file."

```
The vcf files, filtred , process and result files are placed in the following directories :

```bash
$ tree -d
.
├── exome_intersected_tabular
├── exome_vcf
├── filtered_InbreedingCoef_vcf
├── genome_intersected_tabular
├── genome_vcf
└── intersected_vcf
    ├── exome_result
    └── genome_result
```
Directories `exome_result` and `genome_result` were downloaded for further processing in R enviornment

3) Statistical analysis

```R
```
There are # exHet occurring in human coding regions based on exome variants from gnomAd v4. These variants are not evenly scattered throu human exome and there are genomic regions enriched by such variants.  We identified a total of 631 exonic intervals (min size(bp) = , max size(bp) = , mean size(bp) = ) that harbor minumum of two  and maximum of # exHet variants. For further analysis we were focused on regions with minum 100 bp length or at least 5 mutation and  defined variant density score as number of variants divided by genomic interval length (bp). 

The number of intervals per chromosome is/ is not significantly diffrent (? test, p-value = ). Also there is / is not a coorelatiuon between chrosome length and the sum of teh identified intervals.
We showed that these intervals are significantly enriched for exHet variants in comparsion to exonic regions out of these intervals. 

We found that ssuch intervals are more frequcntly associated with gen/gene families among others (? test, p-vale = ?)



