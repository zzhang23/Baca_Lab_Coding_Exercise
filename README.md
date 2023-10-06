# Baca Lab Coding Exercise
## Ze Zhang 10/06/2023

## Task
* Obtain the files with fragment locations for the H3K4me3 cfChIP-seq experiments from this publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7610786/ (they are available in a Zenodo repository and can be downloaded without restricted access). (You only need to download the portion of the data needed for the following steps)
*	Select 10 samples from healthy donors and 10 samples from patients with colorectal cancer.
*	Identify peaks from these files using MACS2.
*	Merge these peaks into a single union set of peaks.
*	Use DESeq2 to identify peaks with stronger H3K4me3 cfChIP-seq signal in colorectal cancer compared to healthy at an FDR of < 0.01.
*	Create a “volcano plot” showing for each peak the -log(FDR) on the y axis and log2-fold difference between healthy and colorectal cancer on the x axis.

## Step 1 Data Download and Cleaning
* Downloaded data from https://zenodo.org/record/3967254.
* Located .tagAlign.gz files for H3K4me3.
* Selected 10 colorectal cancer samples and 10 healthy control samples and stored them in *Samples* folder.
* Created a metadata sheet with sampleID and tissue type [selected_samples.csv](https://github.com/zzhang23/Baca_Lab_Coding_Exercise/files/12831151/selected_samples.csv).
* Unzip all .tagAlign.gz files in *Samples* folder and stored them in *tagAlign* folder 
```
mkdir -p tagAlign
for file in Samples/*.gz; do
    filename=$(basename -- "$file")
    filename_noext="${filename%.*}"
    gzip -d -c "$file" > "tagAlign/$filename_noext"
done
```

## Step 2 MACS2 Peak Calling
* Use macs2 to perform peak calling on. tagAlign files in *tagAlign* folder and stored the outputs in *macs2_output* folder.
```
output_dir="macs2_output"
genome_size="3.2e9" 
mkdir -p "$output_dir"
for file in tagAlign/*; do
    filename=$(basename -- "$file")
    filename_noext="${filename%.*}"
    macs2 callpeak -t "$file" -g "$genome_size" --format AUTO -n "$output_dir/$filename_noext" --nomodel --extsize 147
done
```

## Step3 Data Conversion for DownStream Analysis
* Covert .tagAlign files in *tagAlign* folder to .bed files and stored them in *convertedBED* folder.
```
source_dir="tagAlign"
destination_dir="convertedBED"
mkdir -p "$destination_dir"
for tagalign_file in "$source_dir"/*.tagAlign; do
    base_filename=$(basename "$tagalign_file" .tagAlign)
    output_bed="$destination_dir/$base_filename.bed"
    awk '{print $0 "\t."}' "$tagalign_file" > "$output_bed"    
    echo "Converted $tagalign_file to $output_bed"
done
```
* Convert .bed files in *convertedBED* folder to. bam files and stored them in *convertedBAM* folder.
```
source_dir="convertedBED"
destination_dir="convertedBAM"
mkdir -p "$destination_dir"
genome_file="genome_file.txt"
for bed_file in "$source_dir"/*.bed; do
    base_filename=$(basename "$bed_file" .bed)
    output_bam="$destination_dir/$base_filename.bam"
    bedtools bedtobam -i "$bed_file" -g "$genome_file" > "$output_bam"
    echo "Converted $bed_file to $output_bam"
done
```
* Sort .bam files in *convertedBAM* folder and stored them in *sortedBAM* folder.
```
input_dir="convertedBAM"
output_dir="sortedBAM"
mkdir -p "$output_dir"
for bam_file in "$input_dir"/*.bam; do
    base_filename=$(basename "$bam_file" .bam)
    output_sorted_bam="$output_dir/$base_filename.sorted.bam"
    samtools sort -o "$output_sorted_bam" "$bam_file"
    echo "Sorted $bam_file to $output_sorted_bam"
done
```

## Step 4 Differential Analysis
* Create a sample sheet for differential analysis [samples_sheet.csv](https://github.com/zzhang23/Baca_Lab_Coding_Exercise/files/12831368/samples_sheet.csv).
* Run differential analysis in R to plot a volcano plot highlighting hits with FDR < 0.01 and absolute log fold change > 1.
```
library(DiffBind)
library(DESeq2)
library(ggplot2)
library(ggrepel)


samples<-read.csv("samples_sheet.csv")
directory_peaks <- "macs2_output"
directory_bams <- "sortedBAM"
sampleNames <- samples$SampleID
peakFiles <- file.path(samples$Peaks)
bamReads <- file.path(directory_bams, samples$bamReads)
head(peakFiles)
file.exists(head(peakFiles))
dba <- dba(sampleSheet = samples)

dba_count = dba.count(dba)

ds <- DESeqDataSetFromMatrix(countData = df,
                             colData = samples,
                             design = ~ Tissue)

df<-matrix(NA,nrow = nrow(dba_count$peaks[[1]]),ncol = 20)
rownames(df)<-paste(dba_count$peaks[[1]]$Chr,dba_count$peaks[[1]]$Start,dba_count$peaks[[1]]$End)
colnames(df)<-dba_count$samples$SampleID
for (i in 1:20) {
  df[,i]<-dba_count$peaks[[i]]$Reads
}

dds <- DESeq(dds)
result <- results(dds)

result_table <- as.data.frame(result)
result_table$location<-rownames(result_table)

ggplot(result_table, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(abs(log2FoldChange) > 1 & padj < 0.01, "Significant", "Not Significant")), size = 2) +
  geom_text_repel(data = result_table[order(result_table$padj),][1:3, ], aes(label = location), size = 3)+
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(x = "log2-fold change (CRC vs. Healthy)",
       y = "-log10(FDR)",
       title = "Volcano Plot of Differential Peaks") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +scale_x_continuous(breaks = seq(-5, 5, by = 1))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


```
![VolPlot](https://github.com/zzhang23/Baca_Lab_Coding_Exercise/assets/32206453/31680d41-f40c-4862-bbde-e0259ccedfa8)

