# Baca Lab Coding Exercise
## Ze Zhang 10/06/2023

## Task
* Obtain the files with fragment locations for the H3K4me3 cfChIP-seq experiments from this publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7610786/ (they are available in a Zenodo repository and can be downloaded without restricted access). (You only need to download the portion of the data needed for the following steps)
*	Select 10 samples from healthy donors and 10 samples from patients with colorectal cancer.
*	Identify peaks from these files using MACS2.
*	Merge these peaks into a single union set of peaks.
*	Use DESeq2 to identify peaks with stronger H3K4me3 cfChIP-seq signal in colorectal cancer compared to healthy at an FDR of < 0.01.
*	Create a “volcano plot” showing for each peak the -log(FDR) on the y axis and log2-fold difference between healthy and colorectal cancer on the x axis.

## Step 1
* Downloaded data from https://zenodo.org/record/3967254.
* Located .tagAlign files for H3K4me3.
* Selected 10 colorectal cancer samples and 10 healthy control samples.
* Created a metadata sheet with sampleID and tissue type. 

