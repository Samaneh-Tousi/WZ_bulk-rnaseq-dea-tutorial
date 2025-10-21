# bulk RNA-seq and differential expression analysis (DEA) on VSC
## Step 1 — Access to VSC and Starting Your Session
Before running any analyses, every student must have a valid
VSC (Vlaamse Supercomputer Centrum) account and access to an
interactive compute environment.
### Request your VSC account

Each student must first apply for an account on the Flemish Supercomputer infrastructure:
[Apply for your VSC account](https://docs.vscentrum.be/accounts/vsc_account.html#applying-for-your-vsc-account)
Follow the steps described on that page carefully.
When asked for your **institution**, select **Hasselt University**.
### Request your introduction credits

After obtaining a VSC account, request your initial computing credits (needed to run jobs):
[Request Introduction Credits (KU Leuven / UHasselt)](https://admin.kuleuven.be/icts/onderzoek/hpc/request-introduction-credits)
You will receive a confirmation email once your credits are approved and available.
### Login to the VSC OnDemand portal

Go to: [OnDemand](https://ondemand.hpc.kuleuven.be/)
1. Click **“Partner Organizations: VSC Account.”**
2. Choose **Hasselt University** from the dropdown.
3. Log in with your **VSC username** and **password**.
4. You will be redirected to the **VSC Dashboard**.
### Start an Interactive Shell Session

To get a working Linux terminal environment for this course:
1. Go to the **“Interactive Apps”** tab.
2. Select **“Interactive Shell.”**
3. Configure the following settings:
| Setting                | Value                                          |
|------------------------|------------------------------------------------|
| Cluster                | Genius                                         |
| Partition              | interactive                                    |
| VSC Account            | your personal VSC account                      |
| Walltime               | 2 hours                                        |
| Nodes                  | 1                                              |
| Processors per node    | 8                                              |
| Cores per task         | 8                                              |
| Memory per core        | 7500 MB                                        |
| Reservation            | Bioinfo_course                                 |
| Email notification     | ✓ “I would like to receive an email when the session starts” |
4. Click **Launch** to start your session.
When ready, click **“Connect to VS Code / Jupyter / Shell”** to open your environment.
## Step 2 — Connect to reserved cores, prepare workspace, load modules, and fetch data
After launching your Interactive Shell session on VSC, connect to the terminal,
create your course workspace under $VSC_DATA, load required modules, and download
the SRA data for the microglia study (GSE111972).
### Connect to your reserved cores (shell/SSH)

In the OnDemand interface for your running job, click the **Connect** button at the bottom
(e.g., “Connect to Shell/SSH”). This opens a terminal on the reserved compute node(s).
### Create your course folder under $VSC_DATA

Use the VSC data path to store course files. Replace the folder name if preferred.
```bash
cd "$VSC_DATA"
mkdir -p Bioinfo_course
cd Bioinfo_course
pwd  # verify your working directory
```
### Load required environment modules

Load the interactive environment and SRA Toolkit (version may vary on your cluster).
```bash
module load cluster/genius/interactive
module load SRA-Toolkit/3.0.5-gompi-2021a
```
> If a module is unavailable, run `module avail` to see alternatives or contact support.
### Study data: Microglia in MS (GSE111972)

Source: Transcriptional profiling of human microglia reveals grey–white matter
heterogeneity and multiple sclerosis-associated changes.
GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111972
**Selected single-end runs (subset for the course):**
| sample_id  | GSM        | SRR       | reads_millions | layout     | approx_size | tissue          | condition |
|------------|------------|-----------|----------------|------------|-------------|-----------------|-----------|
| MS_WM_1    | GSM3045818 | SRR6849240| 33.4           | Single-end | 891 MB      | MS white matter | MS        |
| MS_WM_10   | GSM3045819 | SRR6849241| 31.3           | Single-end | 842 MB      | MS white matter | MS        |
| MS_WM_11   | GSM3045820 | SRR6849242| 26.7           | Single-end | 698 MB      | MS white matter | MS        |
| CON_WM_1   | GSM3045833 | SRR6849255| 29.1           | Single-end | 756 MB      | Control white matter | Control |
| CON_WM_10  | GSM3045834 | SRR6849256| 35.0           | Single-end | 905 MB      | Control white matter | Control |
| CON_WM_11  | GSM3045835 | SRR6849257| 27.5           | Single-end | 716 MB      | Control white matter | Control |
> Note: Data are **single-end**; do **not** use paired-end flags when processing.
### Download FASTQ files from SRA (adjust output path!)

Set an **absolute path** where your FASTQs will be saved (edit OUT=... to your own VSC path),
then run the loop. This uses `prefetch` + `fasterq-dump` and compresses with `pigz` if available,
else falls back to `gzip`.
```bash
OUT="$VSC_DATA/Bioinfo_course/MS_microglia_fastq"
mkdir -p "$OUT"
for SRR in SRR6849240 SRR6849241 SRR6849242 SRR6849255 SRR6849256 SRR6849257; do
echo "Downloading $SRR ..."
prefetch "$SRR" \
&& fasterq-dump -O "$OUT" -e 8 "$SRR" \
&& pigz -p 8 "$OUT/${SRR}.fastq" \
|| gzip -9 "$OUT/${SRR}.fastq"
done
```
**Tips**
- Ensure you have sufficient quota in your `$VSC_DATA`.
- For single-end data, `--split-files` is **not** needed.
- If `pigz` is not installed as a module, `gzip` fallback will be used.
### (Optional) Organize or link into your repository

If you have a Git repository for the tutorial, you may symlink FASTQs into `data/raw`
(keeps the repo small, data outside the repo):
```bash
cd /path/to/your/repo
mkdir -p data/raw
ln -s "$OUT"/*.fastq.gz data/raw/
ls -lh data/raw
```
> Keep large data out of Git. Use symlinks or paths in your config files instead.
metadata:
dataset_accession: "GSE111972"
sra_runs:
- sample_id: "MS_WM_1"
GSM: "GSM3045818"
SRR: "SRR6849240"
condition: "MS"
tissue: "MS white matter"
layout: "SE"
reads_millions: 33.4
- sample_id: "MS_WM_10"
GSM: "GSM3045819"
SRR: "SRR6849241"
condition: "MS"
tissue: "MS white matter"
layout: "SE"
reads_millions: 31.3
- sample_id: "MS_WM_11"
GSM: "GSM3045820"
SRR: "SRR6849242"
condition: "MS"
tissue: "MS white matter"
layout: "SE"
reads_millions: 26.7
- sample_id: "CON_WM_1"
GSM: "GSM3045833"
SRR: "SRR6849255"
condition: "Control"
tissue: "Control white matter"
layout: "SE"
reads_millions: 29.1
- sample_id: "CON_WM_10"
GSM: "GSM3045834"
SRR: "SRR6849256"
condition: "Control"
tissue: "Control white matter"
layout: "SE"
reads_millions: 35.0
- sample_id: "CON_WM_11"
GSM: "GSM3045835"
SRR: "SRR6849257"
condition: "Control"
tissue: "Control white matter"
layout: "SE"
reads_millions: 27.5
## Step 3 — Downsampling FASTQ files to accelerate analysis
To make the downstream analysis faster, downsample each FASTQ file to 5 million reads
using `seqtk`. This step reduces runtime and storage usage while preserving biological patterns.
### Load seqtk module
|
```bash
module load seqtk
```
> If unavailable, check with `module avail seqtk` or ask your system admin.
### Downsample to 5 million reads

Define input and output directories (using the standard `$VSC_DATA` structure):
```bash
IN="$VSC_DATA/Bioinfo_course/MS_microglia_fastq"
OUT="$VSC_DATA/Bioinfo_course/MS_microglia_fastq_sub5M"
mkdir -p "$OUT"
```
Then run the loop to process all samples:
```bash
for s in SRR6849240 SRR6849241 SRR6849242 SRR6849255 SRR6849256 SRR6849257; do
echo "Subsampling $s to 5,000,000 reads ..."
seqtk sample -s42 "$IN/${s}.fastq.gz" 5000000 | gzip > "$OUT/${s}.sub5M.fastq.gz"
done
```
**Notes**
- `-s42` sets a random seed (42) to make subsampling reproducible.
- Adjust the number `5000000` if you want a different depth.
- The output FASTQs will be stored in `$OUT` and compressed with gzip.
## Step 4 — Quality control (QC) and trimming of sub-sampled FASTQ files
Perform quality control and adapter trimming on the 5M sub-sampled FASTQ files using fastp.
This step removes low-quality bases and short reads, generating cleaned FASTQs and HTML QC reports.
### Load fastp module

```bash
module load fastp/0.23.2-GCC-10.3.0
```
> Check `module avail fastp` if this exact version is unavailable on your cluster.
### Run fastp on all sub-sampled FASTQs

Define the input and output directories based on your `$VSC_DATA` structure:
```bash
IN="$VSC_DATA/Bioinfo_course/MS_microglia_fastq_sub5M"
OUT="$VSC_DATA/Bioinfo_course/MS_microglia_fastp"
mkdir -p "$OUT"
```
Then run the following loop to trim all samples:
```bash
for fq in "$IN"/SRR*.fastq.gz; do
s=$(basename "${fq%.fastq.gz}")
echo "Running fastp on $s ..."
fastp -i "$fq" \
-o "$OUT/${s}.trimmed.fastq.gz" \
-h "$OUT/${s}_fastp.html" \
-w 8
done
```
**Notes**
- `-i` specifies the input FASTQ file (single-end).
- `-o` sets the output trimmed FASTQ file.
- `-h` generates an interactive HTML report for each sample.
- `-w 8` uses 8 threads for faster processing.
- Output files and reports are stored in `$OUT`.
## Step 5 — Building the reference genome and aligning trimmed FASTQ files
Download the human reference genome (GRCh38), build a STAR index, and prepare for alignment
of trimmed FASTQ files. This step sets up the reference needed for mapping reads to genes.
### Create a reference genome directory

```bash
cd "$VSC_DATA/Bioinfo_course"
mkdir -p Ref_genome
cd Ref_genome
```
### Load the STAR module

```bash
module load STAR/2.7.3a-GCCcore-6.4.0
```
> Check `module avail STAR` if this version differs on your cluster.
### Download the reference genome and annotation (Ensembl GRCh38)

```bash
wget https://ftp.ensembl.org/pub/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/current/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.chr.gtf.gz
gunzip *.gz
```
This will produce:
- `Homo_sapiens.GRCh38.dna.primary_assembly.fa`
- `Homo_sapiens.GRCh38.115.chr.gtf`
### Generate the STAR genome index

```bash
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir "$VSC_DATA/Bioinfo_course/Ref_genome" \
--genomeFastaFiles "$VSC_DATA/Bioinfo_course/Ref_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
--sjdbGTFfile "$VSC_DATA/Bioinfo_course/Ref_genome/Homo_sapiens.GRCh38.115.chr.gtf" \
--sjdbOverhang 100
```
**Notes**
- `--sjdbOverhang` should be read length minus 1 (e.g., 100 for 101 bp reads).
- The index will be stored inside `$VSC_DATA/Bioinfo_course/Ref_genome`.
- This step may take several minutes depending on resources.
### (Preview of next step)

Once the reference is ready, you will align your trimmed FASTQ files
located in `$VSC_DATA/Bioinfo_course/MS_microglia_fastp` against this genome
using STAR in Step 6.
## Step 6 — Aligning trimmed FASTQ files to the reference genome using STAR
Map the trimmed FASTQ reads to the human reference genome (GRCh38)
using STAR, producing sorted BAM files and gene-level read counts.
### Set directories and load STAR

```bash
module load STAR/2.7.3a-GCCcore-6.4.0
IDX="$VSC_DATA/Bioinfo_course/Ref_genome"
IN="$VSC_DATA/Bioinfo_course/MS_microglia_fastp"
OUT="$VSC_DATA/Bioinfo_course/MS_microglia_STAR_aligned"
mkdir -p "$OUT"
```
### Run STAR alignment on all trimmed FASTQ files

```bash
for fq in "$IN"/*.sub5M.trimmed.fastq.gz; do
s=$(basename "$fq" .sub5M.trimmed.fastq.gz)
echo "Aligning $s ..."
STAR --runThreadN 8 \
--genomeDir "$IDX" \
--readFilesIn "$fq" \
--readFilesCommand zcat \
--outFileNamePrefix "$OUT/${s}." \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts
done
```
**Explanation**
- `--runThreadN 8` uses 8 CPU threads per job.
- `--readFilesCommand zcat` decompresses gzipped FASTQ files on the fly.
- `--outFileNamePrefix` defines unique output filenames for each sample.
- `--outSAMtype BAM SortedByCoordinate` writes coordinate-sorted BAMs directly.
- `--quantMode GeneCounts` produces a gene-level count summary (used later in DESeq2).
### Output

STAR will create the following files for each sample inside `$OUT`:
- `<sample>.Aligned.sortedByCoord.out.bam` – aligned reads (sorted BAM)
- `<sample>.ReadsPerGene.out.tab` – gene-level counts
- `<sample>.Log.final.out` – summary alignment metrics
You can inspect alignment quality with:
```bash
less "$OUT"/<sample>.Log.final.out
```
## Step 7 — Verifying library strandness from STAR gene counts
Confirm the library orientation (strandness) using STAR’s ReadsPerGene output.
The NEBNext Ultra Directional RNA Library Prep Kit produces reverse-stranded libraries
(“RF”, “fr-firststrand”, or “-s 2”), but this can be verified empirically.
### Check strandness of aligned samples

Navigate to the STAR output directory:
```bash
cd "$VSC_DATA/Bioinfo_course/MS_microglia_STAR_aligned"
```
Run the following commands to calculate forward/reverse read fractions per sample:
```bash
OUT=strandness_from_STAR.tsv
printf "sample\tforward_counts(col3)\treverse_counts(col4)\tforward_frac\treverse_frac\n" > "$OUT"
for f in *.ReadsPerGene.out.tab; do
s=$(basename "$f" .ReadsPerGene.out.tab)
awk -v S="$s" 'NR>4 {f+=$3; r+=$4} END {t=f+r; if(t==0)t=1;
printf "%s\t%.0f\t%.0f\t%.3f\t%.3f\n", S, f, r, f/t, r/t}' "$f"
done >> "$OUT"
echo "Wrote $OUT"
cat "$OUT"
```
**Explanation**
- Column 3 = forward-strand counts
- Column 4 = reverse-strand counts
- The script calculates the fraction of reads aligning to each strand.
**Interpretation**
- If **reverse_frac ≫ forward_frac**, the library is **reverse-stranded** (expected for NEBNext Ultra Directional).
- If **forward_frac ≫ reverse_frac**, it would indicate a **forward-stranded** library.
- Use this result to set the correct strandness option in downstream counting or DE analysis (e.g., `-s 2` in featureCounts).
## Step 8 — Counting reads per gene using featureCounts (reverse-stranded libraries)
Quantify gene-level expression counts from the aligned BAM files using featureCounts.
Since these samples were prepared with the NEBNext Ultra Directional RNA-Seq kit,
the libraries are reverse-stranded (antisense orientation, `-s 2`).
### Load the Subread module

```bash
module load Subread
```
### Set paths and run featureCounts
|
```bash
GTF="$VSC_DATA/Bioinfo_course/Ref_genome/Homo_sapiens.GRCh38.115.chr.gtf"
ALIGN_DIR="$VSC_DATA/Bioinfo_course/MS_microglia_STAR_aligned"
OUTDIR="$VSC_DATA/Bioinfo_course/MS_microglia_featureCounts"
mkdir -p "$OUTDIR"
featureCounts -T 8 -s 2 -t exon -g gene_id \
-a "$GTF" \
-o "$OUTDIR/featureCounts_counts.txt" \
"$ALIGN_DIR"/*.Aligned.sortedByCoord.out.bam
```
**Explanation**
- `-T 8`: use 8 threads
- `-s 2`: strand-specific, **reverse-stranded** (as verified in Step 7)
- `-t exon`: count features of type “exon”
- `-g gene_id`: group by gene identifier from the GTF file
- `-a`: path to annotation file (GTF)
- `-o`: output count table file
- Input: all sorted BAMs in `$ALIGN_DIR`
### Convert to a clean count matrix (tab-delimited)

```bash
awk 'NR==2{
printf "gene"
for(i=7;i<=NF;i++){
  g=$i
  sub(/^.*\//,"",g)
  sub(/\\.Aligned\\.sortedByCoord\\.out\\.bam$/,"",g)
  printf "\\t" g
}
printf "\\n"
next
}
NR>2{
printf "%s", $1
for(i=7;i<=NF;i++) printf "\\t%s", $i
printf "\\n"
}' \
"$OUTDIR/featureCounts_counts.txt" > "$OUTDIR/featureCounts_counts_matrix.tsv"
```
**Outputs**
- Raw featureCounts summary:
`$VSC_DATA/Bioinfo_course/MS_microglia_featureCounts/featureCounts_counts.txt`
- Clean matrix for DESeq2 or downstream analysis:
`$VSC_DATA/Bioinfo_course/MS_microglia_featureCounts/featureCounts_counts_matrix.tsv`
## Step 9 — Switch to RStudio Server and run DEA (DESeq2) + GSEA
After generating featureCounts outputs, stop the Shell interactive job and start an
RStudio Server session (Interactive Apps) with ~4 cores. Perform DE analysis
(MS vs Control) using DESeq2, annotate genes, save DEGs, plot a volcano, and run GSEA.
### End Shell job and launch RStudio

1) **Terminate** the Shell interactive job.
2) In **OnDemand → Interactive Apps → RStudio Server**, request ~**4 cores**, launch, and connect.
### Run DESeq2 (MS vs Control), annotate, save outputs, and perform GSEA

In RStudio, run the following (adjust only if your paths differ):
```r
suppressPackageStartupMessages({
library(DESeq2); library(readr); library(dplyr); library(ggplot2); library(ggrepel); library(apeglm)
library(biomaRt)
library(clusterProfiler); library(msigdbr); library(org.Hs.eg.db); library(enrichplot)
})
# Paths
counts_path <- file.path(Sys.getenv("VSC_DATA"),
				 "Bioinfo_course/MS_microglia_featureCounts/featureCounts_counts_matrix.tsv")
out_dir <- file.path(Sys.getenv("VSC_DATA"), "Bioinfo_course/MS_microglia_DEA")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# Read counts (tab-delimited from Step 8)
counts <- read_tsv(counts_path)
gene <- counts$gene
cts <- as.data.frame(counts[,-1])
rownames(cts) <- gene
# Ensure sample order matches design
cts <- cts[, c("SRR6849240","SRR6849241","SRR6849242",
	   "SRR6849255","SRR6849256","SRR6849257")]
coldata <- data.frame(
row.names = colnames(cts),
condition = factor(c("MS","MS","MS","Control","Control","Control"))
)
# DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
# Wald results; (optionally use LFC shrinkage for plotting)
res <- results(dds, cooksCutoff = FALSE, independentFiltering = FALSE)
#res_shrunk <- lfcShrink(dds, coef = "condition_MS_vs_Control", type = "apeglm", res = res)
# ---- Annotation via Ensembl (biomaRt) ----
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_ids <- gsub("\\..*", "", rownames(res))
annot <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
	   filters = "ensembl_gene_id", values = gene_ids, mart = mart)
res_tbl <- as.data.frame(res)
res_tbl$ensembl_gene_id <- gsub("\\..*", "", rownames(res_tbl))
res_annot <- merge(res_tbl, annot, by = "ensembl_gene_id", all.x = TRUE)
res_annot <- res_annot[order(res_annot$padj), ]
# Save all results
write_tsv(res_annot,
  file.path(out_dir, "deseq2_results_annotated.tsv"))
# ---- Define DEGs and save ----
lfc_thr  <- 0.5
padj_thr <- 0.05
DEGs <- subset(res_annot, !is.na(padj) & padj < padj_thr & abs(log2FoldChange) > lfc_thr)
write_tsv(DEGs, file.path(out_dir, sprintf("DEGs_MS_vs_Control_padj%.2f_LFC%.2f.tsv", padj_thr, lfc_thr)))
# ---- Volcano plot (top 10 labeled) ----
df <- res_annot
if (!"external_gene_name" %in% names(df)) df$external_gene_name <- NA_character_
df$label <- df$external_gene_name
df$status <- "NotSig"
df$status[!is.na(df$padj) & df$padj < padj_thr & df$log2FoldChange >  lfc_thr] <- "Up"
df$status[!is.na(df$padj) & df$padj < padj_thr & df$log2FoldChange < -lfc_thr] <- "Down"
df$status <- factor(df$status, levels = c("NotSig","Up","Down"))
df$mlog10padj <- -log10(df$padj); df$mlog10padj[!is.finite(df$mlog10padj)] <- NA
sig <- df[!is.na(df$padj) & df$padj < padj_thr & abs(df$log2FoldChange) > lfc_thr, ]
ord <- order(sig$padj, -abs(sig$log2FoldChange))
top10 <- head(sig[ord, ], 10)
p <- ggplot(df, aes(x = log2FoldChange, y = mlog10padj, color = status)) +
geom_point(size = 1.3, alpha = 0.8, na.rm = TRUE) +
geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
ggrepel::geom_text_repel(
data = top10, aes(label = external_gene_name),
size = 3, max.overlaps = Inf, box.padding = 0.5, point.padding = 0.2
) +
scale_color_manual(values = c(NotSig = "grey70", Up = "red", Down = "blue")) +
labs(x = "log2 fold change",
y = expression(-log[10]("adjusted p-value")),
color = "Status",
title = "MS vs Control — Volcano plot",
subtitle = sprintf("FDR < %.02f & |LFC| > %.2f; top 10 labeled", padj_thr, lfc_thr)) +
theme_minimal(base_size = 12) +
theme(legend.position = "right", panel.grid.minor = element_blank())
ggsave(filename = file.path(out_dir, "volcano_DEGs.png"), p, width = 7, height = 5, dpi = 300)
# ---- GSEA (Hallmark gene sets) ----
df2 <- res_annot
symbols <- df2$external_gene_name
entrez  <- mapIds(org.Hs.eg.db, keys = symbols, keytype = "SYMBOL",
		  column = "ENTREZID", multiVals = "first")
df2$ENTREZID <- entrez
score <- df2$stat
if (is.null(score) || all(is.na(score))) score <- df2$log2FoldChange
rank_df <- tibble::tibble(ENTREZID = df2$ENTREZID, score = score) |>
dplyr::filter(!is.na(ENTREZID), is.finite(score)) |>
dplyr::group_by(ENTREZID) |>
dplyr::summarise(score = score[which.max(abs(score))], .groups = "drop")
ranks <- rank_df$score
names(ranks) <- rank_df$ENTREZID
ranks <- sort(ranks, decreasing = TRUE)
m_h <- msigdbr(species = "Homo sapiens", category = "H") |>
dplyr::select(gs_name, entrez_gene)
set.seed(42)
gseaH <- GSEA(ranks, TERM2GENE = m_h,
	  pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, verbose = FALSE)
# Save GSEA results and a summary plot
gsea_tbl <- as.data.frame(gseaH@result)
readr::write_tsv(gsea_tbl, file.path(out_dir, "GSEA_Hallmark_results.tsv"))
# Top pathways plot
png(file.path(out_dir, "GSEA_Hallmark_top4.png"), width = 900, height = 700)
print(gseaplot2(gseaH, geneSetID = 1:4, title = "Top Hallmark pathways"))
dev.off()
# Enrichment map (optional)
gseaH2 <- pairwise_termsim(gseaH)
png(file.path(out_dir, "GSEA_Hallmark_emap.png"), width = 1200, height = 900)
print(emapplot(gseaH2, showCategory = 10))
dev.off()
message("Outputs saved to: ", out_dir)
```
**Saved outputs (in `$VSC_DATA/Bioinfo_course/MS_microglia_DEA`):**
- `deseq2_results_annotated.tsv` — full DESeq2 results with gene names
- `DEGs_MS_vs_Control_padj0.05_LFC0.50.tsv` — significant DEGs (FDR<0.05 & |LFC|>0.5)
- `volcano_DEGs.png` — volcano plot
- `GSEA_Hallmark_results.tsv`, `GSEA_Hallmark_top4.png`, `GSEA_Hallmark_emap.png`