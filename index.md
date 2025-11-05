---
title: Bulk RNA-seq & DEA on VSC
description: A step-by-step, reproducible tutorial for bulk RNA-seq QC, alignment, counting, and DE analysis on Vlaamse Supercomputer Centrum (VSC).
---

# Bulk RNA-seq & Differential Expression Analysis (DEA) on VSC

> **Goal**: Run a complete bulk RNA-seq workflow on the **Vlaamse Supercomputer Centrum (VSC)**:  
> Sequence QC & trimming → alignment → counting → DESeq2 → GSEA.

---

## Contents
[Prerequisites](#prerequisites)  
[Step 1 - Access & Session](#step-1--access--session)  
[Step 2 - Connect, workspace, data](#step-2--connect-workspace-data)  
[Step 3 - Downsample FASTQ](#step-3--downsample-fastq)  
[Step 4 - QC & trimming (fastp)](#step-4--qc--trimming-fastp)  
[Step 5 - Reference genome (STAR index)](#step-5--reference-genome-star-index)  
[Step 6 - Alignment (STAR)](#step-6--alignment-star)  
[Step 7 - Strandness check](#step-7--strandness-check)  
[Step 8 - featureCounts (reverse-stranded)](#step-8--featurecounts-reverse-stranded)
[Step 9 - MultiQC summary report](#step-9--MultiQC summary report)  
[Step 10 - DESeq2 + GSEA (RStudio)](#step-10--deseq2--gsea-rstudio)  

---

## Pipeline at a glance

```mermaid
flowchart LR
  A[FASTQ (SRA)] --> B[Subsample (seqtk)]
  B --> C[QC & trimming (fastp)]
  C --> D[STAR Align]
  D --> E[Library Strandness check]
  E --> F[featureCounts (-s 2)]
  F --> G[DESeq2 (MS vs Control)]
  G --> H[GSEA (Hallmark)]
```  
 
  
# Prerequisites

- VSC account + intro credits or project credits
- Access to OnDemand and Interactive Apps
- Basic shell + R familiarity
- Storage under $VSC_DATA or staging storage

---
### What is an HPC Node?

An **HPC node** (High-Performance Computing node) is a single powerful computer within the larger **VSC cluster**.  
Each node has:
- Multiple **CPU cores** (e.g., 64,72, 96, or more)
- Large **memory (RAM)**
- Access to **high-speed storage** and **shared file systems**

When you start an **interactive session** on VSC, you’re temporarily reserving a portion of one or more nodes (e.g., 1 node with 8 cores and 7.5 GB RAM per core).  
Think of it like booking a workstation on the supercomputer for a limited time to run your analysis.
<img src="assets/HPC_nodes.png" alt="VSC HPC" width="600">

For more information about the VSC system and its usage, please see:
- [VSC Training Resources](https://www.vscentrum.be/vsctraining)
- [HPC Introduction (GitHub)](https://github.com/hpcleuven/HPC-intro/blob/master/HPCintro.pdf)

---

# Step 1 - Access to VSC & Interactive Sessions

[Apply for your VSC account](https://docs.vscentrum.be/accounts/vsc_account.html#applying-for-your-vsc-account), Select **Hasselt University** when prompted.

[Request introduction credits](https://admin.kuleuven.be/icts/onderzoek/hpc/request-introduction-credits)

**If you are a student enrolled in the Bioinformatics course, you do not need to request a VSC account, one has already been created for you and comes with the necessary project credits required to run the analyses.**

Login to OnDemand: [VSC OnDemand Portal](https://ondemand.hpc.kuleuven.be/)

<img src="assets/OnDemand.png" alt="VSC OnDemand login" width="600">

Start **Interactive Shell** with:

- Cluster: *Wice*
- Partition: *interactive*
- Account: *your intro_VSC account* or *your project account*
- Number of hours: *2 h* or more
- Number of nodes: *1*
- Number of processes per node: *8*
- Number of cores per each task: *8*
- Required memory per core in megabytes: *7500 MB*
- Reservation: *Bioinfo_course*
- ✅ I would like to receive an email when the session starts

# Step 2 - Connect, workspace, data

Connect to the running shell job (bottom Connect button). 
Once connected, use these commands to prepare your course workspace and load required software:

```
cd $VSC_DATA
mkdir -p Bioinfo_course && cd Bioinfo_course && pwd

# Load environment modules for this session
module load cluster/genius/interactive
module load SRA-Toolkit/3.0.5-gompi-2021a
```
Explanation:

**cd $VSC_DATA** → moves you into your personal storage directory on the cluster

**mkdir -p Bioinfo_course** → creates a folder for all course-related work

**cd Bioinfo_course** → enters that folder

**pwd** → confirms your current directory

**module load ...** → loads the tools needed for downloading and processing sequencing data

After this step, you’re ready to start downloading the example RNA-seq datasets in the next section.


## Dataset: GSE111972 (microglia, MS vs Control).
[GEO]https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111972

| sample_id | SRR        | condition | tissue               | layout     |
| --------- | ---------- | --------- | -------------------- | ---------- |
| MS_WM_1   | SRR6849240 | MS        | MS white matter      | Single-end |
| MS_WM_10  | SRR6849241 | MS        | MS white matter      | Single-end |
| MS_WM_11  | SRR6849242 | MS        | MS white matter      | Single-end |
| CON_WM_1  | SRR6849255 | Control   | Control white matter | Single-end |
| CON_WM_10 | SRR6849256 | Control   | Control white matter | Single-end |
| CON_WM_11 | SRR6849257 | Control   | Control white matter | Single-end |

Download:

```
OUT="$VSC_DATA/Bioinfo_course/MS_microglia_fastq"
mkdir -p "$OUT"
for SRR in SRR6849240 SRR6849241 SRR6849242 SRR6849255 SRR6849256 SRR6849257; do echo "Downloading $SRR ..."; prefetch "$SRR" && fasterq-dump -O "$OUT" -e 8 "$SRR" && pigz -p 8 "$OUT/${SRR}.fastq" || gzip -9 "$OUT/${SRR}.fastq"; done
```

# Step 3 - Downsample FASTQ

```
module load seqtk
IN="$VSC_DATA/Bioinfo_course/MS_microglia_fastq"
OUT="$VSC_DATA/Bioinfo_course/MS_microglia_fastq_sub5M"
mkdir -p "$OUT"

for s in SRR6849240 SRR6849241 SRR6849242 SRR6849255 SRR6849256 SRR6849257; do
  echo "Subsampling $s to 5,000,000 reads ..."
  seqtk sample -s42 "$IN/${s}.fastq.gz" 5_000_000 | gzip > "$OUT/${s}.sub5M.fastq.gz"
done
```

# Step 4 - QC & trimming (fastp)

```
module load fastp/0.23.2-GCC-10.3.0
IN="$VSC_DATA/Bioinfo_course/MS_microglia_fastq_sub5M"
OUT="$VSC_DATA/Bioinfo_course/MS_microglia_fastp"
mkdir -p "$OUT"

for fq in "$IN"/SRR*.fastq.gz; do
  s=$(basename "${fq%.fastq.gz}")
  echo "fastp $s"
  fastp -i "$fq" \
        -o "$OUT/${s}.trimmed.fastq.gz" \
        -h "$OUT/${s}_fastp.html" \
		-j "$OUT/${s}_fastp.json" \
        -w 8
done
```

# Step 5 - Reference genome (STAR index)

```
cd "$VSC_DATA/Bioinfo_course"
mkdir -p Ref_genome && cd Ref_genome
module load STAR/2.7.3a-GCCcore-6.4.0

wget https://ftp.ensembl.org/pub/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/current/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.chr.gtf.gz
gunzip *.gz

STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir "$VSC_DATA/Bioinfo_course/Ref_genome" \
     --genomeFastaFiles "$VSC_DATA/Bioinfo_course/Ref_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
     --sjdbGTFfile "$VSC_DATA/Bioinfo_course/Ref_genome/Homo_sapiens.GRCh38.115.chr.gtf" \
     --sjdbOverhang 100
	 
```
Note: sjdbOverhang = readLength − 1.

# Step 6 - Alignment (STAR)

```
module load STAR/2.7.3a-GCCcore-6.4.0
IDX="$VSC_DATA/Bioinfo_course/Ref_genome"
IN="$VSC_DATA/Bioinfo_course/MS_microglia_fastp"
OUT="$VSC_DATA/Bioinfo_course/MS_microglia_STAR_aligned"
mkdir -p "$OUT"

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
Check metrics:

```
less "$OUT"/<sample>.Log.final.out
```

# Step 7 - Strandness check
NEBNext Ultra Directional → reverse-stranded expected.

```
cd "$VSC_DATA/Bioinfo_course/MS_microglia_STAR_aligned"

OUT=strandness_from_STAR.tsv
printf "sample\tforward_counts(col3)\treverse_counts(col4)\tforward_frac\treverse_frac\n" > "$OUT"
for f in *.ReadsPerGene.out.tab; do
  s=$(basename "$f" .ReadsPerGene.out.tab)
  awk -v S="$s" 'NR>4 {f+=$3; r+=$4} END {t=f+r; if(t==0)t=1;
    printf "%s\t%.0f\t%.0f\t%.3f\t%.3f\n", S, f, r, f/t, r/t}' "$f"
done >> "$OUT"
cat "$OUT"
```
Interpretation: reverse fraction ≫ forward → use -s 2 downstream.

# Step 8 - featureCounts (reverse-stranded)

```
module load Subread
GTF="$VSC_DATA/Bioinfo_course/Ref_genome/Homo_sapiens.GRCh38.115.chr.gtf"
ALIGN_DIR="$VSC_DATA/Bioinfo_course/MS_microglia_STAR_aligned"
OUTDIR="$VSC_DATA/Bioinfo_course/MS_microglia_featureCounts"
mkdir -p "$OUTDIR"

featureCounts -T 8 -s 2 -t exon -g gene_id \
  -a "$GTF" \
  -o "$OUTDIR/featureCounts_counts.txt" \
  "$ALIGN_DIR"/*.Aligned.sortedByCoord.out.bam

awk 'NR==2{
        printf "gene"
        for(i=7;i<=NF;i++){
          g=$i; sub(/^.*\//,"",g); sub(/\.Aligned\.sortedByCoord\.out\.bam$/,"",g)
          printf "\t" g
        } printf "\n"; next
     }
     NR>2{
        printf "%s", $1; for(i=7;i<=NF;i++) printf "\t%s", $i; printf "\n"
     }' \
  "$OUTDIR/featureCounts_counts.txt" > "$OUTDIR/featureCounts_counts_matrix.tsv"
```

# Step 9 - MultiQC summary report

```
# (Optional) if MultiQC isn't preinstalled on VSC
conda install -c bioconda multiqc

# Define input and output directories
FASTP_DIR="$VSC_DATA/Bioinfo_course/MS_microglia_fastp"
STAR_DIR="$VSC_DATA/Bioinfo_course/MS_microglia_STAR_aligned"
FC_DIR="$VSC_DATA/Bioinfo_course/MS_microglia_featureCounts"
MQC_OUT="$VSC_DATA/Bioinfo_course/MS_microglia_MultiQC"

# Create output directory
mkdir -p "$MQC_OUT"

# Run MultiQC across all stages
multiqc -o "$MQC_OUT" -n multiqc_all "$FASTP_DIR" "$STAR_DIR" "$FC_DIR"

```

Outputs:

A single HTML summary:
$VSC_DATA/Bioinfo_course/MS_microglia_MultiQC/multiqc_all.html

A directory of summary data (multiqc_data/)

💡 Tip: Download and open the multiqc_all.html file locally to interactively browse read qualities, trimming stats, alignment rates, and featureCounts summaries.


# Step 10 - DESeq2 & GSEA (RStudio)

Switch app: Stop the shell job. Launch RStudio Server (4 cores).

```
suppressPackageStartupMessages({
  library(DESeq2); library(readr); library(dplyr); library(ggplot2); library(ggrepel); library(apeglm)
  library(biomaRt); library(clusterProfiler); library(msigdbr); library(org.Hs.eg.db); library(enrichplot)
})

counts_path <- file.path(Sys.getenv("VSC_DATA"),
                         "Bioinfo_course/MS_microglia_featureCounts/featureCounts_counts_matrix.tsv")
out_dir <- file.path(Sys.getenv("VSC_DATA"), "Bioinfo_course/MS_microglia_DEA")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# counts
counts <- read_tsv(counts_path)
gene <- counts$gene
cts <- as.data.frame(counts[,-1]); rownames(cts) <- gene
cts <- cts[, c("SRR6849240","SRR6849241","SRR6849242",
               "SRR6849255","SRR6849256","SRR6849257")]
coldata <- data.frame(row.names = colnames(cts),
                      condition = factor(c("MS","MS","MS","Control","Control","Control")))

dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, cooksCutoff = FALSE, independentFiltering = FALSE)

# annotate
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_ids <- gsub("\\..*", "", rownames(res))
annot <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
               filters = "ensembl_gene_id", values = gene_ids, mart = mart)
res_tbl <- as.data.frame(res); res_tbl$ensembl_gene_id <- gsub("\\..*", "", rownames(res_tbl))
res_annot <- merge(res_tbl, annot, by = "ensembl_gene_id", all.x = TRUE)
res_annot <- res_annot[order(res_annot$padj), ]
readr::write_tsv(res_annot, file.path(out_dir, "deseq2_results_annotated.tsv"))

# DEGs
lfc_thr <- 0.5; padj_thr <- 0.05
DEGs <- subset(res_annot, !is.na(padj) & padj < padj_thr & abs(log2FoldChange) > lfc_thr)
readr::write_tsv(DEGs, file.path(out_dir, sprintf("DEGs_MS_vs_Control_padj%.2f_LFC%.2f.tsv", padj_thr, lfc_thr)))

# Volcano
df <- res_annot
if (!"external_gene_name" %in% names(df)) df$external_gene_name <- NA_character_
df$status <- "NotSig"
df$status[!is.na(df$padj) & df$padj < padj_thr & df$log2FoldChange >  lfc_thr] <- "Up"
df$status[!is.na(df$padj) & df$padj < padj_thr & df$log2FoldChange < -lfc_thr] <- "Down"
df$mlog10padj <- -log10(df$padj); df$mlog10padj[!is.finite(df$mlog10padj)] <- NA
sig <- df[!is.na(df$padj) & df$padj < padj_thr & abs(df$log2FoldChange) > lfc_thr, ]
top10 <- head(sig[order(sig$padj, -abs(sig$log2FoldChange)), ], 10)

p <- ggplot(df, aes(x = log2FoldChange, y = mlog10padj, color = status)) +
  geom_point(size = 1.3, alpha = 0.8, na.rm = TRUE) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
  ggrepel::geom_text_repel(data = top10, aes(label = external_gene_name), size = 3) +
  scale_color_manual(values = c(NotSig = "grey70", Up = "red", Down = "blue")) +
  labs(x = "log2 fold change", y = expression(-log[10]("adjusted p-value")),
       color = "Status", title = "MS vs Control — Volcano") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "volcano_DEGs.png"), p, width = 7, height = 5, dpi = 300)

# GSEA (Hallmark)
symbols <- res_annot$external_gene_name
entrez  <- mapIds(org.Hs.eg.db, keys = symbols, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
res_annot$ENTREZID <- entrez
score <- if (!all(is.na(res_annot$stat))) res_annot$stat else res_annot$log2FoldChange
rank_df <- tibble::tibble(ENTREZID = res_annot$ENTREZID, score = score) |>
  dplyr::filter(!is.na(ENTREZID), is.finite(score)) |>
  dplyr::group_by(ENTREZID) |>
  dplyr::summarise(score = score[which.max(abs(score))], .groups = "drop")
ranks <- sort(setNames(rank_df$score, rank_df$ENTREZID), decreasing = TRUE)

m_h <- msigdbr(species = "Homo sapiens", category = "H") |>
  dplyr::select(gs_name, entrez_gene)

set.seed(42)
gseaH <- GSEA(ranks, TERM2GENE = m_h, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, verbose = FALSE)
readr::write_tsv(as.data.frame(gseaH@result), file.path(out_dir, "GSEA_Hallmark_results.tsv"))

png(file.path(out_dir, "GSEA_Hallmark_top4.png"), width = 900, height = 700)
print(enrichplot::gseaplot2(gseaH, geneSetID = 1:4, title = "Top Hallmark pathways"))
dev.off()

gseaH2 <- pairwise_termsim(gseaH)
png(file.path(out_dir, "GSEA_Hallmark_emap.png"), width = 1200, height = 900)
print(emapplot(gseaH2, showCategory = 10))
dev.off()
```
Outputs ($VSC_DATA/Bioinfo_course/MS_microglia_DEA):
deseq2_results_annotated.tsv, DEGs_MS_vs_Control_padj0.05_LFC0.50.tsv, volcano_DEGs.png,
GSEA_Hallmark_results.tsv, GSEA_Hallmark_top4.png, GSEA_Hallmark_emap.png

  
	 







  