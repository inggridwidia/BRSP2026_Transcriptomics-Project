# PROJECT: Gene Expression Analysis of Intestinal Adaptation Following Gastric Bypass in Type 2 Diabetes Patients
# DATASET: GSE281144 
# PLATFORM: Microarray (Affymetrix Human Transcriptome Array 2.0)
# DESCRIPTION: Differential Expression Genes (DEGs) analysis comparing 
# Pre-Op vs Post-Op conditions using R/Bioconductor.


# ------------------------------------------------------------------------------
# PART A. ENVIRONMENT SETUP
# Ensure required bioinformatics libraries are installed and available
# ------------------------------------------------------------------------------

# Install BiocManager to access the Bioconductor repository 
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install core packages (GEOquery for data retrieval, limma for differential expression analysis) 
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE)

# Install CRAN packages for visualization, data manipulation, and dimensionality reduction 
install.packages(c("pheatmap", "ggplot2", "dplyr", "umap"))

# Load required libraries
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(AnnotationDbi)
library(umap)



# ------------------------------------------------------------------------------
# PART B. DATA ACQUISITION FROM GEO
# This step retrieves the 'ExpressionSet' containing both raw data and metadata.
# ------------------------------------------------------------------------------

# Download dataset directly using the GSE281144 accession ID
gset <- getGEO("GSE281144", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

# ExpressionSet components:
# - exprs() : Gene expression matrix
# - pData() : Sample metadata
# - fData() : Feature metadata (probes / genes)



# ------------------------------------------------------------------------------
# PART C. DATA PRE-PROCESSING & NORMALIZATION
# Raw microarray data often exhibits heteroscedasticity.
# Log2 transformation is performed to normalize the data distribution.
# ------------------------------------------------------------------------------

# exprs(): Extract gene expression matrix
# Rows    = Probes/Genes
# Columns = Samples
ex <- exprs(gset) 

# Determine if the dataset requires Log transformation
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

#IF statement:
# If LogTransform = TRUE, apply log2 transformation.
# Values <= 0 are replaced with NA to avoid undefined log results.
if (LogTransform) {
  ex[ex <= 0] <- NA 
  ex <- log2(ex)    
}



# ------------------------------------------------------------------------------
# PART D. DATA FILTERING & GROUP DEFINITION (SUBSETTING)
# Subset dataset to retain diabetic samples only
# Define Pre-Op vs Post-Op conditions for comparison
# ------------------------------------------------------------------------------

# Identify the specific metadata column containing biological information
View(pData(gset))

# Identify samples with diabetic status in metadata (characteristics_ch1.1)
is_diabetic <- grepl("diabetic", pData(gset)$characteristics_ch1.1, ignore.case = TRUE)

# Filter the gset object to include only diabetic samples
# Update 'gset' and the expression matrix 'ex' with the filtered diabetic samples
gset <- gset[, is_diabetic]
ex <- exprs(gset)

# Extract time-point information (Pre-Op vs. Post-Op) from sample titles
# baseline = Pre-Op, other time points (1 month/6 months) = Post-Op
titles <- pData(gset)$title
group_info <- ifelse(grepl("baseline", titles, ignore.case = TRUE), "PreOp_DM", "PostOp_DM")

# Standardize group names and convert into factors for statistical modeling
groups <- make.names(group_info)
gset$group <- factor(groups)
nama_grup <- levels(gset$group)

# Define contrast: Post-Op (Treatment) vs. Pre-Op (Control) comparison
# Verify group ordering for contrast definition
print(nama_grup)

grup_post <- nama_grup[1] 
grup_pre  <- nama_grup[2]



# ------------------------------------------------------------------------------
# PART E. DIFFERENTIAL EXPRESSION STATISTICAL ANALYSIS (LIMMA)
# Using Linear Models (limma) to identify genes 
# significant expression changes between groups
# ------------------------------------------------------------------------------

# Construct design matrix representing the experimental group structure
design <- model.matrix(~0 + gset$group)
colnames(design) <- levels(gset$group)

# Fit linear model to estimate expression differences across groups
fit <- lmFit(ex, design)

# Create Contrast: Define the specific comparison (Post-Op minus Pre-Op)
contrast_formula <- paste(grup_post, "-", grup_pre, sep = "")
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

# Apply empirical Bayes (eBayes) moderation to stabilize variance estimates
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract differentially expressed genes with FDR-adjusted significance
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.05
)
head(topTableResults)



# ------------------------------------------------------------------------------
# PART F. BIOLOGICAL ANNOTATION (GENE MAPPING)
# Map probe IDs to gene symbols and gene names for biological interpretation.
# Use feature metadata from GEO (fData) for gene annotation.
# ------------------------------------------------------------------------------

# Retrieve feature metadata (gene info) directly from the gset object
feature_data <- fData(gset)

# Verify the metadata column containing biological information (Gene Symbols)
colnames(fData(gset))

# Extract gene information from 'gene_assignment' field
# Format contains accession, gene symbol, and gene description 
# Extract gene symbols (second field in annotation string)
extracted_symbols <- sub("^[^//]* // ([^//]*) // .*$", "\\1", feature_data$gene_assignment)

# Extract gene names (third field in annotation string)
extracted_names <- sub("^[^//]* // [^//]* // ([^//]*) // .*$", "\\1", feature_data$gene_assignment)

# Construct annotation mapping (Probe -> Symbol & Name)
anno_map <- data.frame(
  PROBEID = feature_data$ID,
  SYMBOL = extracted_symbols,
  GENENAME = extracted_names,
  stringsAsFactors = FALSE
)

# Merge differential expression results with gene annotation
topTableResults <- merge(topTableResults, anno_map, 
                         by.x = "row.names", 
                         by.y = "PROBEID", 
                         all.x = TRUE)
colnames(topTableResults)[1] <- "PROBEID"

# Clean missing annotation values
topTableResults$SYMBOL[topTableResults$SYMBOL == "---"] <- NA
topTableResults$GENENAME[topTableResults$GENENAME == "---"] <- NA

# Preview annotated results
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])



# ------------------------------------------------------------------------------
# PART G. QUALITY CONTROL & DATA EXPLORATION (INITIAL VISUALIZATION)
# Verify data quality and assess sample clustering patterns prior to downstream analysis
# ------------------------------------------------------------------------------

# 1. BOXPLOT: Assess expression distribution across samples
# Identify potential outliers and batch effects
# Assign colors based on sample groups
group_colors <- as.numeric(gset$group) 

boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Expression Value Distribution: Post-Op vs. Pre-Op",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.7
)

# 2. DENSITY PLOT: Visualize overall expression distribution across samples
# Confirm that log2 normalization effectively aligned the distributions across groups.
# Reshape expression data and group info into a long-format data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Global Expression Density Plot",
    x = "Log2 Expression",
    y = "Density"
  )

# 3. UMAP: Low-Dimensional Visualization
# Evaluate sample clustering based on global gene expression patterns
# Transpose expression matrix for sample-level analysis
umap_input <- t(ex) 

# Perform UMAP dimensionality reduction
umap_result <- umap(umap_input)

# Store coordinates in a data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

# Generate UMAP visualization
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP: Pengelompokan Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )



# ------------------------------------------------------------------------------
# PART H. DATA VISUALIZATION (RESULTS OVERVIEW)
# ------------------------------------------------------------------------------

# 1. VOLCANO PLOT: Visualize global differential expression patterns (Log Fold Change vs. Significance) 
volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

# Classify gene status based on differential expression thresholds
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.05] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.05] <- "DOWN"

# Generate volcano plot visualization
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Differential Gene Expression: Post-Op vs Pre-Op", 
       subtitle = "Dataset GSE281144 - Gastric Bypass Adaptation")

# 2. HEATMAP: Visualize expression patterns of the top 50 differentially expressed genes
# Select top 50 genes ranked by adjusted p-value
# Remove invalid gene annotations (NA or placeholder values)
topTableAnnotated <- subset(topTableResults, !is.na(SYMBOL) & SYMBOL != "" & SYMBOL != "---")
topTableAnnotated <- topTableAnnotated[order(topTableAnnotated$adj.P.Val), ]
top50 <- head(topTableAnnotated, 50)

# Extract expression matrix for the selected genes using PROBEID
mat_heatmap <- ex[top50$PROBEID, ]

# Assign gene symbols as row names
gene_label <- top50$SYMBOL
rownames(mat_heatmap) <- gene_label

# Remove missing values and zero-variance genes to ensure valid clustering
mat_heatmap <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, ]
mat_heatmap <- mat_heatmap[apply(mat_heatmap, 1, var) > 0, ]

# Define column annotations (sample groups)
annotation_col <- data.frame(
  Group = gset$group
)
rownames(annotation_col) <- colnames(mat_heatmap)

# Generate clustered heatmap
pheatmap(
  mat_heatmap,
  scale = "row",                
  annotation_col = annotation_col,
  show_colnames = FALSE,       
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Heatmap: Top 50 Identified Genes (Diabetic Post-Op vs Pre-Op)"
)



# ------------------------------------------------------------------------------
# PART I. FINAL DATA EXPORT
# Save the analysis results to CSV format for further reporting.
# ------------------------------------------------------------------------------
# Extract significant genes for downstream GO/KEGG analysis
sig_genes <- subset(topTableAnnotated, adj.P.Val < 0.05)$SYMBOL

# Export gene list as text file
write.table(sig_genes, "daftar_gen_signifikan.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Save full DEG results table as CSV
write.csv(topTableResults, "GSE281144_DEG_Results.csv")
message("Analysis complete: results exported successfully")



# ------------------------------------------------------------------------------
# PART J. GENE ONTOLOGY (GO) ENRICHMENT ANALYSIS
# Focus on Gene Ontology Biological Processes (BP)
# BP analysis provides insight into functional alterations associated with Post-Op conditions
# ------------------------------------------------------------------------------

# Install required Bioconductor packages for GO enrichment analysis
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"), update = FALSE)
n

# Load libraries for enrichment analysis and visualization
library(clusterProfiler)
library(org.Hs.eg.db) # Human genome annotation database (Homo sapiens)
library(enrichplot)
library(ggplot2)

# Filter significant genes with annotated gene symbols
sig_genes_df <- subset(topTableResults, adj.P.Val < 0.05 & !is.na(SYMBOL) & SYMBOL != "" & SYMBOL != "---")

# Convert gene symbols to Entrez IDs
# Ensure compatibility with GO and KEGG databases
genes_entrez <- bitr(
  sig_genes_df$SYMBOL, 
  fromType = "SYMBOL", 
  toType   = "ENTREZID", 
  OrgDb    = org.Hs.eg.db
)
cat("Total unique genes included in GO analysis:", nrow(genes_entrez), "\n")

# Perform GO enrichment analysis (Biological Process)
# Identify enriched biological processes associated with the gene set
go_results <- enrichGO(
  gene          = genes_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE  # Convert IDs back to Symbols for plot readability
)

# Generate dotplot visualization
# Display top 10 enriched biological processes
dotplot(go_results, showCategory = 10) + 
  labs(title = "GO Enrichment: Biological Processes",
       subtitle = "Overall functional changes after Gastric Bypass",
       x = "Gene Ratio",
       y = "Biological Process")



# ------------------------------------------------------------------------------
# PART K. KEGG PATHWAY ANALYSIS
# KEGG enrichment analysis to identify relevant metabolic and signaling pathways
# ------------------------------------------------------------------------------

# Perform KEGG enrichment analysis
# Use previously generated Entrez ID list
kegg_enrich <- enrichKEGG(
  gene         = genes_entrez$ENTREZID,
  organism     = 'hsa', 
  pvalueCutoff = 0.05
)

# Convert KEGG results to data frame
kegg_table <- as.data.frame(kegg_enrich)
View(kegg_table) # View results interactively (RStudio only)

# Generate KEGG enrichment barplot
# Display top enriched pathways
library(ggplot2)
barplot(kegg_enrich, showCategory = 18) +
  scale_y_discrete(
    labels = function(x) paste0(kegg_enrich@result$ID[match(x, kegg_enrich@result$Description)], ": ", x)
  ) +
  labs(
    title    = "KEGG Pathway Enrichment Analysis",
    subtitle = "All 18 significantly affected pathways (p.adj < 0.05)",
    x        = "Gene Count",
    y        = "Pathway Description (ID: Name)"
  ) +
  theme_minimal()

# Generate pathway visualization outputs (PNG format)
# Red = Up-regulated, Green = Down-regulated.
# Install and load pathview package if needed
BiocManager::install("pathview")
library(pathview)

# Prepare logFC data indexed by Entrez IDs
kegg_logFC <- sig_genes_df$logFC
names(kegg_logFC) <- genes_entrez$ENTREZID[match(sig_genes_df$SYMBOL, genes_entrez$SYMBOL)]

# Pathview: Vitamin digestion and absorption (hsa04977)
pathview(
  gene.data  = kegg_logFC,
  pathway.id = "hsa04977", 
  species    = "hsa",
  limit      = list(gene = 2, cpd = 1)
)

# Pathview: PI3K-Akt signaling (hsa04151)
pathview(
  gene.data  = kegg_logFC,
  pathway.id = "hsa04151", 
  species    = "hsa",
  limit      = list(gene = 2, cpd = 1)
)

# Pathview: PPAR signaling (hsa03320)
pathview(
  gene.data  = kegg_logFC,
  pathway.id = "hsa03320", 
  species    = "hsa",
  limit      = list(gene = 2, cpd = 1)
)
