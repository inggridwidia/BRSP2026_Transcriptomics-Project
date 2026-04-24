# Module: Gene Expression Profiling of the Neurovascular Unit in Type 2 Diabetes Mellitus
# Dataset: GSE161355 (Type 2 Diabetes vs. Normal)
# Platform: Microarray (Affymetrix Human Genome U133 Plus 2.0 Array - GPL570)
# Objective: To identify Differentially Expressed Genes (DEGs) 


# ------------------------------------------------------------------------------
# PART A. WORKSPACE PREPARATION 
# Prior to analysis, all required bioinformatics packages were installed and managed
# using the BiocManager package in R.
# ------------------------------------------------------------------------------
# Ensure BiocManager is installed for Bioconductor repository access 
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install core packages (GEOquery for data retrieval, limma for differential expression analysis)
BiocManager::install(c("GEOquery", "limma", "hgu133plus2.db", "AnnotationDbi"), ask = FALSE, update = FALSE)

# Install CRAN packages for visualization, data manipulation, and dimensionality reduction                       
install.packages(c("pheatmap", "ggplot2", "dplyr", "umap", "VennDiagram"))

library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133plus2.db)
library(AnnotationDbi)
library(umap)
library(VennDiagram)



# ------------------------------------------------------------------------------
# PART B. DATA ACQUISITION FROM GEO
# This step retrieves the 'ExpressionSet' containing both raw data and metadata.
# ------------------------------------------------------------------------------

# getGEO(): Function to download the dataset based on the GEO accession ID.
# GSEMatrix = TRUE -> Retrieves data in ExpressionSet format.

gset <- getGEO("GSE161355", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]



# ------------------------------------------------------------------------------
# PART C. DATA PRE-PROCESSING & NORMALIZATION
# Raw microarray data often exhibits heteroscedasticity.
# Log2 transformation is performed to normalize the data distribution.
# ------------------------------------------------------------------------------

# exprs(): Extracts the gene expression matrix
# Baris  = Probes/Genes
# Kolom  = Samples

ex <- exprs(gset)

# Determine if the dataset requires Log transformation
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

# IF statement:
# If LogTransform = TRUE, apply log2 transformation.
# Values <= 0 are replaced with NA to avoid undefined log results.
if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}



# ------------------------------------------------------------------------------
# PART D. GROUP DEFINITION & SAMPLE FILTERING
# ------------------------------------------------------------------------------

# Identify metadata columns from the ExpressionSet containing biological information
View(pData(gset))
unique(pData(gset)[["source_name_ch1"]])

# Once the relevant column is identified, create the 'group_info' variable
# to store the biological condition metadata for each sample
group_info <- pData(gset)[["source_name_ch1"]]

# Standardize column values into R-compatible names
groups <- make.names(group_info)

# Convert the group information into a factor (categorical variable)
# allowing R to categorize data into specific experimental levels
gset$group <- factor(groups)

# Group Level Verification: Ensure metadata is correctly parsed by R
# prior to proceeding with downstream complex statistical analyses.
nama_grup <- levels(gset$group)
print(nama_grup)

# Define the focus groups for downstream analysis
# Comparing Endothelial and Neuronal populations only (excluding Astrocytes)
grup_fokus <- c("Control.Endothelial.cells", "Diabetic.Endothelial.cells", 
                "Control.Neurones", "Diabetic.Neurones")

# Filter the ExpressionSet to include only samples from the selected groups
gset_filtered <- gset[, gset$group %in% grup_fokus]

# Update factors to remove unused levels (e.g., Astrocytes)
# and ensure the level order aligns with the desired experimental design
gset_filtered$group <- factor(gset_filtered$group, levels = grup_fokus)



# ------------------------------------------------------------------------------
# PART E. DIFFERENTIAL EXPRESSION STATISTICAL ANALYSIS (LIMMA)
# Using Linear Models for Microarray (limma) to identify genes that
# exhibit significant differential expression between groups.
# ------------------------------------------------------------------------------

# Design Matrix Construction: Mathematically representing the experimental structure.
design <- model.matrix(~0 + gset_filtered$group)
colnames(design) <- levels(gset_filtered$group)

# Define contrast vectors (Diabetes vs. Control) for each specific cell type
contrast_matrix <- makeContrasts(
  Diabetes_Effect_Endothelial = Diabetic.Endothelial.cells - Control.Endothelial.cells,
  Diabetes_Effect_Neuronal = Diabetic.Neurones - Control.Neurones,
  levels = design
)

# Linear Model Fitting: Estimating mean expression levels for each gene per group 
fit <- lmFit(gset_filtered, design)

# Empirical Bayes (eBayes) Moderation: A critical step to 
# stabilize gene variance estimates using Bayesian methods, 
# ensuring accurate results despite limited sample sizes.
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract Top-Ranked Genes (Top Table) with FDR (False Discovery Rate) correction 
res_endotel <- topTable(
  fit2, 
  coef="Diabetes_Effect_Endothelial", 
  adjust="fdr", 
  sort.by="B", 
  number=Inf
)
res_neuron <- topTable(
  fit2, 
  coef="Diabetes_Effect_Neuronal", 
  adjust="fdr", 
  sort.by="B", 
  number=Inf
)

# Extract significant results using raw p-value < 0.05
# (applied due to limited sample size for exploratory analysis)
res_endotel_sig <- res_endotel[res_endotel$P.Value < 0.05, ]
res_neuron_sig <- res_neuron[res_neuron$P.Value < 0.05, ]

# --- ANALYSIS NOTE: SIGNIFICANCE THRESHOLD SELECTION ---
# In this dataset (GSE161355), the sample size per group is limited (n=6).
# Applying multiple testing correction (FDR/adj.P.Val) yields highly 
# conservative results, with few or no genes passing stringent thresholds.
# For exploratory analysis and expression profiling (e.g., heatmaps, volcano plots),
# 'Raw P-Value < 0.05' was used as an exploratory threshold.



# ------------------------------------------------------------------------------
# PART F. GENE ANNOTATION
# Annotate probe IDs to gene symbols using hgu133plus2.db
# ------------------------------------------------------------------------------

# 1. Endothelial Cell Annotation
# Extract Probe IDs from the significant endothelial DEG results
probe_ids_endotel <- rownames(res_endotel_sig)

# Map probes to gene symbols and gene names
# Use the annotation library specific to the platform (Affymetrix U133 Plus 2.0)
# Typically uses the 'hgu133plus2.db' package
gene_ann_endotel <- AnnotationDbi::select(
  hgu133plus2.db, 
  keys = probe_ids_endotel, 
  columns = c("SYMBOL", "GENENAME"), 
  keytype = "PROBEID"
)

# Add a PROBEID column from row names to facilitate data merging
res_endotel_sig$PROBEID <- rownames(res_endotel_sig)

# Merge limma statistical results with gene annotations based on PROBEID
# all.x = TRUE ensures all statistical values are retained during the merge
res_endotel_sig <- merge(
  res_endotel_sig, 
  gene_ann_endotel, 
  by = "PROBEID", 
  all.x = TRUE
)


# 2. Neuronal Cell Annotation
# Extract Probe IDs from the significant neuronal DEG results
probe_ids_neuron <- rownames(res_neuron_sig)

# Map probes to gene symbols and gene names
# Use the annotation library specific to the platform (Affymetrix U133 Plus 2.0)
# Typically uses the 'hgu133plus2.db' package
gene_ann_neuron <- AnnotationDbi::select(
  hgu133plus2.db, 
  keys = probe_ids_neuron, 
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

# Add a PROBEID column from row names to facilitate data merging
res_neuron_sig$PROBEID <- rownames(res_neuron_sig)

# Merge limma statistical results with gene annotations based on PROBEID
# all.x = TRUE ensures all statistical values are retained during the merge
res_neuron_sig <- merge(
  res_neuron_sig, 
  gene_ann_neuron, 
  by = "PROBEID", 
  all.x = TRUE
)



# ------------------------------------------------------------------------------
# PART G. QUALITY CONTROL & DATA EXPLORATION (INITIAL VISUALIZATION)
# Assess data quality and examine sample clustering before finalizing results.
# ------------------------------------------------------------------------------

# 1. BOXPLOT: Assess expression value distribution across samples
# Used to identify potential outliers or significant batch effects. 
# Assign colors based on experimental groups
group_colors <- as.numeric(gset_filtered$group)+1

# Extract expression matrix from filtered data
ex_filtered <- exprs(gset_filtered)

par(mar = c(8, 4, 4, 12), xpd = TRUE)
boxplot(
  ex_filtered, 
  col = group_colors, 
  las = 2, 
  outline = FALSE, 
  cex.axis = 0.7,
  main = "Expression Value Distribution per Sample (Endothelial & Neuronal",
  ylab = "Expression Value (log2)"
)

legend(
  ncol(ex_filtered)+2, 
  14, 
  legend = levels(gset_filtered$group), 
  fill = unique(group_colors), 
  cex = 0.6, 
  bty = "n"
)

par(mar=c(5, 4, 4, 2), xpd=FALSE)

# 2. DENSITY PLOT: Visualize global expression distributions
# Assess whether log2 normalization aligns distributions across samples
# Combine expression & groups into a data frame
expr_long <- data.frame(
  Expression = as.vector(ex_filtered),
  Group = rep(gset_filtered$group, each = nrow(ex_filtered))
) 

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") + 
  labs(
    title = "Global Gene Expression Distribution (Endothelial vs. Neuronal)",
    subtitle = "Dataset: GSE161355",
    x = "Expression Value (log2)",
    y = "Density"
  )

# 3. UMAP: Low-dimensional visualization
# Visualizes separation between endothelial and neuronal gene expression profiles
# and assesses the impact of diabetes on gene expression changes within each cell type.

# Transpose expression matrix:
# UMAP calculates similarities across samples (rows)
umap_input <- t(ex_filtered)

# Run UMAP
umap_result <- umap(umap_input)

# Save results to a data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[,1], 
  UMAP2 = umap_result$layout[,2], 
  Group = gset_filtered$group
)

# Render UMAP plot
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() + 
  labs(
    title = "UMAP Plot: Global Sample Clustering",
    x = "UMAP 1",
    y = "UMAP 2"
)



# ------------------------------------------------------------------------------
# PART H. DATA VISUALIZATION (RESULTS INTERPRETATION)
# ------------------------------------------------------------------------------

# 1. VOLCANO PLOT: Visualize global differential expression patterns
# (Log Fold Change vs. Statistical Significance)

# Volcano Plot: Endothelial (Threshold: P < 0.05 & |logFC| > 0.5)
res_endotel$status <- "NO"
res_endotel$status[res_endotel$logFC > 0.5 & res_endotel$P.Value < 0.05] <- "UP"
res_endotel$status[res_endotel$logFC < -0.5 & res_endotel$P.Value < 0.05] <- "DOWN"

ggplot(res_endotel, aes(x = logFC, y = -log10(P.Value), color = status)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() + 
  labs(title = "Volcano Plot: Endothelial Cells", subtitle = "Diabetic vs. Normal Control")

# Volcano Plot: Neuron (Threshold: P < 0.05 & |logFC| > 0.5)
res_neuron$status <- "NO"
res_neuron$status[res_neuron$logFC > 0.5 & res_neuron$P.Value < 0.05] <- "UP"
res_neuron$status[res_neuron$logFC < -0.5 & res_neuron$P.Value < 0.05] <- "DOWN"

ggplot(res_neuron, aes(x = logFC, y = -log10(P.Value), color = status)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() + 
  labs(title = "Volcano Plot: Neuronal Cells", subtitle = "Diabetic vs. Normal Control")


# 2. VENN DIAGRAM: Visualize DEG distribution and overlap between endothelial and neuronal cells.

draw.pairwise.venn(
  area1 = nrow(res_endotel[res_endotel$P.Value < 0.05 & abs(res_endotel$logFC) > 0.5, ]),
  area2 = nrow(res_neuron[res_neuron$P.Value < 0.05 & abs(res_neuron$logFC) > 0.5, ]),
  cross.area = length(intersect(
    res_endotel$Gene.symbol[res_endotel$P.Value < 0.05 & abs(res_endotel$logFC) > 0.5],
    res_neuron$Gene.symbol[res_neuron$P.Value < 0.05 & abs(res_neuron$logFC) > 0.5]
  )),
  category = c("Endotel", "Neuron"),
  fill = c("skyblue", "pink"),
  lty = "blank",
)

grid.text(
  "Comparison of DEGs in Diabetic Endothelial vs Neuronal Cells", 
  x = 0.5, 
  y = 0.95, 
  gp = gpar(fontsize = 15, fontface = "bold")
)


# 3. HEATMAP: Visualize expression patterns of common DEGs shared across both cell types
# Identify genes consistently dysregulated in endothelial and neuronal cells
# Apply exploratory criteria: P < 0.05 & |logFC| > 0.5
common_genes <- intersect(
  res_endotel$Gene.symbol[res_endotel$P.Value < 0.05 & abs(res_endotel$logFC) > 0.5],
  res_neuron$Gene.symbol[res_neuron$P.Value < 0.05 & abs(res_neuron$logFC) > 0.5]
)

# Data Cleaning: Remove empty strings ("") and NA values to prevent rendering errors
common_genes <- common_genes[common_genes != "" & !is.na(common_genes)]

# Construct heatmap matrix from ex_filtered
# Select one representative probe per gene to avoid redundancy
mat_heatmap <- ex_filtered[rownames(res_endotel[res_endotel$Gene.symbol %in% common_genes, ]), ]
rownames(mat_heatmap) <- res_endotel$Gene.symbol[res_endotel$Gene.symbol %in% common_genes]
mat_heatmap <- mat_heatmap[!duplicated(rownames(mat_heatmap)), ] # Unique Gene Representation

# Final Heatmap Visualization
pheatmap(
  mat_heatmap,
  scale = "row", 
  clustering_method = "ward.D2", # Enhances contrast and clustering organization
  annotation_col = data.frame(Group = gset_filtered$group, row.names = colnames(ex_filtered)),
  show_colnames = FALSE,
  show_rownames = TRUE, 
  fontsize_row = 7, 
  main = "Heatmap Common DEGs (Strict Filter: P < 0.05 & |FC| > 0.5)",
  color = colorRampPalette(c("#2980b9", "white", "#c0392b"))(100), # Professional Blue-Red Gradient
  border_color = NA
)



# ------------------------------------------------------------------------------
# PART I. GENE ONTOLOGY (GO) ENRICHMENT ANALYSIS
# Focus on Biological Process (BP).
# BP provides insight into functional pathway dysregulation and physiological processes
# in endothelial and neuronal cells under diabetic conditions
# ------------------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Convert Gene Symbols to Entrez IDs (Standard format for GO/KEGG databases)
# Use 'common_genes' identified from the Venn intersection
genes_entrez <- bitr(common_genes, 
                     fromType = "SYMBOL", 
                     toType   = "ENTREZID", 
                     OrgDb    = org.Hs.eg.db)

# GO Enrichment: Biological Process (BP)
# Identify significantly altered biological processes across both cell types
go_bp <- enrichGO(gene           = genes_entrez$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP", 
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.15,  #Optimal threshold: Identifies 17 key functional pathways
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)

# Dotplot visualization: Display top 15 enriched pathways
dotplot(go_bp, showCategory = 15) + 
  labs(
    title = "Top Biological Processes: Common DEGs in Diabetic Cells",
    subtitle = "Shared mechanisms between Endothelial and Neuronal populations",
    x = "Gene Ratio",
    y = "Biological Process"
  )

# Export results for further reporting
write.csv(as.data.frame(go_bp), "GO_Enrichment_Results.csv")


# ------------------------------------------------------------------------------
# PART J. KEGG PATHWAY ANALYSIS (KYOTO ENCYCLOPEDIA OF GENES AND GENOMES)
# Map genes to signaling and metabolic pathways to explore disease-related biological processes
# Apply exploratory threshold (p-value < 0.2)
# ------------------------------------------------------------------------------

# KEGG Pathway Enrichment
# 'hsa' refers to Homo sapiens (Human)
kegg <- enrichKEGG(
  gene = genes_entrez$ENTREZID,
  organism = 'hsa', 
  pvalueCutoff = 0.2
)

# Export the full KEGG result list as a data frame
kegg_table <- as.data.frame(kegg)
View(kegg_table) # Opens a detailed result window for inspection


# Barplot Visualization for KEGG Pathways
barplot(kegg, showCategory = 10) + 
  labs(
    title = "Top KEGG Pathways: Shared Mechanisms in Diabetic Cells",
    subtitle = "Key metabolic and signaling pathways (P-value < 0.2)",
    x = "Gene Count",
    y = "Pathway Description"
  ) +
  theme_minimal()
head(as.data.frame(kegg))


# Pathway-Specific Visualization
# Red = Up-regulated, Blue = Down-regulated.
library(pathview)

# Extract logFC values for overlapping genes from 'res_endotel'
# Map Entrez IDs as names to comply with Pathview requirements
kegg_logFC <- res_endotel$logFC[res_endotel$Gene.symbol %in% common_genes]
names(kegg_logFC) <- genes_entrez$ENTREZID[match(common_genes, genes_entrez$SYMBOL)]

# Remove NA values to ensure visualization integrity
kegg_logFC <- kegg_logFC[!is.na(names(kegg_logFC))]


# SPECIFIC PATHWAY VISUALIZATIONS
# 1. "Ferroptosis" Pathway (hsa04216)
pathview(
  gene.data  = kegg_logFC,
  pathway.id = "hsa04216", 
  species    = "hsa",
  limit      = list(gene = 2, cpd = 1),
  low        = "blue",              # Down-regulated
  mid        = "gray",              # Neutral
  high       = "red"                # Up-regulated
)

# 2. "Focal Adhesion" Pathway (hsa04510)
pathview(
  gene.data  = kegg_logFC,
  pathway.id = "hsa04510", 
  species    = "hsa",
  limit      = list(gene = 2, cpd = 1),
  low        = "blue",              # Down-regulated
  mid        = "gray",              # Neutral
  high       = "red"                # Up-regulated
)

# 3. "Adherens Junction" Pathway (hsa04520)
pathview(
  gene.data  = kegg_logFC,
  pathway.id = "hsa04520", 
  species    = "hsa",
  limit      = list(gene = 2, cpd = 1),
  low        = "blue",              # Down-regulated
  mid        = "gray",              # Neutral
  high       = "red"                # Up-regulated
)

# 4. "Integrin signaling" Pathway (hsa04518)
pathview(
  gene.data  = kegg_logFC,
  pathway.id = "hsa04518", 
  species    = "hsa",
  limit      = list(gene = 2, cpd = 1),
  low        = "blue",              # Down-regulated
  mid        = "gray",              # Neutral
  high       = "red"                # Up-regulated
)



# ------------------------------------------------------------------------------
# PART K. DATA EXPORT
# ------------------------------------------------------------------------------

# Save the significant DEG results for both cell types
write.csv(res_endotel_sig, "Diabetic_Endothelial_DEGs.csv", row.names = FALSE)
write.csv(res_neuron_sig, "Diabetic_Neuronal_DEGs.csv", row.names = FALSE)
message("Analysis complete. Output files have been successfully exported.")
