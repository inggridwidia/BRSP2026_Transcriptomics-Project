# Gene Expression Profiling of the Neurovascular Unit in Type 2 Diabetes Mellitus


# 1. Introduction
Type 2 Diabetes Mellitus (T2DM) is a major risk factor for neurodegenerative diseases, including dementia and Alzheimer’s disease. This damage is thought to arise from peripheral insulin resistance, which contributes to dysfunction within the brain’s neurovascular unit. This project aims to characterize synchronous gene expression changes in cortical neurons and endothelial cells to identify shared mechanisms of damage associated with chronic hyperglycemia.

# 2. Methods
This analysis is based on the publicly available dataset GSE161355 retrieved from the NCBI Gene Expression Omnibus (GEO), utilizing the Affymetrix Human Genome U133 Plus 2.0 Array platform. 

The analysis was conducted in R following a systematic workflow comprising these primary stages:
### 2.1. Data Preprocessing 
Microarray data normalization and Log2 transformation.

```R
# Data Acquisition and Normalization
gset <- getGEO("GSE161355", GSEMatrix = TRUE)[[1]]
ex <- exprs(gset)
ex[ex <= 0] <- NA
ex <- log2(ex) # Log2 transformation for normal distribution
```

### 2.2. Sample Filtering by Cell Type
The GSE161355 dataset contains various cell types (endothelial cells, neurons, and astrocytes). Filtering was performed to retain endothelial and neuronal populations, enabling a focused comparison between diabetic and control conditions.

```R
# Selecting specific groups (Focusing on Endothelial and Neuronal cells)
grup_focus <- c(
    "Control.Endothelial.cells", 
    "Diabetic.Endothelial.cells", 
    "Control.Neurones", 
    "Diabetic.Neurones"
)

# Filter ExpressionSet and reset factor levels
gset_filtered <- gset[, gset$group %in% grup_fokus]
gset_filtered$group <- factor(gset_filtered$group, levels = grup_focus)
```

### 2.3. Differential Expression Analysis (`limma`)
Differentially expressed genes (DEGs) were identified using linear modeling implemented in the `limma` package. Given the limited sample size (n = 6), a significance threshold of raw p-value < 0.05 was applied as an exploratory cutoff to increase sensitivity for pathway discovery, while acknowledging the potential for an elevated false-positive rate.

```R
# Statistical Modeling with limma
design <- model.matrix(~0 + gset_filtered$group)
fit <- lmFit(gset_filtered, design)
contrast_matrix <- makeContrasts(
    Diabetes_Efek_Endotel = Diabetic.Endothelial.cells - Control.Endothelial.cells,
    Diabetes_Efek_Neuron = Diabetic.Neurones - Control.Neurones,
    levels = design
)
fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))

# Extracting final results for each cell type
res_endotel <- topTable(fit2, coef="Diabetes_Efek_Endotel", adjust="fdr", number=Inf)
res_neuron <- topTable(fit2, coef="Diabetes_Efek_Neuron", adjust="fdr", number=Inf)

# Filtering using Raw P-Value < 0.05 (Exploratory Threshold)
res_endotel_sig <- res_endotel[res_endotel$P.Value < 0.05, ]
res_neuron_sig <- res_neuron[res_neuron$P.Value < 0.05, ]
```

### 2.4. _Genes Annotation_
Probe IDs were mapped to gene symbols using the `hgu133plus2.db` database to facilitate biological interpretation.

```R
# Gene Mapping and Annotation
gene_ann <- AnnotationDbi::select(
    hgu133plus2.db, 
    keys = probe_ids, 
    columns = c("SYMBOL", "GENENAME"), 
    keytype = "PROBEID"
)
```

### 2.5. Data Visualization (Volcano Plot, Venn Diagram, & Heatmap)
Data visualization was performed to assess the distribution of gene expression values and to evaluate global and sample-level clustering patterns.

***1. Volcano Plot***
```R
# Volcano Plot: Endotel (Threshold: P < 0.05 & |logFC| > 0.5)
res_endotel$status <- "NO"
res_endotel$status[res_endotel$logFC > 0.5 & res_endotel$P.Value < 0.05] <- "UP"
res_endotel$status[res_endotel$logFC < -0.5 & res_endotel$P.Value < 0.05] <- "DOWN"
ggplot(res_endotel, aes(x = logFC, y = -log10(P.Value), color = status)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  theme_minimal() + labs(title = "Volcano Plot: Endothelial Cells"
)

# Volcano Plot: Neuron (Threshold: P < 0.05 & |logFC| > 0.5)
res_neuron$status <- "NO"
res_neuron$status[res_neuron$logFC > 0.5 & res_neuron$P.Value < 0.05] <- "UP"
res_neuron$status[res_neuron$logFC < -0.5 & res_neuron$P.Value < 0.05] <- "DOWN"
ggplot(res_neuron, aes(x = logFC, y = -log10(P.Value), color = status)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  theme_minimal() + labs(title = "Volcano Plot: Neuronal Cells"
)
```
***2. Venn Diagram***
```R
# Venn Diagram: Comparing overlapping DEGs across cell types
draw.pairwise.venn(
  area1 = nrow(res_endotel_sig),
  area2 = nrow(res_neuron_sig),
  cross.area = length(common_genes),
  category = c("Endotel", "Neuron"),
  fill = c("skyblue", "pink")
)
```
***3. Heatmap***
```R
# Heatmap for Common DEGs
common_genes <- intersect(
  res_endotel$Gene.symbol[res_endotel$P.Value < 0.05 & abs(res_endotel$logFC) > 0.5],
  res_neuron$Gene.symbol[res_neuron$P.Value < 0.05 & abs(res_neuron$logFC) > 0.5]
)

# Expression matrix
common_genes <- common_genes[common_genes != "" & !is.na(common_genes)]
mat_heatmap <- ex_filtered[rownames(res_endotel[res_endotel$Gene.symbol %in% common_genes, ]), ]
rownames(mat_heatmap) <- res_endotel$Gene.symbol[res_endotel$Gene.symbol %in% common_genes]
mat_heatmap <- mat_heatmap[!duplicated(rownames(mat_heatmap)), ] # Unique Gene Representation

# Visualization with pheatmap
pheatmap(
  mat_heatmap, 
  scale = "row", 
  clustering_method = "ward.D2", 
  annotation_col = data.frame(Group = gset_filtered$group),
  color = colorRampPalette(c("blue", "white", "red"))(100)
)
```

### 2.6. _Enrichment Analysis_ 
Common DEGs between the two cell types were identified via intersection analysis, followed by functional enrichment using the Gene Ontology and KEGG databases to identify enriched pathways.

```R
# GO Enrichment (Biological Process)
go_bp <- enrichGO(gene = genes_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP")

# KEGG Pathway Mapping & Visualization
kegg <- enrichKEGG(gene = genes_entrez$ENTREZID, organism = 'hsa')
pathview(gene.data = kegg_logFC, pathway.id = "hsa04216", species = "hsa")
```

The complete R script used for this analysis (including `limma`, `ggplot2` visualization, and `clusterProfiler` enrichment) is available in the main directory of this repository (file `Coding GSE161355.R`)


# 3. Results and Discussion
### 3.1. Identification of Common DEGs

Integrative analysis using Venn diagrams successfully identified 77 overlapping genes (common DEGs) that were significantly dysregulated in both endothelial and neuronal cells.

![Venn_Diagram.png](plot-result/Venn_Diagram.png)

Figure 1. Venn diagram illustrating significant gene overlaps endothelial and neuronal cells.

![Heatmap.png](plot-result/Heatmap.png)

Figure 2. Heatmap of the top 77 common DEGs showing distinct clustering between Normal and Diabetic groups.

### 3.2. Functional Pathway Analysis (GO & KEGG)

![GO_Analysis.png](plot-result/GO_Analysis.png)

Figure 3. Gene Ontology (GO) enrichment analysis of the 77 common DEGs revealed significant involvement in biological processes related to cell adhesion, including positive regulation of cell–cell adhesion and cell–matrix adhesion.

![KEGG_Enrich.png](plot-result/KEGG_Enrich.png)

Figure 4. KEGG pathway analysis identified seven enriched pathways (p-value < 0.2).

![ferroptosis_hsa04216.pathview.png](plot-result/ferroptosis_hsa04216.pathview.png)

Figure 5. _Ferroptosis_ (hsa04216)

![focal-adhesion_hsa04510.pathview.png](plot-result/focal-adhesion_hsa04510.pathview.png)

Figure 6. _Focal Adhesion_ (hsa04510)

![Adherens-Junction_hsa04520.pathview.png](plot-result/Adherens-Junction_hsa04520.pathview.png)

Figure 7. _Adherens Junction_ (hsa04520)

![Integrin-signaling_hsa04518.pathview.png](plot-result/Integrin-signaling_hsa04518.pathview.png)

Figure 8. _Integrin signaling_ (hsa04518)

# 4. Conclusion
Neurovascular unit impairment in Type 2 Diabetes reflects tightly coupled interactions between vascular dysfunction and neuronal degeneration, suggesting that effective therapeutic strategies should simultaneously target both compartments.
