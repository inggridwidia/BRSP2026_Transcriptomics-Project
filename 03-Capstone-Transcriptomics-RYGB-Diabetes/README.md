# Analisis Transkriptomik Pasca _Roux-en-Y Gastric Bypass_ pada Pasien Diabetes Mellitus Tipe 2

# 1. Pendahuluan
Operasi _Roux-en-Y Gastric Bypass_ (RYGB) tidak hanya bertujuan untuk penurunan berat badan, tetapi juga dikenal efektif dalam memperbaiki kontrol glikemik pada pasien Diabetes Mellitus Tipe 2 (DMT2). Proyek ini bertujuan untuk mengidentifikasi profil ekspresi gen yang berubah secara signifikan (_Differentially Expressed Genes_) serta memetakan perubahan jalur metabolisme dan pensinyalan seluler pasca-prosedur.

# 2. Metode
Analisis ini menggunakan dataset publik GSE281144 yang diperoleh dari database NCBI _Gene Expression Omnibus_ (GEO), menggunakan platform _Affymetrix Human Transcriptome Array_ 2.0.

Analisis dilakukan secara sistematis menggunakan bahasa pemrograman R dengan tahapan utama:
### 2.1. _Preprocessing data_ 
Normalisasi data microarray dan transformasi Log2.

```R
# Data Acquisition and Normalization
gset <- getGEO("GSE281144", GSEMatrix = TRUE)[[1]]
ex <- exprs(gset)
ex[ex <= 0] <- NA 
ex <- log2(ex) # Log2 transformation for normal distribution
```

### 2.2. _Differential Analysis_
Identifikasi DEG menggunakan model linier (`limma`) dengan kriteria Adjusted P-Value < 0.05 untuk membandingkan kelompok Post-Op terhadap Pre-Op.

```R
# Statistical Modeling with limma
design <- model.matrix(~0 + gset$group)
fit <- lmFit(ex, design)
contrast_matrix <- makeContrasts(PostOp_DM - PreOp_DM, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2) # Empirical Bayes moderation

# Get top significant genes (FDR < 0.05)
topTableResults <- topTable(fit2, adjust = "fdr", number = Inf, p.value = 0.05)
```

### 2.3. Genes Annotation
Pemetaan Probe ID ke Gene Symbol menggunakan simbol gen resmi yang disertakan langsung oleh GEO

```R
# Mengambil metadata fitur (gen) langsung dari objek gset
feature_data <- fData(gset)
extracted_symbols <- sub("^[^//]* // ([^//]*) // .*$", "\\1", feature_data$gene_assignment)
extracted_names <- sub("^[^//]* // [^//]* // ([^//]*) // .*$", "\\1", feature_data$gene_assignment)

# Mapping probe
anno_map <- data.frame(
  PROBEID = feature_data$ID,
  SYMBOL = extracted_symbols,
  GENENAME = extracted_names,
  stringsAsFactors = FALSE
)

# Gabungkan ke hasil statistik topTable
topTableResults <- merge(topTableResults, anno_map, by = "PROBEID", all.x = TRUE)
topTableResults$SYMBOL[topTableResults$SYMBOL == "---"] <- NA
topTableResults$GENENAME[topTableResults$GENENAME == "---"] <- NA
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])
```

### 2.4. _Data Visualization (Volcano Plot & Heatmap)_
Visualisasi dilakukan untuk memvalidasi distribusi gen secara global dan melihat pola pengelompokan ekspresi gen antar kondisi klinis (Pre-Op vs Post-Op).
**1. _Volcano Plot_**
```R
# Data Preparation (Dataset GSE281144)
volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

# Klasifikasi Status Ekspresi (Threshold: |logFC| > 1 & adj.P.Val < 0.05)
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.05] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.05] <- "DOWN"

# Visualisasi menggunakan ggplot2
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Differential Gene Expression: Post-Op vs Pre-Op", 
       subtitle = "Dataset GSE281144 - Gastric Bypass Adaptation")
```
**2. _Heatmap (Top 50 DEGs)_**
```R
# Seleksi 50 gen paling signifikan (Top 50 DEGs)
topTableAnnotated <- subset(topTableResults, !is.na(SYMBOL) & SYMBOL != "" & SYMBOL != "---")
topTableAnnotated <- topTableAnnotated[order(topTableAnnotated$adj.P.Val), ]
top50 <- head(topTableAnnotated, 50)

# Matriks ekspresi dengan label Gene Symbol
mat_heatmap <- ex[top50$PROBEID, ]
rownames(mat_heatmap) <- top50$SYMBOL

# Visualisasi Heatmap dengan Klastering Hierarkis
pheatmap(
  mat_heatmap,
  scale = "row",                
  annotation_col = data.frame(Group = gset$group, row.names = colnames(mat_heatmap)),
  show_colnames = FALSE,       
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_method = "complete",
  main = "Heatmap: Top 50 Identified Genes (Diabetic Post-Op vs Pre-Op)"
)
```

### 2.5. _Enrichment Analysis_ 
Interpretasi biologis menggunakan database _Gene Ontology_ (GO) dan _Kyoto Encyclopedia of Genes and Genomes_ (KEGG) untuk memetakan perubahan jalur fungsional.

```R
# GO Enrichment (Biological Process)
go_up <- enrichGO(gene = up_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP")
go_down <- enrichGO(gene = down_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP")

# KEGG Pathway Mapping and Visualization
kegg_enrich <- enrichKEGG(gene = all_entrez_df$ENTREZID, organism = 'hsa')
pathview(gene.data = kegg_logFC, pathway.id = "hsa04151", species = "hsa")
```

Skrip R lengkap yang digunakan untuk analisis ini (termasuk pemrosesan `limma`, visualisasi `ggplot2`, dan pengayaan `clusterProfiler`) tersedia di folder utama repositori ini (file `Coding GSE281144.R`)

# 3. Hasil dan Pembahasan
### 3.1. Visualisasi Ekspresi Gen (DEG)

Analisis menunjukkan pergeseran transkriptomik yang masif pasca-operasi, ditandai dengan pemisahan klaster yang jelas antara pasien sebelum dan sesudah operasi.

![Volcano_Plot.png](plot-result/Volcano_Plot.png)

Gambar 1. Volcano plot menampilkan distribusi signifikansi gen. Titik merah (Up-regulated) dan biru (Down-regulated) menunjukkan gen dengan perubahan ekspresi secara signifikan (FDR < 0.05).

![Heatmap_TOP50_DEGs.png](plot-result/Heatmap_TOP50_DEGs.png)

Gambar 2. Heatmap 50 DEG teratas menunjukkan pola ekspresi yang konsisten secara sistemik pada kelompok Post-Op dibandingkan Pre-Op.

### 3.2. Analisis Fungsional dan Jalur (GO & KEGG)

![GO_UpRegulated.Genes.png](plot-result/GO_UpRegulated.Genes.png)
![GO_DownRegulated.Genes.png](plot-result/GO_DownRegulated.Genes.png)

Gambar 3. Analisis _Gene Ontology_ (GO) mengindikasikan bahwa gen-gen yang terekspresi diferensial terlibat dalam proses imunologis serta regulasi homeostasis metabolik.

![KEGG_Pathway.png](plot-result/KEGG_Pathway.png)

Gambar 4. Analisis _Kyoto Encyclopedia of Genes and Genomes_ (KEGG) mengidentifikasi 18 jalur yang signifikan secara statistik (adj.P.Val < 0.05).

Terdapat tiga jalur kunci yang menjelaskan perbaikan klinis pada pasien DMT2:
1. Modulasi glukoneogenesis (_PI3K-Akt signaling_)

![hsa04151.pathview.png](plot-result/hsa04151.pathview.png)

Gambar 5. _PI3K-Akt signaling_ (hsa04151): ditemukan penurunan signifikan pada gen PEPCK (_Phosphoenolpyruvate carboxykinase_), yang menyebabkan terjadinya penekanan proses glukoneogenesis di hati

2. Reprogramming metabolisme lipid (_PPAR signaling_)

![hsa03320.pathview.png](plot-result/hsa03320.pathview.png)

Gambar 6. _PPAR signaling_ (hsa03320): penurunan ekspresi transporter lipid (_Apo-AI_ dan _PLTP_) dan _Perilipin_ menunjukkan adaptasi mobilisasi lemak tubuh.

3. Adaptasi penyerapan nutrisi (_vitamin digestion and absorption_)

![hsa04977.pathview.png](plot-result/hsa04977.pathview.png)

Gambar 7. _Vitamin digestion and absorption_ (hsa04977): penurunan ekspresi gen transporter mikronutrien (_PCFT_ dan _RFC_) memberikan bukti molekuler terkait risiko malabsorpsi folat pasca-bypass.

# 4. Kesimpulan
Efektivitas operasi RYGB dalam mengatasi diabetes melibatkan mekanisme multi-jalur. Selain restrukturisasi respon imun, operasi ini secara spesifik menekan jalur glukoneogenesis dan mengubah dinamika transportasi nutrisi. Temuan ini menegaskan pentingnya manajemen nutrisi jangka panjang bagi pasien pasca-operasi untuk mengantisipasi defisiensi vitamin.
