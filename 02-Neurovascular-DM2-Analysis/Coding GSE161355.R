# Modul: Profil Ekspresi Gen Unit Neurovaskular pada Diabetes Melitus Tipe 2
# Dataset: GSE161355 (Diabetes Tipe 2 vs Normal)
# Platform: Microarray (Affymetrix Human Genome U133 Plus 2.0 Array - GPL570)
# Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG) 


# ------------------------------------------------------------------------------
# PART A. PERSIAPAN LINGKUNGAN KERJA 
# Sebelum analisis dimulai, kita harus memastikan semua 'library' khusus 
# bioinformatika tersedia dengan menggunakan BiocManager untuk mengelola paket.
# ------------------------------------------------------------------------------
# Memastikan BiocManager terinstal untuk akses ke repository Bioconductor 
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Instalasi paket utama (GEOquery untuk data GEO, limma untuk analisis statistik)
BiocManager::install(c("GEOquery", "limma", "hgu133plus2.db", "AnnotationDbi"), ask = FALSE, update = FALSE)
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
# PART B. PENGAMBILAN DATA DARI GEO
# Langkah ini mengambil 'ExpressionSet' yang berisi data mentah dan metadata.
# ------------------------------------------------------------------------------

# getGEO(): fungsi untuk mengunduh dataset berdasarkan ID GEO.
# GSEMatrix = TRUE -> data diambil dalam format ExpressionSet.

gset <- getGEO("GSE161355", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]



# ------------------------------------------------------------------------------
# PART C. DATA PRE-PROCESSING & NORMALISASI
# Data microarray mentah seringkali memiliki varians yang tidak stabil.
# Kita melakukan transformasi log2 untuk menormalisasi distribusi data.
# ------------------------------------------------------------------------------

# exprs(): mengambil matriks ekspresi gen
# Baris  = probe/gen
# Kolom  = sampel

ex <- exprs(gset)

# Cek apakah data membutuhkan Log transformasi
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

#IF statement:
#Jika LogTransform = TRUE, maka lakukan log2
#Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}



# ------------------------------------------------------------------------------
# PART D. DEFINISI KELOMPOK DAN FILTERING SAMPEL
# ------------------------------------------------------------------------------

# Dari ExpressionSet, identifikasi kolom metadata yang berisi informasi biologis
View(pData(gset))
unique(pData(gset)[["source_name_ch1"]])

# Setelah kolom yang relevan diketahui, buat variabel 'group_info'
# untuk menyimpan informasi kondisi biologis setiap sampel
group_info <- pData(gset)[["source_name_ch1"]]

# Mengonversi isi kolom menjadi format yang valid sebagai nama di R
groups <- make.names(group_info)

# Mengubah data menjadi faktor (variabel kategorikal)
# sehingga R dapat mengelompokkan data ke dalam level tertentu
gset$group <- factor(groups)

# Pengecekan Level Grup: Memastikan metadata terbaca dengan benar oleh R
# sebelum melangkah ke tahap analisis statistik yang lebih kompleks.
nama_grup <- levels(gset$group)
print(nama_grup)

# Pilih kelompok yang akan dianalisis
# Hanya ingin membandingkan Endotel dan Neuron (tanpa Astrosit)
grup_fokus <- c("Control.Endothelial.cells", "Diabetic.Endothelial.cells", 
                "Control.Neurones", "Diabetic.Neurones")

# Memfilter ExpressionSet agar hanya mencakup sampel dari kelompok yang dipilih
gset_filtered <- gset[, gset$group %in% grup_fokus]

# Memperbarui faktor untuk menghapus level yang tidak digunakan (misalnya Astrocytes)
# serta memastikan urutan level sesuai dengan yang diinginkan
gset_filtered$group <- factor(gset_filtered$group, levels = grup_fokus)



# ------------------------------------------------------------------------------
# PART E. ANALISIS STATISTIK DIFFERENTIAL EXPRESSION (LIMMA)
# Menggunakan Linear Models for Microarray (limma) untuk menemukan gen yang 
# berubah secara signifikan di antara kedua kelompok.
# ------------------------------------------------------------------------------

# Membuat Matriks Desain: Merepresentasikan struktur eksperimen secara matematis.
design <- model.matrix(~0 + gset_filtered$group)
colnames(design) <- levels(gset_filtered$group)

# Menentukan perbandingan (Diabetes vs Control) untuk masing-masing tipe sel
contrast_matrix <- makeContrasts(
  Diabetes_Efek_Endotel = Diabetic.Endothelial.cells - Control.Endothelial.cells,
  Diabetes_Efek_Neuron = Diabetic.Neurones - Control.Neurones,
  levels = design
)

# Fit Model: Menghitung rata-rata ekspresi tiap gen pada setiap grup 
fit <- lmFit(gset_filtered, design)

# eBayes: Langkah krusial untuk menstabilkan estimasi varians gen menggunakan 
# metode statistik Bayesian, sehingga hasil lebih akurat meski jumlah sampel terbatas
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Mengambil tabel hasil akhir (Top Table) dengan koreksi FDR (False Discovery Rate) 
res_endotel <- topTable(
  fit2, 
  coef="Diabetes_Efek_Endotel", 
  adjust="fdr", 
  sort.by="B", 
  number=Inf
)
res_neuron <- topTable(
  fit2, 
  coef="Diabetes_Efek_Neuron", 
  adjust="fdr", 
  sort.by="B", 
  number=Inf
)

# Ekstraksi hasil (Menggunakan Raw P-Value < 0.05 karena sampel terbatas)
res_endotel_sig <- res_endotel[res_endotel$P.Value < 0.05, ]
res_neuron_sig <- res_neuron[res_neuron$P.Value < 0.05, ]

# --- CATATAN ANALISIS: PEMILIHAN AMBANG BATAS SIGNIFIKANSI ---
# Pada dataset ini (GSE161355), jumlah sampel per kelompok terbatas (n=6).
# Penggunaan koreksi Multiple Testing (FDR/adj.P.Val) menghasilkan nilai yang 
# sangat konservatif (mendekati 1), sehingga tidak ada gen yang lolos filter ketat.
# Untuk tujuan eksplorasi biologis dan pembuatan profil ekspresi (Heatmap/Volcano), 
# saya menggunakan 'Raw P-Value < 0.05' sebagai ambang batas. 



# ------------------------------------------------------------------------------
# PART F. ANOTASI NAMA GEN
# Melakukan mapping dari Probe ID ke Gene Symbol menggunakan database hgu133plus2.db.
# ------------------------------------------------------------------------------

# 1. Anotasi untuk Endotel
# Mengambil ID probe dari hasil signifikan Endotel
probe_ids_endotel <- rownames(res_endotel_sig)

# Mapping probe -> gene symbol & gene name
# Gunakan library anotasi yang sesuai dengan platform (Affymetrix U133 Plus 2.0)
# Biasanya menggunakan package 'hgu133plus2.db'
gene_ann_endotel <- AnnotationDbi::select(
  hgu133plus2.db, 
  keys = probe_ids_endotel, 
  columns = c("SYMBOL", "GENENAME"), 
  keytype = "PROBEID"
)

# Menambahkan kolom PROBEID dari rownames untuk keperluan penggabungan
res_endotel_sig$PROBEID <- rownames(res_endotel_sig)

# Menggabungkan hasil analisis limma dengan anotasi gen
# berdasarkan PROBEID (probe ID)
# all.x = TRUE memastikan semua hasil statistik tetap dipertahankan
res_endotel_sig <- merge(
  res_endotel_sig, 
  gene_ann_endotel, 
  by = "PROBEID", 
  all.x = TRUE
)


# 2. Anotasi untuk Neuron
# Mengambil ID probe dari hasil signifikan Neuron
probe_ids_neuron <- rownames(res_neuron_sig)

# Mapping probe -> gene symbol & gene name
# Gunakan library anotasi yang sesuai dengan platform (Affymetrix U133 Plus 2.0)
# Biasanya menggunakan package 'hgu133plus2.db'
gene_ann_neuron <- AnnotationDbi::select(
  hgu133plus2.db, 
  keys = probe_ids_neuron, 
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

# Menambahkan kolom PROBEID dari rownames untuk keperluan penggabungan
res_neuron_sig$PROBEID <- rownames(res_neuron_sig)

# Menggabungkan hasil analisis limma dengan anotasi gen
# berdasarkan PROBEID (probe ID)
# all.x = TRUE memastikan semua hasil statistik tetap dipertahankan
res_neuron_sig <- merge(
  res_neuron_sig, 
  gene_ann_neuron, 
  by = "PROBEID", 
  all.x = TRUE
)



# ------------------------------------------------------------------------------
# PART G. QUALITY CONTROL & EKSPLORASI DATA (VISUALISASI AWAL)
# Sebelum masuk ke hasil akhir, kita perlu memastikan kualitas data dan 
# melihat bagaimana sampel terpisah secara alami.
# ------------------------------------------------------------------------------

# 1. BOXPLOT: Mengecek distribusi nilai ekspresi antar sampel
# Penting untuk memastikan tidak ada sampel yang memiliki rentang nilai yang sangat berbeda (batch effect).
# Set warna berdasarkan grup
group_colors <- as.numeric(gset_filtered$group)+1

# Mengambil matriks ekspresi dari data yang sudah difilter
ex_filtered <- exprs(gset_filtered)

par(mar = c(8, 4, 4, 12), xpd = TRUE)
boxplot(
  ex_filtered, 
  col = group_colors, 
  las = 2, 
  outline = FALSE, 
  cex.axis = 0.7,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel (Endotel & Neuron)",
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

# 2. DENSITY PLOT: Melihat sebaran global nilai ekspresi
# Plot ini memastikan bahwa proses normalisasi log2 telah berhasil menyamakan distribusi antar grup.
# Gabungkan ekspresi & grup ke data frame
expr_long <- data.frame(
  Expression = as.vector(ex_filtered),
  Group = rep(gset_filtered$group, each = nrow(ex_filtered))
) 

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") + 
  labs(
    title = "Distribusi Global Ekspresi Gen (Endotel vs Neuron)",
    subtitle = "Dataset: GSE161355",
    x = "Expression Value (log2)",
    y = "Density"
  )

# 3. UMAP: Visualisasi Dimensi Rendah
# UMAP memastikan bahwa profil genetik sel Endotel dan Neuron memang terpisah secara alami, 
# serta melihat seberapa kuat pengaruh Diabetes terhadap perubahan ekspresi gen di masing-masing tipe sel.

# Transpose matriks ekspresi:
# UMAP bekerja dengan menghitung kemiripan antar sampel (baris)
umap_input <- t(ex_filtered)

# Jalankan UMAP
umap_result <- umap(umap_input)

# Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[,1], 
  UMAP2 = umap_result$layout[,2], 
  Group = gset_filtered$group
)

# Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() + 
  labs(
    title = "UMAP Plot: Pengelompokan Sampel Global",
    x = "UMAP 1",
    y = "UMAP 2"
)



# ------------------------------------------------------------------------------
# PART H. VISUALISASI DATA (INTERPRETASI HASIL)
# ------------------------------------------------------------------------------

# 1. VOLCANO PLOT: Memberikan gambaran global perubahan gen (Log Fold Change vs Signifikansi)

# Volcano Plot: Endotel (Threshold: P < 0.05 & |logFC| > 0.5)
res_endotel$status <- "NO"
res_endotel$status[res_endotel$logFC > 0.5 & res_endotel$P.Value < 0.05] <- "UP"
res_endotel$status[res_endotel$logFC < -0.5 & res_endotel$P.Value < 0.05] <- "DOWN"

ggplot(res_endotel, aes(x = logFC, y = -log10(P.Value), color = status)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() + 
  labs(title = "Volcano Plot: Endothelial Cells", subtitle = "Diabetic vs Normal Control")

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
  labs(title = "Volcano Plot: Neuronal Cells", subtitle = "Diabetic vs Normal Control")


# 2. VENN DIAGRAM: Memberikan gambaran mengenai jumlah gen DEGs 
# pada masing-masing tipe sel (endotel & neuron), serta jumlah gen yang overlap di antara keduanya.

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


# 3. HEATMAP: Melihat pola ekspresi Common Genes (Gen yang sama-sama muncul di kedua sel)
# Sebelum membuat Heatmap, identifikasi daftar simbol gen yang sama-sama muncul di kedua sel (Common Genes)
# Menggunakan kriteria P < 0.05 & |logFC| > 0.5
common_genes <- intersect(
  res_endotel$Gene.symbol[res_endotel$P.Value < 0.05 & abs(res_endotel$logFC) > 0.5],
  res_neuron$Gene.symbol[res_neuron$P.Value < 0.05 & abs(res_neuron$logFC) > 0.5]
)

# Membersihkan list: Buang baris kosong "" dan NA agar heatmap tidak error
common_genes <- common_genes[common_genes != "" & !is.na(common_genes)]

# Membuat matriks heatmap langsung dari ex_filtered
# Mengambil hanya 1 probe per gen (unique) agar plot tidak terlalu rapat/numpuk
mat_heatmap <- ex_filtered[rownames(res_endotel[res_endotel$Gene.symbol %in% common_genes, ]), ]
rownames(mat_heatmap) <- res_endotel$Gene.symbol[res_endotel$Gene.symbol %in% common_genes]
mat_heatmap <- mat_heatmap[!duplicated(rownames(mat_heatmap)), ] # Unique Gene Representation

# Visualisasi Heatmap Final
pheatmap(
  mat_heatmap,
  scale = "row", 
  clustering_method = "ward.D2", # Kelompok warna lebih kontras & rapi
  annotation_col = data.frame(Group = gset_filtered$group, row.names = colnames(ex_filtered)),
  show_colnames = FALSE,
  show_rownames = TRUE, 
  fontsize_row = 7, 
  main = "Heatmap Common DEGs (Strict Filter: P < 0.05 & |FC| > 0.5)",
  color = colorRampPalette(c("#2980b9", "white", "#c0392b"))(100), # Pro Blue-Red Gradient
  border_color = NA
)



# ------------------------------------------------------------------------------
# PART I. ANALISIS GENE ONTOLOGY (GO)
# Analisis GO difokuskan pada Biological Process (BP).
# Alasan: BP memberikan gambaran langsung mengenai perubahan jalur fungsional 
# dan proses sistemik (seperti inflamasi atau metabolisme) yang terjadi pada 
# sel endotel dan neuron akibat paparan kondisi diabetik kronis.
# ------------------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Konversi Gene Symbol ke Entrez ID (Standar database GO/KEGG)
# Menggunakan data 'common_genes' hasil irisan Venn sebelumnya
genes_entrez <- bitr(common_genes, 
                     fromType = "SYMBOL", 
                     toType   = "ENTREZID", 
                     OrgDb    = org.Hs.eg.db)

# Pengayaan GO: Biological Process (BP)
# Menemukan proses biologis utama yang paling terpengaruh di kedua populasi sel
go_bp <- enrichGO(gene           = genes_entrez$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP", 
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.15,  #Optimal threshold: Menghasilkan 17 jalur fungsional utama
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)

# Visualisasi dotplot (Menampilkan 15 jalur paling signifikan)
dotplot(go_bp, showCategory = 15) + 
  labs(
    title = "Top Biological Processes: Common DEGs in Diabetic Cells",
    subtitle = "Shared mechanisms between Endothelial and Neuronal populations",
    x = "Gene Ratio",
    y = "Biological Process"
  )
write.csv(as.data.frame(go_bp), "GO_Enrichment_Results.csv")


# ------------------------------------------------------------------------------
# PART J. ANALISIS JALUR KEGG (KYOTO ENCYCLOPEDIA OF GENES AND GENOMES)
# Analisis KEGG digunakan untuk memetakan gen-gen ke dalam jalur sinyal metabolik 
# dan patofisiologi penyakit yang terintegrasi secara global.
# Threshold: P-Value < 0.2 digunakan untuk mengeksplorasi jalur fungsional utama.
# ------------------------------------------------------------------------------

# Pengayaan Jalur KEGG
# 'hsa' merujuk pada Homo sapiens (Manusia)
kegg <- enrichKEGG(
  gene = genes_entrez$ENTREZID,
  organism = 'hsa', 
  pvalueCutoff = 0.2
)

# Melihat seluruh daftar hasil KEGG dalam bentuk tabel
kegg_table <- as.data.frame(kegg)
View(kegg_table) # Akan membuka jendela baru berisi daftar lengkap


# Visualisasi Barplot untuk Jalur KEGG
barplot(kegg, showCategory = 10) + 
  labs(
    title = "Top KEGG Pathways: Shared Mechanisms in Diabetic Cells",
    subtitle = "Key metabolic and signaling pathways (P-value < 0.2)",
    x = "Gene Count",
    y = "Pathway Description"
  ) +
  theme_minimal()
head(as.data.frame(kegg))


# Visualisasi Spesifik: Pathview (Peta Berwarna Merah-Hijau)
# Merah = Gen Naik (Up), Hijau = Gen Turun (Down).
library(pathview)

# Ambil data logFC untuk gen irisan dari res_endotel
# Kita gunakan Entrez ID sebagai nama (karena pathview butuh Entrez ID)
kegg_logFC <- res_endotel$logFC[res_endotel$Gene.symbol %in% common_genes]
names(kegg_logFC) <- genes_entrez$ENTREZID[match(common_genes, genes_entrez$SYMBOL)]

# Hapus jika ada NA
kegg_logFC <- kegg_logFC[!is.na(names(kegg_logFC))]


# VISUALISASI PATHWAYS
# 1. Visualisasi Jalur "Ferroptosis" (hsa04216)
pathview(
  gene.data  = kegg_logFC,
  pathway.id = "hsa04216", 
  species    = "hsa",
  limit      = list(gene = 2, cpd = 1),
  low        = "blue",              # Down-regulated
  mid        = "gray",              # Netral
  high       = "red"                # Up-regulated
)

# 2. Visualisasi Jalur "Focal Adhesion" (hsa04510)
pathview(
  gene.data  = kegg_logFC,
  pathway.id = "hsa04510", 
  species    = "hsa",
  limit      = list(gene = 2, cpd = 1),
  low        = "blue",              # Down-regulated
  mid        = "gray",              # Netral
  high       = "red"                # Up-regulated
)

# 3. Visualisasi Jalur "Adherens Junction" (hsa04520)
pathview(
  gene.data  = kegg_logFC,
  pathway.id = "hsa04520", 
  species    = "hsa",
  limit      = list(gene = 2, cpd = 1),
  low        = "blue",              # Down-regulated
  mid        = "gray",              # Netral
  high       = "red"                # Up-regulated
)

# 4. Visualisasi Jalur "Integrin signaling" (hsa04518)
pathview(
  gene.data  = kegg_logFC,
  pathway.id = "hsa04518", 
  species    = "hsa",
  limit      = list(gene = 2, cpd = 1),
  low        = "blue",              # Down-regulated
  mid        = "gray",              # Netral
  high       = "red"                # Up-regulated
)



# ------------------------------------------------------------------------------
# PART K. MENYIMPAN HASIL
# ------------------------------------------------------------------------------

write.csv(res_endotel_sig, "Hasil_DEGs_Endotel_Diabetes.csv", row.names = FALSE)
write.csv(res_neuron_sig, "Hasil_DEGs_Neuron_Diabetes.csv", row.names = FALSE)
message("Analisis selesai! File hasil telah disimpan.")