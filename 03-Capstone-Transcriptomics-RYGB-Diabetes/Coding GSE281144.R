# PROJECT: Analisis Ekspresi Gen Adaptasi Usus Pasca Gastric Bypass Pada Pasien Diabetes Tipe 2
# DATASET: GSE281144 
# PLATFORM: Microarray (Affymetrix Human Transcriptome Array 2.0)
# DESCRIPTION: Analisis Diferensial Ekspresi Gen (DEG) untuk membandingkan 
# kondisi Pre-Op vs Post-Op (Gabungan 1 & 6 bulan) menggunakan R/Bioconductor.


# ------------------------------------------------------------------------------
# PART A. PERSIAPAN LINGKUNGAN (ENVIRONMENT SETUP)
# Sebelum analisis dimulai, kita harus memastikan semua 'library' khusus 
# bioinformatika tersedia dengan menggunakan BiocManager untuk mengelola paket.
# ------------------------------------------------------------------------------

# Memastikan BiocManager terinstal untuk akses ke repository Bioconductor 
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# GEOquery: Digunakan untuk parsing data dari NCBI Gene Expression Omnibus
# limma: Menggunakan model linear untuk identifikasi ekspresi gen diferensial 
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE)

# Paket visualisasi untuk membuat plot yang estetis (ggplot2, pheatmap) 
install.packages(c("pheatmap", "ggplot2", "dplyr", "umap"))

library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(AnnotationDbi)
library(umap)



# ------------------------------------------------------------------------------
# PART B. IMPORT DATA DARI DATABASE GEO
# Langkah ini mengambil 'ExpressionSet' yang berisi data mentah dan metadata.
# ------------------------------------------------------------------------------

# getGEO mendownload data langsung menggunakan ID GSE281144
gset <- getGEO("GSE281144", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

#ExpressionSet berisi:
# - exprs() : matriks ekspresi gen
# - pData() : metadata sampel
# - fData() : metadata fitur (probe / gen)



# ------------------------------------------------------------------------------
# PART C. DATA PRE-PROCESSING & NORMALISASI
# Data microarray mentah seringkali memiliki varians yang tidak stabil.
# Kita melakukan transformasi log2 untuk menormalisasi distribusi data.
# ------------------------------------------------------------------------------

# exprs(): mengambil matriks ekspresi gen
# Baris  = probe/gen
# Kolom  = sampel
ex <- exprs(gset) 

# Cek apakah data sudah dalam skala log atau masih linear
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

#IF statement:
#Jika LogTransform = TRUE, maka lakukan log2
#Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
if (LogTransform) {
  ex[ex <= 0] <- NA # Menghindari nilai negatif/nol sebelum log
  ex <- log2(ex)    # Transformasi log2 agar distribusi data lebih 'normal'
}



# ------------------------------------------------------------------------------
# PART D. PENYARINGAN DATA & DEFINISI KELOMPOK (SUBSETTING)
# Pada bagian ini, kita melakukan subsetting untuk hanya mengambil pasien 
# dengan status 'Diabetic' dan membandingkan kondisi sebelum vs sesudah operasi.
# ------------------------------------------------------------------------------

#Dari ExpressionSet, perlu tahu mana kolom spesifik mana yang berisikan info biologis
View(pData(gset))

# Identifikasi baris yang memiliki status 'diabetic'
# Kita mencari kata 'diabetic' pada kolom characteristics_ch1.1
is_diabetic <- grepl("diabetic", pData(gset)$characteristics_ch1.1, ignore.case = TRUE)

# Filter Objek gset (hanya mengambil kolom sampel diabetes)
# Pastikan untuk update objek 'gset' yang hanya berisi sampel pasien diabetes
# Pastikan juga update matriks ekspresi 'ex' agar hanya berisi sampel diabetes yang sudah difilter
gset <- gset[, is_diabetic]
ex <- exprs(gset)

# Mengambil informasi waktu (Pre-Op vs Post-Op) dari kolom 'title'
# baseline = Pre-Op, sisanya (1 month/6 month) = Post-Op
titles <- pData(gset)$title
group_info <- ifelse(grepl("baseline", titles, ignore.case = TRUE), "PreOp_DM", "PostOp_DM")

# Membersihkan format teks menjadi faktor agar dikenali oleh model statistik R
groups <- make.names(group_info)
gset$group <- factor(groups)
nama_grup <- levels(gset$group)

# Kontras: Membandingkan efek SETELAH operasi terhadap SEBELUM operasi
# Post-Op (Treatment) vs Pre-Op (Control)
# Menentukan urutan yang akan dibandingkan 
print(nama_grup)

grup_post <- nama_grup[1] 
grup_pre  <- nama_grup[2]



# ------------------------------------------------------------------------------
# PART E. ANALISIS STATISTIK DIFFERENTIAL EXPRESSION (LIMMA)
# Menggunakan Linear Models for Microarray (limma) untuk menemukan gen yang 
# berubah secara signifikan di antara kedua kelompok.
# ------------------------------------------------------------------------------

# Membuat Matriks Desain: Merepresentasikan struktur eksperimen secara matematis
design <- model.matrix(~0 + gset$group)
colnames(design) <- levels(gset$group)

# Fit Model: Menghitung rata-rata ekspresi tiap gen pada setiap grup 
fit <- lmFit(ex, design)

# Create Contrast: Mendefinisikan perbandingan spesifik (Post-Op minus Pre-Op)
contrast_formula <- paste(grup_post, "-", grup_pre, sep = "")
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

# eBayes: Langkah krusial untuk menstabilkan estimasi varians gen menggunakan 
# metode statistik Bayesian, sehingga hasil lebih akurat meski jumlah sampel terbatas
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Mengambil tabel hasil akhir (Top Table) dengan koreksi FDR (False Discovery Rate) 
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.05
)
head(topTableResults)



# ------------------------------------------------------------------------------
# PART F. ANOTASI BIOLOGIS (GENE MAPPING)
# Data awal berupa ID Probe (kode teknis chip). Kita mapping ke Symbol Gen resmi.
# Karena kendala versi pada library '.db' spesifik, kita menggunakan data anotasi 
# yang disertakan langsung oleh GEO (fData).
# ------------------------------------------------------------------------------

# Mengambil metadata fitur (gen) langsung dari objek gset
feature_data <- fData(gset)

# Mapping menggunakan data dari GEO 
# Cek informasi biologis (simbol gen) ternyata disimpan di dalam kolom bernama apa?
colnames(fData(gset))

# Simbol gen ternyata tersimpan di dalam kolom bernama 'gene_assignment'
# Kolom 'gene_assignment' ini biasanya berisi teks panjang 
# Menggabungkan Accession Number, Simbol Gen, dan Deskripsi. 
# Kita perlu sedikit melakukan "pembedahan" teks 

# Mengambil SYMBOL (biasanya di bagian kedua setelah //)
extracted_symbols <- sub("^[^//]* // ([^//]*) // .*$", "\\1", feature_data$gene_assignment)

# Mengambil GENE NAME (biasanya di bagian ketiga setelah //)
extracted_names <- sub("^[^//]* // [^//]* // ([^//]*) // .*$", "\\1", feature_data$gene_assignment)

# Mapping probe -> gene symbol & gene name
anno_map <- data.frame(
  PROBEID = feature_data$ID,
  SYMBOL = extracted_symbols,
  GENENAME = extracted_names,
  stringsAsFactors = FALSE
)

# Gabungkan ke hasil statistik topTable
topTableResults <- merge(topTableResults, anno_map, by = "PROBEID", all.x = TRUE)

# Pembersihan: Jika tidak ada nama gen (ditandai "---"), ubah jadi NA
topTableResults$SYMBOL[topTableResults$SYMBOL == "---"] <- NA
topTableResults$GENENAME[topTableResults$GENENAME == "---"] <- NA

# Cek hasil anotasi
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])



# ------------------------------------------------------------------------------
# PART G. QUALITY CONTROL & EKSPLORASI DATA (VISUALISASI AWAL)
# Sebelum masuk ke hasil akhir, kita perlu memastikan kualitas data dan 
# melihat bagaimana sampel terpisah secara alami.
# ------------------------------------------------------------------------------

# 1. BOXPLOT: Mengecek distribusi nilai ekspresi antar sampel
# Penting untuk memastikan tidak ada sampel yang memiliki rentang nilai yang sangat berbeda (batch effect).
# Set warna berdasarkan grup
group_colors <- as.numeric(gset$group) 

boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Distribusi Nilai Ekspresi: Post-Op vs Pre-Op",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.7
)

# 2. DENSITY PLOT: Melihat sebaran global nilai ekspresi
# Plot ini memastikan bahwa proses normalisasi log2 telah berhasil menyamakan distribusi antar grup.
# Gabungkan ekspresi & grup ke data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Kerapatan Distribusi (Density Plot)",
    x = "Log2 Expression",
    y = "Density"
  )

# 3. UMAP: Visualisasi Dimensi Rendah
# UMAP membantu kita melihat apakah sampel Pre-Op dan Post-Op mengelompok (cluster) 
# berdasarkan kemiripan ekspresi gen global mereka.

# Transpose matriks ekspresi:
# UMAP bekerja dengan menghitung kemiripan antar sampel (baris)
umap_input <- t(ex) 

# Jalankan UMAP
umap_result <- umap(umap_input)

# Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

#Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP: Pengelompokan Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )



# ------------------------------------------------------------------------------
# PART H. VISUALISASI DATA (INTERPRETASI HASIL)
# ------------------------------------------------------------------------------

# 1. VOLCANO PLOT: Memberikan gambaran global perubahan gen (Log Fold Change vs Signifikansi) 
volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

#Klasifikasi status gen
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.05] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.05] <- "DOWN"

#Visualisasi
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Differential Gene Expression: Post-Op vs Pre-Op", 
       subtitle = "Dataset GSE281144 - Gastric Bypass Adaptation")

# 2. HEATMAP: Melihat pola ekspresi 50 gen paling signifikan di setiap sampel (Top 50 DEGs)
# Pilih 50 gen paling signifikan berdasarkan adj.P.Val
# Kita hanya ambil baris yang SYMBOL-nya tidak NA dan bukan "---"
topTableAnnotated <- subset(topTableResults, !is.na(SYMBOL) & SYMBOL != "" & SYMBOL != "---")
topTableAnnotated <- topTableAnnotated[order(topTableAnnotated$adj.P.Val), ]
top50 <- head(topTableAnnotated, 50)

# Ambil matriks ekspresi untuk gen terpilih menggunakan PROBEID
mat_heatmap <- ex[top50$PROBEID, ]

# Karena sudah kita filter, kita langsung pakai SYMBOL sebagai label (tanpa fallback)
gene_label <- top50$SYMBOL
rownames(mat_heatmap) <- gene_label

# Membersihkan data dari NA atau varians nol agar heatmap tidak error
mat_heatmap <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, ]
mat_heatmap <- mat_heatmap[apply(mat_heatmap, 1, var) > 0, ]

# Anotasi kolom (kelompok sampel)
annotation_col <- data.frame(
  Group = gset$group
)
rownames(annotation_col) <- colnames(mat_heatmap)

# Visualisasi heatmap 
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
# PART I. EXPORT DATA HASIL AKHIR
# Menyimpan tabel hasil analisis ke format CSV untuk pelaporan lebih lanjut.
# ------------------------------------------------------------------------------
# Mengambil daftar simbol gen yang signifikan (misal p-adj < 0.05)
# untuk dimasukkan ke analisis GO/KEGG online
sig_genes <- subset(topTableAnnotated, adj.P.Val < 0.05)$SYMBOL

# Simpan ke file teks agar mudah di-copy paste
write.table(sig_genes, "daftar_gen_signifikan.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.csv(topTableResults, "GSE281144_DEG_Results.csv")
message("Analisis selesai. File CSV telah siap digunakan!")



# ------------------------------------------------------------------------------
# PART J. ANALISIS GO 
# Analisis Gene Ontology difokuskan pada Biological Process (BP)
# Alasan: BP memberikan gambaran langsung mengenai perubahan jalur fungsional 
# dan proses sistemik yang terjadi pada pasien diabetes pasca-gastric bypass.
# ------------------------------------------------------------------------------

# Load library (pastikan sudah terinstal)
# Install jika belum ada
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"), update = FALSE)
n

library(clusterProfiler)
library(org.Hs.eg.db) # Database untuk gen manusia
library(enrichplot)
library(ggplot2)

# Mengambil gen yang signifikan dan memiliki SYMBOL
sig_genes_df <- subset(topTableResults, adj.P.Val < 0.05 & !is.na(SYMBOL) & SYMBOL != "" & SYMBOL != "---")

# Memisahkan gen yang NAIK (Up-regulated) dan TURUN (Down-regulated)
# Ini penting karena jalur biologisnya bisa sangat berbeda
up_genes <- sig_genes_df$SYMBOL[sig_genes_df$logFC > 1]
down_genes <- sig_genes_df$SYMBOL[sig_genes_df$logFC < -1]

# Cek jumlah gen
cat("Jumlah gen Up:", length(up_genes), "\n")
cat("Jumlah gen Down:", length(down_genes), "\n")
write.table(sig_genes_df$SYMBOL, "gen_list.txt", row.names=F, col.names=F, quote=F)


# 1. Analisis GO (Biological Process) untuk Up-Regulated Genes
# Menjalankan Analisis Enrichment GO (Biological Process)
# Kita gunakan OrgDb "org.Hs.eg.db" karena ini data manusia (Homo sapiens)
# Karena jumlah gen yang NAIK hanya sedikit, ambil gen signifikan tanpa filter logFC yang ketat
up_genes_more <- sig_genes_df$SYMBOL[sig_genes_df$logFC > 0.1]

# Konversi Symbol ke Entrez ID untuk Gen yang NAIK (Up)
up_entrez <- bitr(
  up_genes_more, 
  fromType="SYMBOL", 
  toType="ENTREZID", 
  OrgDb=org.Hs.eg.db
)
cat("Jumlah gen Up terkonversi:", nrow(up_entrez), "\n")

# Jalankan GO dengan p-value yang lebih ramah (0.1)
go_up <- enrichGO(
  gene          = up_entrez_more$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)
# Visualisasi dotplot untuk gen UP
dotplot(go_up_final, showCategory = 10) + 
  labs(title = "GO Enrichment: Up-Regulated Processes",
       subtitle = "Processes activated after Gastric Bypass in Diabetic Patients",
       x = "Gene Ratio",
       y = "Biological Process"
  )

# 2. Analisis GO (Biological Process) untuk Down-Regulated Genes
# Menjalankan Analisis Enrichment GO (Biological Process)
# Kita gunakan OrgDb "org.Hs.eg.db" karena ini data manusia (Homo sapiens)
# Konversi Symbol ke Entrez ID untuk Gen yang TURUN (Down)
down_entrez <- bitr(
  down_genes, 
  fromType = "SYMBOL", 
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db
)
cat("Jumlah gen Up terkonversi:", nrow(down_entrez), "\n")

# Jalankan GO untuk gen DOWN
go_down <- enrichGO(
  gene          = down_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Visualisasi Dotplot Down
dotplot(go_down, showCategory = 10) + 
  labs(title = "GO Enrichment: Down-Regulated Processes",
       subtitle = "Biological functions suppressed after Gastric Bypass",
       x = "Gene Ratio",
       y = "Biological Process"
  )

dotplot(go_up_final, showCategory = 10) + 
  labs(title = "GO Enrichment: Up-Regulated Processes",
       subtitle = "Processes activated after Gastric Bypass in Diabetic Patients",
       x = "Gene Ratio",
       y = "Biological Process"
  )



# ------------------------------------------------------------------------------
# PART K. ANALISIS KEGG PATHWAY (GABUNGAN UP & DOWN REGULATED)
# Tujuan: Melihat interaksi gen dalam peta jalur metabolisme yang utuh.
# ------------------------------------------------------------------------------

# Menyiapkan Daftar Gen Gabungan (Up + Down)
# Kita menggunakan semua gen yang signifikan (adj.P.Val < 0.05) 
# agar mendapatkan gambaran sirkuit biologis yang lengkap.
all_genes_symbol <- sig_genes_df$SYMBOL

# Konversi SYMBOL ke ENTREZ ID
# Database KEGG secara teknis hanya mengenali format ENTREZID (angka).
all_entrez_df <- bitr(
  all_genes_symbol, 
  fromType = "SYMBOL", 
  toType   = "ENTREZID", 
  OrgDb    = org.Hs.eg.db
)

# Catatan: Cek jumlah gen yang berhasil terkonversi
cat("Jumlah gen yang masuk ke analisis KEGG:", nrow(all_entrez_df), "\n")

# Menjalankan Enrichment KEGG
# 'hsa' adalah kode untuk Homo sapiens (manusia).
kegg_enrich <- enrichKEGG(
  gene         = all_entrez_df$ENTREZID,
  organism     = 'hsa', 
  pvalueCutoff = 0.05
)

# Melihat seluruh daftar hasil KEGG dalam bentuk tabel
kegg_table <- as.data.frame(kegg_enrich)
View(kegg_table) # Akan membuka jendela baru berisi daftar lengkap

# Visualisasi Barplot KEGG
# Menunjukkan jalur metabolisme mana yang paling banyak berubah (paling signifikan).
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

# Visualisasi Spesifik: Pathview (Peta Berwarna Merah-Hijau)
# Langkah ini akan menghasilkan file gambar .PNG di folder kerja R kamu.
# Merah = Gen Naik (Up), Hijau = Gen Turun (Down).
# Instal jika belum ada: 
BiocManager::install("pathview")
library(pathview)

# Menyiapkan data logFC yang sudah dipasangkan dengan Entrez ID
# Kita ambil logFC dari sig_genes_df dan beri nama sesuai Entrez ID-nya
kegg_logFC <- sig_genes_df$logFC
names(kegg_logFC) <- all_entrez_df$ENTREZID[match(sig_genes_df$SYMBOL, all_entrez_df$SYMBOL)]


# Visualisasi Jalur "Vitamin Digestion and Absorption" (hsa04977)
pathview(
  gene.data  = kegg_logFC,
  pathway.id = "hsa04977", 
  species    = "hsa",
  limit      = list(gene = 2, cpd = 1)
)

# Visualisasi Jalur "PI3K-Akt signaling" (hsa04151)
pathview(
  gene.data  = kegg_logFC,
  pathway.id = "hsa04151", 
  species    = "hsa",
  limit      = list(gene = 2, cpd = 1)
)

# Visualisasi Jalur "PPAR signaling" (hsa03320)
pathview(
  gene.data  = kegg_logFC,
  pathway.id = "hsa03320", 
  species    = "hsa",
  limit      = list(gene = 2, cpd = 1)
)



