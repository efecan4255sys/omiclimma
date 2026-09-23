# omiclimma

`omiclimma`, mikrodizi ve RNA-seq sayım verileri için `limma` ve `limma-voom` tabanlı diferansiyel gen ekspresyonu (DGE) analizlerini kolaylaştıran bir R paketidir.

## Kurulum

Paketi GitHub üzerinden yüklemek için:
## Temel Fonksiyonlar

- **`run_dge()`**: Log-dönüştürülmüş mikrodizi veya normalize ekspresyon verileri için standart limma hattı.
- **`run_voom_dge()`**: Ham RNA-seq sayım verileri için TMM normalizasyonu ve `limma-voom` analizi.

## Hızlı Başlangıç (RNA-seq voom Örneği)

# Örnek sayım matrisi ve gruplar
counts <- matrix(rpois(6000, lambda = 40), nrow = 1000, ncol = 6)
rownames(counts) <- paste0("Gene_", 1:1000)
colnames(counts) <- c("Ctrl1", "Ctrl2", "Ctrl3", "Tumor1", "Tumor2", "Tumor3")

gruplar <- c("Control", "Control", "Control", "Tumor", "Tumor", "Tumor")

# Analizi çalıştır
sonuclar <- run_voom_dge(counts, gruplar)
head(sonuclar)
