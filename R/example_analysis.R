library(limma)

# 1. 500 gen ve 6 örnekli (3 Kontrol, 3 Tümör) simüle veri oluşturalım
set.seed(123)
expr_data <- matrix(rnorm(500 * 6, mean = 8, sd = 1.2), nrow = 500, ncol = 6)
rownames(expr_data) <- paste0("Gene_", 1:500)
colnames(expr_data) <- c("Ctrl1", "Ctrl2", "Ctrl3", "Tumor1", "Tumor2", "Tumor3")

# İlk 30 genin tümörde ifadesini yapay olarak artıralım/azaltalım (sinyal olsun diye)
expr_data[1:15, 4:6] <- expr_data[1:15, 4:6] + 2.5   # Up-regulated genler
expr_data[16:30, 4:6] <- expr_data[16:30, 4:6] - 2.5 # Down-regulated genler

# 2. Deney Tasarımı (Design Matrix)
group <- factor(c("Control", "Control", "Control", "Tumor", "Tumor", "Tumor"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# 3. Modeli Oturtma
fit <- lmFit(expr_data, design)

# 4. Kontrast: Tümör vs Kontrol
contrast_matrix <- makeContrasts(Tumor_vs_Control = Tumor - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)

# 5. Empirical Bayes Düzeltmesi
fit2 <- eBayes(fit2)

# 6. Tüm sonuçları bir tabloya çekelim
results <- topTable(fit2, coef = "Tumor_vs_Control", number = Inf, adjust.method = "BH")
head(results, 10)
