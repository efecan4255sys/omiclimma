#' RNA-seq Ham Sayim Verileri Icin limma-voom Analizi
#'
#' @param counts_mat Ham RNA-seq sayim matrisi (satirlar: genler, sutunlar: ornekler)
#' @param group_labels Orneklerin gruplari (vektor, orn: c('Control', 'Treatment'))
#' @return topTable formatinda diferansiyel ekspresyon tablosu
#' @export
run_voom_dge <- function(counts_mat, group_labels) {
  # 1. Faktor ve tasarim matrisi
  group <- factor(group_labels)
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  # 2. DGEList olustur ve dusuk sayimli genleri filtrele
  dge <- edgeR::DGEList(counts = counts_mat, group = group)
  keep <- edgeR::filterByExpr(dge, design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  # 3. Kutuphane boyutlarini normalize et (TMM yontemi)
  dge <- edgeR::normLibSizes(dge)

  # 4. voom donusumu ile sayim varyansini modelle
  v <- limma::voom(dge, design, plot = FALSE)

  # 5. Lineer model ve Bayes duzeltmesi
  fit <- limma::lmFit(v, design)
  contrasts_str <- paste(levels(group)[2], levels(group)[1], sep = " - ")
  contrast_mat <- limma::makeContrasts(contrasts = contrasts_str, levels = design)
  fit_con <- limma::contrasts.fit(fit, contrast_mat)
  fit_eb <- limma::eBayes(fit_con)

  # 6. Sonuclari dondur
  limma::topTable(fit_eb, number = Inf, adjust.method = "BH")
}
