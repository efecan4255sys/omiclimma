#' Diferansiyel Gen Ekspresyonu Analizi
#'
#' @param expr_mat Gen ekspresyon matrisi (log2 olmali)
#' @param groups Orneklerin gruplari (vektor)
#' @return topTable sonucu
#' @export
run_dge <- function(expr_mat, groups) {
  group_fac <- factor(groups)
  design <- stats::model.matrix(~ 0 + group_fac)
  colnames(design) <- levels(group_fac)

  fit <- limma::lmFit(expr_mat, design)
  contrast <- limma::makeContrasts(contrasts = paste(levels(group_fac)[2], levels(group_fac)[1], sep = " - "), levels = design)
  fit_con <- limma::contrasts.fit(fit, contrast)
  fit_eb <- limma::eBayes(fit_con)

  limma::topTable(fit_eb, number = Inf, adjust.method = "BH")
}
