#' ---
#' title: "Comparing patterns of gene expression: tumor v. normal colon tissue"
#' author: Gail Hobbs
#' date: "April 5, 2020"
#' geometry: margin=2cm
#' number_sections: true
#' ---
#' <style type="text/css">
#' .main-container {
#' max-width: 800px !important;
#' font-size: 18px;
#' }
#' code.r{
#' font-size: 18px;
#' }
#' pre {
#' font-size: 18px
#' }
#' h1.title {
#'   font-size: 30px;
#' }
#' h1 {
#'   font-size: 24px;
#' }
#' h2 {
#'   font-size: 18px;
#' }
#' h3 {
#'   font-size: 12px;
#' }
#' </style>

#+ setup, include = FALSE
knitr::opts_chunk$set(comment=NA, warning=FALSE, message=FALSE, size=12)

#'
#' # FWER (Holm) and FDR (Benjamini–Yekutieli)
#'
#'
#'One method for identifying differentially expressed genes is to perform a two-sided Student t-test
#'for each gene comparing the two groups, which is what is done in the corresponding paper:
#'
#' [Broad patterns of gene expression revealed by clustering analysis of tumor and normal colon tissues probed by oligonucleotide arrays (Alon 1999)](https://www.pnas.org/content/96/12/6745.full)
#' 
#' First, we run a multiple t-test on the data, and again with respect to FWER and FDR.

multiple.t.test <- function(A, grp_labels, c = "none") {
  rownames(A) <- grp_labels
  m = ncol(A)
  pval = numeric(m)
  ind = as.integer(grp_labels)

  for (i in 1:m) {
   tumor = A[(ind == 1), i]
   normal = A[(ind == 2), i]
   pval[i] = t.test(tumor, normal)$p.value
  }
  
  pval.corrected = p.adjust(pval, c)
  return(pval.corrected)
}

dat = load("alon1999.rda")
expression_data = gene_expression$expression
labels = gene_expression$status

pval = multiple.t.test(expression_data, labels) # raw
pval.holm = multiple.t.test(expression_data, labels, "holm")  # FWER
pval.by = multiple.t.test(expression_data, labels, "BY") # FDR

plot(sort(pval),
     type = "l",
     xlab = "",
     ylab = "sorted p-values")
abline(h = 0.10, lty = 2) # 0.10 level
abline(h = 0, lty = 3)
abline(h = 1, lty = 3)
lines(sort(pval.holm), type = "l", col = 2)
lines(sort(pval.by), type = "l", col = 3)
legend("bottomright",
       c("raw", "Holm", "BY"),
       col = 1:3, lty = 1, bg = "white")

#' Rejections at the 10% level, 10% FWER, and 10% FDR, respectively:
sum(pval <= 0.10)
sum(pval.holm <= 0.10)
sum(pval.by <= 0.10)

#' Running a multiple t-test on the data without adjustment at a 10% level, we reject 627 null hypotheses.
#' Given the size of this dataset, we would like to reduce our considerations much further;
#' 627/2000 is too large (30% of dataset). Therefore, a p-value adjustment is recommended.
#' We can see that controlling the FWER at 10% with the Holm method results in the smallest number
#' of rejections, less than half of the rejections from controlling the FDR with the Benjamini–Yekutieli's (BY) method.
#' 
#' Controlling FWER signifcantly narrows our pool of genes to study for differences in cancerous vs. noncancerous tissue.


#' 
#' # Association between gene expression levels
#'
#' ## Normal
#'

m = nrow(expression_data)
n = ncol(expression_data)

expression_data.normal = expression_data[(labels == "normal"), ]
spearman.normal = cor(expression_data.normal, method = "spearman")

#'
#' ## Tumor
#'

expression_data.tumor = expression_data[(labels == "tumor"), ]
spearman.tumor = cor(expression_data.tumor, method = "spearman")


#' To identify pairs of genes that behave differently in normal and tumor subjects at this stage,
#' we can simply order the pairs according to this difference. For illustration purposes, we extract the pair of
#' genes yielding the maximum difference in absolute value. On normal subjects the correlation is negative equal to
#' about 0.45, while on tumor subjects the correlation is positive equal to about 0.73. The fact that the correlation
#' is nontrivial in both cases and changes sign could be indicative that these genes are (jointly) involved in the
#' development of this particular type of tumor. All that said, it would be better to account for the multiplicity of it
#' all, as we looked at 2,000 pairs of genes.

spearman.diff = abs(spearman.normal - spearman.tumor)
which(spearman.diff == max(spearman.diff), arr.ind = TRUE)

spearman.normal[151, 91]
spearman.tumor[151, 91]
