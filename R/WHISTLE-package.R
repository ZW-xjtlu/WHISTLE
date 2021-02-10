#' WHISTLE package for the annotation of hand-crafted transcriptomic features used in the prediction of epitranscriptomic modifications.
#'
#' @description
#'
#' \bold{WHISTLE} is an R package for annotating hand-crafted transcriptomic features that are powerful in the machine learning prediction of the epitranscriptomics markers.
#' The package provides the computational representation of the known topological and functional genomic features related to the mRNA modifications.
#' Users are required to provide the range based site information in \code{\link{GRanges}} format.
#' Then, WhisleR will automatically retrieve transcriptomic features from the corresponding databases to return a summarized feature matrix.
#' In addition, the feature set annotated can be expanded and tailored on the specific application domain according to user's requirement.
#'
#' For specific example of the package usage, see the package vignette, by typing \code{browseVignettes("WHISTLE")}.
#'
#' The package is still under active development, and we welcome all biologists and bioinformaticians for all kinds of communications and collaborations. Please feel free to contact Dr. Zhen Wei <zhen.wei@xjtlu.edu.cn> if you have any questions.
#'
#'
#' @author
#' Zhen Wei
#'
#' @references
#'
#' WHISTLE reference:
#'
#' Chen K, Wei Z, Zhang Q, Wu X, Rong R, Lu Z, Su J, de Magalhães JP, Rigden DJ, Meng J: \bold{WHISTLE: a high-accuracy map of the human N6-methyladenosine (m6A) epitranscriptome predicted using a machine learning approach.} Nucleic Acids Research 2019.
#'
#' @examples
#' # For detailed information of on usage, please check the main function with:
#' ?transcriptomicFeatures
#'
#' @docType package
#' @name WHISTLE-package
NULL
