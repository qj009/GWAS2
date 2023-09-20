#' SUDO letter coded genotype matrix creation function
#' @description
#' This function allows you to convert numeric coded genptype matrix into a sudo letter coded genotype matrix. It can only be applied for sample which are pure homozyous.

#' @details
#' Additional details...
#'
#'
#' @param gencode.a Genotype matrix with numeric code: -1/1. The rows represent SNPs. The columns represent samples.

#' @returns
#' It outputs sudo letter coded genotype matrix in list format.


#' @keywords
#' @export
#' @examples
#' SUDO.GG(gencode.a)

SUDO.GG <- function(gencode.a){
  gencode.a[which(gencode.a==-1)]="A"
  gencode.a[which(gencode.a==1)]="C"
  return(gencode.a)
}
