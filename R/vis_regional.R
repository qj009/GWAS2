#' GWAS visualization function: vis

#' @description
#' This function allows you to make manhattan plot containing SNP-based GWAS, initial haplotype scan and Haplotype-based GWAS result.

#' @details
#' Input file p value column should not contain missing value and 0. If the original p value is 0 (which is too small to calculate by the software), please compute p value manually or set a standard small p value for plot the figure.
#'

#' @param T.Name Phenotype name, shown in the plot name;
#' @param snp_file SNP-based GWAS result. The table contain four column: chromosome, base pair position, p value and SNP ID;
#' @param hapi_file Initial di-SNP GWAS result. The table contain five column: chromosome, start base pair position, end base pair position,p value and haplotype ID;
#' @param hap_file Haplotype-based GWAS final result. The table contain four column: chromosome, start base pair position, end base pair position, p value and haplotype ID;
#' @param sig_line Significant p value threshold;
#' @param ylim The range of y-axis. The default is NULL for default scale range. It also can be a numeric vector of length two providing limits of the sale;
#' @param snp_sig_size Point size of significant SNP signals, default is 8;
#' @param snp_sig_alpha Point transparency of significant SNP signals, default is 0.8;
#' @param snp_nosig_size Point size of non-significant SNP signals, default is 1;
#' @param snp_no_sig_alpha Point transparency of non-significant SNP signals, default is 0.3;
#' @param hapi_sig_lineend Line end style (round, butt, square) of initial haplotype significant signals, default is "round";
#' @param hapi_sig_linewidth Line width of initial haplotype significant signals, default is 5;
#' @param hapi_sig_alpha Line transparency of initial haplotype significant signals, default is 1;
#' @param hapi_nosig_lineend Line end style (round, butt, square) of initial haplotype non-significant signals, default is "round";
#' @param hapi_nosig_linewidth Line width of initial haplotype non-significant signals, default is 1;
#' @param hapi_nosig_alpha Line transparence of initial haplotype non-significant signals, default is 0.3;
#' @param hap_sig_lineend Line end style (round, butt, square) of haplotype significant signals, default is "round";
#' @param hap_sig_linewidth Line width of haplotype significant signals, default is 3;
#' @param hap_sig_alpha Line transparence of haplotype significant signals, default is 1;
#' @param hap_nosig_lineend Line end style (round, butt, square) of haplotype non-significant signals, default is "round";
#' @param hap_nosig_linewidth Line width of haplotype non-significant signals, , default is 1;
#' @param hap_nosig_alpha Line transparence of haplotype non-significant signals, default is 0.3;
#' @param snp_color SNP signal color, default is "#93b5cf";
#' @param hapi_color Initial haplotype signal color, default is "#ffa60f";
#' @param hap_color Haplotype signal color, default is "#f03752";
#' @param GWAS1_name SNP GWAS label, default is "SNP";
#' @param GWAS2_name Initial haplotype GWAS label, default is "Initial haplotype";
#' @param GWAS3_name Haplotype GWAS label, default is "Haplotype";
#' @param xtitle x-axis title;
#' @returns manhattan plot containing SNP-based GWAS, initial haplotype scan and Haplotype-based GWAS result.


#' @keywords Manhattan plot
#' @export




vis_regional <- function(T.Name, snp_file, hapi_file, hap_file, sig_line=5e-08, ylim=NULL,
                snp_sig_size = 8, snp_sig_alpha = 0.8, snp_nosig_size=1, snp_no_sig_alpha=0.8,
                hapi_sig_lineend="round", hapi_sig_linewidth=5,hapi_sig_alpha = 1,
                hapi_nosig_lineend="round", hapi_nosig_linewidth=1,hapi_nosig_alpha =1,
                hap_sig_lineend="round", hap_sig_linewidth=3,hap_sig_alpha = 1,
                hap_nosig_lineend="round", hap_nosig_linewidth=1,hap_nosig_alpha =1,
                snp_color = "#93b5cf", hapi_color= "#ffa60f", hap_color = "#f03752", GWAS1_name = "SNP", GWAS2_name = "Initial haplotype", GWAS3_name = "Haplotype",xtitle ){
  print(paste0("pheno type is ",T.Name))

  # Example SNP data (replace this with your own SNP data)
  snp <-snp_file
  colnames(snp) <- c("CHROM","POS","P","ID")
  snp_plot <- snp[which(snp$P!=0),]
  snp_plot <- snp_plot %>% mutate(sig = ifelse(P<sig_line,"sig","nosig"))

  # Example di-snp data (replace this with your own interval data)
  hapi_plot <- hapi_file
  colnames(hapi_plot) <- c("CHROM","START","END","P","ID")
  hapi_plot <- hapi_plot[which(hapi_plot$P!=0),]
  hapi_plot <- hapi_plot %>% mutate(sig = ifelse(P<sig_line,"sig","nosig"))


  #haplotype
  hap_plot <- hap_file
  colnames(hap_plot) <- c("CHROM","START","END","P","ID")
  hap_plot <- hap_plot[which(hap_plot$P!=0),]
  hap_plot <- hap_plot %>% mutate(sig = ifelse(P<sig_line,"sig","nosig"))


  # Create Manhattan plot with SNPs and intervals
  #
  p <- ggplot() +
    # snp sig
    geom_point(data = snp_plot %>% dplyr::filter(sig=="sig"), aes(x = POS, y = -log10(P),color = GWAS1_name), size = snp_sig_size, alpha = snp_sig_alpha,show.legend=TRUE) +
    # snp no sig
    geom_point(data = snp_plot %>% dplyr::filter(sig=="nosig"), aes(x = POS, y = -log10(P),color = GWAS1_name), size = snp_nosig_size, alpha = snp_no_sig_alpha,show.legend	=TRUE) +
    # di-SNP sig
    geom_segment(data = hapi_plot %>% dplyr::filter(sig=="sig"), aes(x = START, xend = END, y = -log10(P), yend = -log10(P),color = GWAS2_name), lineend = hapi_sig_lineend, linewidth= hapi_sig_linewidth,alpha = hapi_sig_alpha,show.legend	=TRUE) +
    # di-SNP nosig
    geom_segment(data = hapi_plot %>% dplyr::filter(sig=="nosig"), aes(x = START, xend = END, y = -log10(P), yend = -log10(P),color = GWAS2_name), lineend = hapi_nosig_lineend,linewidth = hapi_nosig_linewidth, alpha = hapi_nosig_alpha,show.legend	=TRUE) +
    # hap sig
    geom_segment(data = hap_plot %>% dplyr::filter(sig=="sig"), aes(x = START, xend = END, y = -log10(P), yend = -log10(P),color = GWAS3_name), lineend = hap_sig_lineend,linewidth = hap_sig_linewidth,alpha = hap_sig_alpha, show.legend	=TRUE) +
    # hap no sig
    geom_segment(data = hap_plot %>% dplyr::filter(sig=="nosig"), aes(x = START, xend = END, y = -log10(P), yend = -log10(P),color = GWAS3_name), lineend = hap_nosig_lineend,linewidth = hap_nosig_linewidth,alpha = hap_nosig_alpha,show.legend	=TRUE) +
    # add significant line
    geom_hline(
      yintercept = -log10(sig_line), color = "red",
      linetype = "dashed") +

    # add x-axis label
   # scale_x_continuous(labels = paste0("chr",axis_set$CHROM),breaks = axis_set$center) +
    # set ylim
    scale_y_continuous(expand = c(0, 0), limits = ylim) +
    # add legend
    scale_color_manual(breaks = c(GWAS1_name, GWAS2_name, GWAS3_name),values = c(snp_color, hapi_color,hap_color), name = "Data Type") +
    # add plot title
    ggtitle(paste0("Manhattan plot of ",GWAS1_name, ", ", GWAS2_name," and ", GWAS3_name," GWAS for ", T.Name)) +
    # add x and y title
    labs(x = xtitle, y = expression(-log[10]*P)) +
    theme_minimal()

  return(p)
}
