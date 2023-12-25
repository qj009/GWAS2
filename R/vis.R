#' GWAS visualization function: vis
#' @description
#'
#' @details
#' Input file p value column should not contain missing value and 0. If the original p value is 0 (which is too small to calculate by the software), please compute p value manually or set a standard small p value for plot the figure.
#'

#' @param T.Name:
#' @param snp_file:
#' @param hapi_file:
#' @param hap_file:
#' @param delim:
#' @param sig_line:
#' @returns:


#' @keywords
#' @export
#' @examples
#'

vis <- function(T.Name, snp_file, hapi_file, hap_file, delim = "\t", sig_line=5e-08, ylim,
                snp_sig_size = 8, snp_sig_alpha = 0.8, snp_nosig_size=1, snp_no_sig_alpha=0.3,
                hapi_sig_lineend="round", hapi_sig_linewidth=5,hapi_sig_alpha = 1,
                hapi_nosig_lineend="round", hapi_nosig_linewidth=1,hapi_nosig_alpha =0.3,
                hap_sig_lineend="round", hap_sig_linewidth=3,hap_sig_alpha = 1,
                hap_nosig_lineend="round", hap_nosig_linewidth=1,hap_nosig_alpha =0.3,
                snp_color = "#93b5cf", hapi_color= "#ffa60f", hap_color = "#f03752"){
  print(paste0("pheno type is ",T.Name))

  # Example SNP data (replace this with your own SNP data)
  snp <-read_delim(snp_file,delim = delim)
  colnames(snp) <- c("CHROM","POS","P","ID")
  snp_plot <- snp[which(snp$P!=0),]
  snp_plot <- snp_plot %>% mutate(sig = ifelse(P<sig_line,"sig","nosig"))

  snp_data_cum <- snp_plot %>%
    group_by(CHROM) %>%
    dplyr::summarise(max_bp = max(POS)) %>%
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
    dplyr::select(CHROM, bp_add)

  snp_plot <- snp_plot %>%
    inner_join(snp_data_cum, by = "CHROM") %>%
    mutate(POS_cum = POS + bp_add)

  snp_plot_sub <- snp_plot[which(snp_plot$CHROM %in% c(1,2)),]

  # Example di-snp data (replace this with your own interval data)
  hapi <- read_delim(hapi_file,delim = delim)
  hapi_plot <- hapi %>% dplyr::select(chr, start,end, p)
  colnames(hapi_plot) <- c("CHROM","START","END","P")
  #hapi_plot$CHROM <- paste0("chr",hapi_plot$CHROM)
  hapi_plot <- hapi_plot[which(hapi_plot$P!=0),]
  # set.seed(2022)
  # jagger_disnp <- runif(nrow(hapi_plot), 1e-9,1e-8)
  # hapi_plot$P <- hapi_plot$P + jagger_disnp
  hapi_plot <- hapi_plot %>% mutate(sig = ifelse(P<sig_line,"sig","nosig"))

  hapi_data_cum <- hapi_plot %>%
    group_by(CHROM) %>%
    dplyr::summarise(max_bp = max(END)) %>%
    mutate(bp_add = lag(cumsum(max_bp),default = 0)) %>%
    dplyr::select(CHROM, bp_add)

  hapi_plot <- hapi_plot %>%
    inner_join(hapi_data_cum, by = "CHROM") %>%
    mutate(START_cum = START + bp_add, END_cum = END + bp_add)

  hapi_plot_sub <- hapi_plot[which(hapi_plot$CHROM %in% c(1,2)),]

  #haplotype
  hap <- read_delim(hap_file,delim = delim)
  hap_plot <- hap %>% dplyr::select(chr, start,end, p)
  colnames(hap_plot) <- c("CHROM","START","END","P")
  #hap_plot$CHROM <- paste0("chr",hap_plot$CHROM)
  hap_plot <- hap_plot[which(hap_plot$P!=0),]
  hap_plot <- hap_plot %>% mutate(sig = ifelse(P<sig_line,"sig","nosig"))
  hap_data_cum <- hap_plot %>%
    group_by(CHROM) %>%
    dplyr::summarise(max_bp = max(END)) %>%
    mutate(bp_add = lag(cumsum(max_bp),default = 0)) %>%
    dplyr::select(CHROM, bp_add)

  hap_plot <- hap_plot %>%
    inner_join(hap_data_cum, by = "CHROM") %>%
    mutate(START_cum = START + bp_add, END_cum = END + bp_add)
  hap_plot_sub <- hap_plot[which(hap_plot$CHROM %in% c(1,2)),]

  # Calculate midpoints for chromosome regions

  axis_set <- snp_plot |>
    group_by(CHROM) |>
    dplyr::summarize(center = mean(POS_cum))

  # Create Manhattan plot with SNPs and intervals
  #
  p <- ggplot() +
    # snp sig
    geom_point(data = snp_plot %>% filter(sig=="sig"), aes(x = POS_cum, y = -log10(P),color = "SNP"), size = snp_sig_size, alpha = snp_sig_alpha,show.legend=TRUE) +
    # snp no sig
    geom_point(data = snp_plot %>% filter(sig=="nosig"), aes(x = POS_cum, y = -log10(P),color = "SNP"), size = snp_nosig_size, alpha = snp_no_sig_alpha,show.legend	=TRUE) +
    # di-SNP sig
    geom_segment(data = hapi_plot %>% filter(sig=="sig"), aes(x = START_cum, xend = END_cum, y = -log10(P), yend = -log10(P),color = "di-SNP"), lineend = hapi_sig_lineend, linewidth= hapi_sig_linewidth,alpha = hapi_sig_alpha,show.legend	=TRUE) +
    # di-SNP nosig
    geom_segment(data = hapi_plot %>% filter(sig=="nosig"), aes(x = START_cum, xend = END_cum, y = -log10(P), yend = -log10(P),color = "di-SNP"), lineend = hapi_nosig_lineend,linewidth = hapi_nosig_linewidth, alpha = hapi_nosig_alpha,show.legend	=TRUE) +
    # hap sig
    geom_segment(data = hap_plot %>% filter(sig=="sig"), aes(x = START_cum, xend = END_cum, y = -log10(P), yend = -log10(P),color = "Haplotype"), lineend = hap_sig_lineend,linewidth = hap_sig_linewidth,alpha = hap_sig_alpha, show.legend	=TRUE) +
    # hap no sig
    geom_segment(data = hap_plot %>% filter(sig=="nosig"), aes(x = START_cum, xend = END_cum, y = -log10(P), yend = -log10(P),color = "Haplotype"), lineend = hap_nosig_lineend,linewidth = hap_nosig_linewidth,alpha = hap_nosig_alpha,show.legend	=TRUE) +
    # add significant line
    geom_hline(
      yintercept = -log10(sig_line), color = "red",
      linetype = "dashed") +
    # add x-axis label
    scale_x_continuous(label = paste0("chr",axis_set$CHROM),breaks = axis_set$center) +
    # set ylim
    scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
    # add legend
    scale_color_manual(breaks = c("SNP", "di-SNP", "Haplotype"),values = c(snp_color, hapi_color,hap_color), name = "Data Type") +
    # add plot title
    ggtitle(paste0("Manhattan plot of SNP, di-SNP, Haplotype of ", T.Name)) +
    # add x and y title
    labs(x = "Chromsome", y = expression(-log[10]*P)) +
    theme_minimal()

   return(p)
}
