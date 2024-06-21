vis_sep <- function(T.Name, snp_file, hapi_file, hap_file, sig_line=5e-08, ylim=NULL,
                snp_sig_size = 3, snp_sig_alpha = 0.8, snp_nosig_size=1, snp_no_sig_alpha=0.3,
                hapi_sig_lineend="round", hapi_sig_linewidth=5,hapi_sig_alpha = 1,
                hapi_nosig_lineend="round", hapi_nosig_linewidth=1,hapi_nosig_alpha =0.3,
                hap_sig_lineend="round", hap_sig_linewidth=8,hap_sig_alpha = 1,
                hap_nosig_lineend="round", hap_nosig_linewidth=1,hap_nosig_alpha =0.3,
                snp_color = "#93b5cf", hapi_color= "#ffa60f", hap_color = "#f03752", GWAS1_name = "SNP", GWAS2_name = "Initial haplotype", GWAS3_name = "Haplotype"){
  print(paste0("pheno type is ",T.Name))

  # Example SNP data (replace this with your own SNP data)
  snp <-snp_file
  colnames(snp) <- c("CHROM","POS","P","ID")
  snp_plot <- snp[which(snp$P!=0),]
  snp_plot <- snp_plot %>% mutate(sig = ifelse(P<sig_line,"sig","nosig"))

  snp_data_cum <- snp_plot %>%
    group_by(CHROM) %>%
    dplyr::summarise(max_bp = max(POS)) %>%
    mutate(bp_add = dplyr::lag(cumsum(max_bp), default = 0)) %>%
    dplyr::select(CHROM, bp_add)

  snp_plot <- snp_plot %>% inner_join(snp_data_cum, by = "CHROM") %>%
    mutate(POS_cum = POS + bp_add)

  # Example di-snp data (replace this with your own interval data)
  hapi_plot <- hapi_file
  colnames(hapi_plot) <- c("CHROM","START","END","P","ID")
  hapi_plot <- hapi_plot[which(hapi_plot$P!=0),]
  hapi_plot <- hapi_plot %>% mutate(sig = ifelse(P<sig_line,"sig","nosig"))

  # hapi_data_cum <- hapi_plot %>%
  #   group_by(CHROM) %>%
  #   dplyr::summarise(max_bp = max(END)) %>%
  #   mutate(bp_add = dplyr::lag(cumsum(max_bp),default = 0)) %>%
  #   dplyr::select(CHROM, bp_add)

  # hapi_plot <- hapi_plot %>%
  #   inner_join(hapi_data_cum, by = "CHROM") %>%
  #   mutate(START_cum = START + bp_add, END_cum = END + bp_add)
  hapi_plot <- hapi_plot %>%
    inner_join(snp_data_cum, by = "CHROM") %>%
    mutate(START_cum = START + bp_add, END_cum = END + bp_add)

  #haplotype
  hap_plot <- hap_file
  colnames(hap_plot) <- c("CHROM","START","END","P","ID")
  hap_plot <- hap_plot[which(hap_plot$P!=0),]
  hap_plot <- hap_plot %>% mutate(sig = ifelse(P<sig_line,"sig","nosig"))
  # hap_data_cum <- hap_plot %>%
  #   group_by(CHROM) %>%
  #   dplyr::summarise(max_bp = max(END)) %>%
  #   mutate(bp_add = dplyr::lag(cumsum(max_bp),default = 0)) %>%
  #   dplyr::select(CHROM, bp_add)

  # hap_plot <- hap_plot %>%
  #  inner_join(hap_data_cum, by = "CHROM") %>%
  #  mutate(START_cum = START + bp_add, END_cum = END + bp_add)
  hap_plot <- hap_plot %>%
    inner_join(snp_data_cum, by = "CHROM") %>%
    mutate(START_cum = START + bp_add, END_cum = END + bp_add)

  # Calculate midpoints for chromosome regions

  axis_set <- snp_plot %>%
    group_by(CHROM) %>%
    dplyr::summarize(center = mean(POS_cum))

  # Cacluate chromsome bountries
  axis_sep <- snp_plot %>% group_by(CHROM) %>% dplyr::summarize(boundaries  = max(POS_cum))

  # Create Manhattan plot with SNPs and intervals
  #
  p <- ggplot() +
    # snp no sig
    geom_point(data = snp_plot %>% dplyr::filter(sig=="nosig"), aes(x = POS_cum, y = -log10(P),color = GWAS1_name), size = snp_nosig_size, alpha = snp_no_sig_alpha,show.legend	=TRUE) +
    # di-SNP nosig
    geom_segment(data = hapi_plot %>% dplyr::filter(sig=="nosig"), aes(x = START_cum, xend = END_cum, y = -log10(P), yend = -log10(P),color = GWAS2_name), lineend = hapi_nosig_lineend,linewidth = hapi_nosig_linewidth, alpha = hapi_nosig_alpha,show.legend	=TRUE) +
    # hap no sig
    geom_segment(data = hap_plot %>% dplyr::filter(sig=="nosig"), aes(x = START_cum, xend = END_cum, y = -log10(P), yend = -log10(P),color = GWAS3_name), lineend = hap_nosig_lineend,linewidth = hap_nosig_linewidth,alpha = hap_nosig_alpha,show.legend	=TRUE) +
    # hap sig
    geom_segment(data = hap_plot %>% dplyr::filter(sig=="sig"), aes(x = START_cum, xend = END_cum, y = -log10(P), yend = -log10(P),color = GWAS3_name), lineend = hap_sig_lineend,linewidth = hap_sig_linewidth,alpha = hap_sig_alpha, show.legend	=TRUE) +
    # di-SNP sig
    geom_segment(data = hapi_plot %>% dplyr::filter(sig=="sig"), aes(x = START_cum, xend = END_cum, y = -log10(P), yend = -log10(P),color = GWAS2_name), lineend = hapi_sig_lineend, linewidth= hapi_sig_linewidth,alpha = hapi_sig_alpha,show.legend	=TRUE) +
    # snp sig
    geom_point(data = snp_plot %>% dplyr::filter(sig=="sig"), aes(x = POS_cum, y = -log10(P),color = GWAS1_name), size = snp_sig_size, alpha = snp_sig_alpha,show.legend=TRUE) +
    # add significant line
    geom_hline(
      yintercept = -log10(sig_line), color = "red",
      linetype = "dashed") +
    # add chromosome separation line
    geom_vline(
      xintercept = axis_sep$boundaries, color = "grey",
      linetype = "dotted", alpha = 0.5) +
    # add x-axis label
    scale_x_continuous(labels = paste0("chr",axis_set$CHROM),breaks = axis_set$center) +
    # set ylim
    scale_y_continuous(expand = c(0, 0), limits = ylim) +
    # add legend
    scale_color_manual(breaks = c(GWAS1_name, GWAS2_name, GWAS3_name),values = c(snp_color, hapi_color,hap_color), name = "Data Type") +
    # add plot title
    ggtitle(paste0("Manhattan plot of ",GWAS1_name, ", ", GWAS2_name," and ", GWAS3_name," GWAS for ", T.Name)) +
    # add x and y title
    labs(x = "Chromsome", y = expression(-log[10]*P)) +
    theme_minimal()

  return(p)
}
