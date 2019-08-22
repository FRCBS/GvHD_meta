
## =======================================================================================
## Modified from:
## https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/GGplot2-Manhattan-Plot-Function
## =======================================================================================


library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

gg.manhattan <- function(df, threshold, hlight, sig, sugg, col, ylims, pointsize, title) {

  # SNP CHR BP P

  # format df
  df.tmp <- df %>%     
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>%   
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate( is_annotate=ifelse(P < threshold, "yes", "no"))
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=pointsize) +
    scale_color_manual(values = rep(col, 22)) +

    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center, expand = c(0.02, 0.02) ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome", y = expression('-log'[10]*'(p)')) +
    
    # add genome-wide sig and sugg lines
    geom_hline(yintercept = -log10(sig), linetype="solid", col='red', alpha=0.8) +
    geom_hline(yintercept = -log10(sugg), linetype="dashed") +
    
    # Add highlighted points
    #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.5) +
    
    # Custom the theme:
    theme_bw(base_size = 22) +
    theme( 
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour="black", size=.4)
    )
}

