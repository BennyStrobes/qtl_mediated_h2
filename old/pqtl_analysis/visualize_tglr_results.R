args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(dplyr)
library(reshape)
library(stringr)
library(reshape2)
library(ggbeeswarm)
library(RColorBrewer)
options(warn=1)
options(bitmapType='cairo')


figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}




make_se_barplot_showing_mediated_heritabilities <- function(df, nm_variant_anno) {
	ordered_traits <- as.character(df$trait_name)

	ord <- order(df$frac_h2_med)
	df$trait_name = factor(df$trait_name, levels=ordered_traits[ord])

	print(mean(df$frac_h2_med))
	
	p<-ggplot(df) +
  		geom_bar(aes(x=trait_name, y=frac_h2_med),stat="identity",position=position_dodge(), fill="cyan4")+figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  		theme(legend.position="bottom") +
  		geom_errorbar( aes(x=trait_name, ymin=frac_h2_med_lb, ymax=frac_h2_med_ub), width=.2,colour="grey50", position=position_dodge()) +
  		labs(x="", y="Fraction h2 mediated", title=nm_variant_anno) +
  		theme(axis.text.x=element_blank())
  	return(p)
}




###################
# Command line args
####################
tglr_results_dir = args[1]
visualize_tglr_results_dir = args[2]



nm_variant_anno = "baselineLD_no_qtl"
summary_file <- paste0(tglr_results_dir, "cross_trait_med_h2_summary_susie_inf_pmces_", nm_variant_anno, ".txt")
df <- read.table(summary_file, header=TRUE, sep="\t")
pp <- make_se_barplot_showing_mediated_heritabilities(df, nm_variant_anno)
output_file <- paste0(visualize_tglr_results_dir, "frac_h2_med_per_trait_", nm_variant_anno, ".pdf")
ggsave(pp, file=output_file, width=7.2, height=5.5, units="in")
print(output_file)


nm_variant_anno = "genotype_intercept"
summary_file <- paste0(tglr_results_dir, "cross_trait_med_h2_summary_susie_inf_pmces_", nm_variant_anno, ".txt")
df <- read.table(summary_file, header=TRUE, sep="\t")
pp <- make_se_barplot_showing_mediated_heritabilities(df, nm_variant_anno)
output_file <- paste0(visualize_tglr_results_dir, "frac_h2_med_per_trait_", nm_variant_anno, ".pdf")
ggsave(pp, file=output_file, width=7.2, height=5.5, units="in")
print(output_file)
