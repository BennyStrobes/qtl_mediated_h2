args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(dplyr)
library(reshape)
library(stringr)
library(RColorBrewer)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


make_power_med_h2_se_barplot <- function(df) {
 	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,200,300,500,1000, "Inf"))
	p<-ggplot(data=df, aes(x=eqtl_sample_size, y=power)) +
  		geom_bar(fill="skyblue",stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9), colour="orange")  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Power") +
  		theme(legend.position="bottom") 
  	return(p)
}

make_med_h2_per_eqtl_category_se_barplot <- function(summary_file, qtl_method) {
	df <- read.table(summary_file, header=TRUE, sep="\t")
	print(head(df))
	df$eqtl_category = as.character(df$eqtl_category)
	arr <- c()
	for (ii in 1:length(df$eqtl_category)) {
		df$eqtl_category[ii] = paste0(qtl_method, "_", ii)
		arr <- c(arr, paste0(qtl_method, "_", ii))
	}
	df$eqtl_category = factor(df$eqtl_category, levels=arr)

	p<-ggplot(data=df, aes(x=eqtl_category, y=mean)) +
  		geom_bar(fill="skyblue",stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=lb, ymax=ub), width=.4, position=position_dodge(.9), colour="orange")  +
  		figure_theme() +
  		labs(x="Gene Category", y="Normalized Mediated Heritability", title=qtl_method) +
  		theme(legend.position="bottom")+
  		 theme(axis.text.x = element_text(angle = 90))
  	return(p)
}


##################
# Load in data
##################
tglr_results_dir = args[1]
viz_results_dir = args[2]

qtl_method="distance_stratification"
result_summary_file = paste0(tglr_results_dir, "normalize_per_gene_heritability_meta_analyzed_cross_traits_", qtl_method, ".txt")
pp <- make_med_h2_per_eqtl_category_se_barplot(result_summary_file, qtl_method)
output_file <- paste0(viz_results_dir, "make_med_h2_per_eqtl_category_se_barplot_", qtl_method, ".pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")


qtl_method="qtl_rank_stratification"
result_summary_file = paste0(tglr_results_dir, "normalize_per_gene_heritability_meta_analyzed_cross_traits_", qtl_method, ".txt")
pp <- make_med_h2_per_eqtl_category_se_barplot(result_summary_file, qtl_method)
output_file <- paste0(viz_results_dir, "make_med_h2_per_eqtl_category_se_barplot_", qtl_method, ".pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")

