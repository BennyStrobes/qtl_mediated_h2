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

make_med_h2_se_barplot <- function(df) {
	df$eqtl_ss = factor(df$eqtl_ss)
	df$method = factor(df$method,levels=c("true_gene","blup_pred_gene", "pred_gene", "ridge_regr_pred_gene", "marginal_gene"))
	p<-ggplot(data=df, aes(x=eqtl_ss, y=med_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=med_h2_lb, ymax=med_h2_ub), width=.4, position=position_dodge(.9), colour="black")  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="h2-med") +
  		theme(legend.position="right") +
  		geom_hline(yintercept=.03, linetype=2)
  	return(p)
}




#######################
# Command line args
########################
visualize_med_h2_results_dir = args[1]






med_h2_summary_file <- paste0(visualize_med_h2_results_dir, "med_h2_summary.txt")
df <- read.table(med_h2_summary_file,header=TRUE,sep="\t")


pp <- make_med_h2_se_barplot(df)
output_file <- paste0(visualize_med_h2_results_dir, "med_h2_summary.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")







