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

make_med_h2_selection_diff_se_barplot <- function(df) {
	df$simulated_selection = factor(df$simulated_selection, levels=c("neutral", "negative"))
	df$method = factor(df$method,levels=c("ldsc_raw_ge","ldsc_standardize_ge"))
	p<-ggplot(data=df, aes(x=simulated_selection, y=med_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=med_h2_lb, ymax=med_h2_ub), width=.4, position=position_dodge(.9), colour="black")  +
  		figure_theme() +
  		labs(x="Simulated gene-trait selection", y="h2-med") +
  		theme(legend.position="right") +
  		geom_hline(yintercept=.09, linetype=2)
  	return(p)
}

make_med_h2_noise_diff_se_barplot <- function(df) {
	df$simulated_noise = factor(df$simulated_noise)
	df$method = factor(df$method,levels=c("no_correction","correction_1", "correction2"))
	p<-ggplot(data=df, aes(x=method, y=med_h2, fill=simulated_noise)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=med_h2_lb, ymax=med_h2_ub), width=.4, position=position_dodge(.9), colour="black")  +
  		figure_theme() +
  		labs(x="Attenuation correction method", y="h2-med") +
  		theme(legend.position="right") +
  		geom_hline(yintercept=.09, linetype=2)
  	return(p)
}

make_med_h2_noise_diff_se_barplot_only_uncorrected <- function(df) {
	df$simulated_noise = factor(df$simulated_noise)
	df = df[as.character(df$method) == "no_correction", ]
	p<-ggplot(data=df, aes(x=simulated_noise, y=med_h2)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=med_h2_lb, ymax=med_h2_ub), width=.4, position=position_dodge(.9), colour="black")  +
  		figure_theme() +
  		labs(x="Simulated noise", y="h2-med") +
  		theme(legend.position="right") +
  		geom_hline(yintercept=.09, linetype=2)
  	return(p)
}



#######################
# Command line args
########################
visualize_med_h2_results_dir = args[1]


# Plot simulated noise analysis
med_h2_noise_diff_summary_file <- paste0(visualize_med_h2_results_dir, "med_h2_summary_simulated_noise_correction.txt")
df_noise <- read.table(med_h2_noise_diff_summary_file,header=TRUE,sep="\t")

pp <- make_med_h2_noise_diff_se_barplot(df_noise)
output_file <- paste0(visualize_med_h2_results_dir, "med_h2_noise_diff_summary.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")


pp <- make_med_h2_noise_diff_se_barplot_only_uncorrected(df_noise)
output_file <- paste0(visualize_med_h2_results_dir, "med_h2_noise_diff_only_uncorrected_summary.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")



# Plot selection analysis
med_h2_selection_diff_summary_file <- paste0(visualize_med_h2_results_dir, "med_h2_summary_selection_differences.txt")
df_selection <- read.table(med_h2_selection_diff_summary_file,header=TRUE,sep="\t")

pp <- make_med_h2_selection_diff_se_barplot(df_selection)
output_file <- paste0(visualize_med_h2_results_dir, "med_h2_selection_diff_summary.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")







