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
	df$eqtl_ss = factor(df$eqtl_ss)
	df$method = factor(df$method,levels=c("no_attenuation_correction","attenuation_correction"))
	p<-ggplot(data=df, aes(x=eqtl_ss, y=med_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=med_h2_lb, ymax=med_h2_ub), width=.4, position=position_dodge(.9), colour="black")  +
  		figure_theme() +
  		labs(x="eQTL sample size", y="h2-med", fill="") +
  		theme(legend.position="bottom") +
  		geom_hline(yintercept=.09, linetype=2)
  	return(p)
}

make_med_h2_noise_diff_se_barplot2 <- function(df) {
	df$eqtl_ss = factor(df$eqtl_ss)
	p<-ggplot(data=df, aes(x=eqtl_ss, y=med_h2)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=med_h2_lb, ymax=med_h2_ub), width=.4, position=position_dodge(.9), colour="black")  +
  		figure_theme() +
  		labs(x="eQTL sample size", y="h2-med") +
  		geom_hline(yintercept=.09, linetype=2)
  	return(p)
}

make_med_h2_se_barplot <- function(df) {
	df$eqtl_ss = factor(df$eqtl_ss)
	p<-ggplot(data=df, aes(x=eqtl_ss, y=med_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=med_h2_lb, ymax=med_h2_ub), width=.4, position=position_dodge(.9), colour="black")  +
  		figure_theme() +
  		labs(x="eQTL sample size", y="h2-med", fill="") +
  		theme(legend.position="bottom") +
  		geom_hline(yintercept=.09, linetype=2)
  	return(p)
}

make_med_h2_se_barplot_with_limits <- function(df) {
	df$eqtl_ss = factor(df$eqtl_ss)
	p<-ggplot(data=df, aes(x=eqtl_ss, y=med_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=med_h2_lb, ymax=med_h2_ub), width=.4, position=position_dodge(.9), colour="black")  +
  		figure_theme() +
  		labs(x="eQTL sample size", y="h2-med", fill="") +
  		ylim(0,.2) +
  		theme(legend.position="bottom") +
  		geom_hline(yintercept=.09, linetype=2)
  	return(p)
}


make_noise_est_se_barplot <- function(df) {
	df$eqtl_ss = factor(df$eqtl_ss)
	df$noise_type = factor(df$noise_type,levels=c("true_noise","est_noise"))
	p<-ggplot(data=df, aes(x=eqtl_ss, y=noise_mean, fill=noise_type)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=noise_lb, ymax=noise_ub), width=.4, position=position_dodge(.9), colour="black")  +
  		figure_theme() +
  		labs(x="eQTL sample size", y="Average noise", fill="") +
  		theme(legend.position="bottom") 
  	return(p)
}

make_med_h2_est_w_w_o_atten <- function(df) {
	df$eqtl_ss = factor(df$eqtl_ss)
	df$method = factor(df$method,levels=c("no_attenuation_correction","attenuation_correction"))
	p<-ggplot(data=df, aes(x=eqtl_ss, y=med_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=med_h2_lb, ymax=med_h2_ub), width=.4, position=position_dodge(.9), colour="black")  +
  		figure_theme() +
  		labs(x="eQTL sample size", y="h2-med", fill="") +
  		theme(legend.position="bottom") +
  		geom_hline(yintercept=.09, linetype=2)
  	return(p)
}



#######################
# Command line args
########################
visualize_med_h2_results_dir = args[1]


med_h2_summary_file <- paste0(visualize_med_h2_results_dir, "med_h2_summary.txt")
df_med_h2 <- read.table(med_h2_summary_file, header=TRUE, sep="\t")
df_med_h2 = df_med_h2[as.character(df_med_h2$method) != "med_h2_est_prediction_gene_ld_scores_blup",]
df_med_h2$method <- recode(df_med_h2$method, med_h2_sim_gene_ld_scores="simulated", med_h2_est_gene_ld_scores="estimated", med_h2_est_gene_ld_scores_dissattenuated_true_noise="dissatten", med_h2_est_gene_ld_scores_dissattenuated_est_noise="dissatten_e_noise", med_h2_est_posterior_gene_ld_scores_known_gene_h2_known_residvar="posterior", med_h2_est_posterior_gene_ld_scores_blup="posterior_blup")
df_med_h2$method = factor(df_med_h2$method, levels=c("simulated", "estimated", "dissatten", "dissatten_e_noise", "posterior", "posterior_blup"))
pp <- make_med_h2_se_barplot(df_med_h2)
output_file <- paste0(visualize_med_h2_results_dir, "multi_eqtl_med_h2_se_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")


pp <- make_med_h2_se_barplot_with_limits(df_med_h2)
output_file <- paste0(visualize_med_h2_results_dir, "multi_eqtl_w_axis_limits_med_h2_se_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")

df_med_h2 = df_med_h2[as.character(df_med_h2$method) != "dissatten",]
df_med_h2 = df_med_h2[as.character(df_med_h2$method) != "dissatten_e_noise",]
pp <- make_med_h2_se_barplot(df_med_h2)
output_file <- paste0(visualize_med_h2_results_dir, "multi_eqtl_no_dissatten_med_h2_se_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")





# Noise estimate
if (FALSE) {
noise_estimation_summary_file <- paste0(visualize_med_h2_results_dir, "noise_estimates.txt")
df_noise <- read.table(noise_estimation_summary_file,header=TRUE,sep="\t")
pp <- make_noise_est_se_barplot(df_noise)
output_file <- paste0(visualize_med_h2_results_dir, "noise_est_se_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")
}

if (FALSE) {
# Attenuation correction
med_h2_summary_file <- paste0(visualize_med_h2_results_dir, "med_h2_summary_est_noise_correction.txt")
df_med_h2 <- read.table(med_h2_summary_file,header=TRUE,sep="\t")
pp <- make_med_h2_est_w_w_o_atten(df_med_h2)
output_file <- paste0(visualize_med_h2_results_dir, "med_h2_w_wo_atten_se_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")
}


if (FALSE) {

# Plot eqtl ss noise analysis
med_h2_noise_diff_summary_file <- paste0(visualize_med_h2_results_dir, "med_h2_summary_known_noise_correction.txt")
df_noise <- read.table(med_h2_noise_diff_summary_file,header=TRUE,sep="\t")


pp <- make_med_h2_noise_diff_se_barplot(df_noise)
output_file <- paste0(visualize_med_h2_results_dir, "med_h2_eqtl_ss_noise_summary.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")
print(output_file)
}

if (FALSE) {
df2 = df_noise[as.character(df_noise$method)=="no_attenuation_correction",]

pp <- make_med_h2_noise_diff_se_barplot2(df2)
output_file <- paste0(visualize_med_h2_results_dir, "med_h2_eqtl_ss_noise_summary2.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")
print(output_file)
}







