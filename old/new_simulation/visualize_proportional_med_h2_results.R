args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(RColorBrewer)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}





make_se_barplot_for_proportional_mediation_estimates <- function(summary_file) {
	df <- read.table(summary_file, header=TRUE, sep="\t")
	df2 = df[as.character(df$heritability_type)=="med",]


	df2$mean_cis_h2 = factor(df2$mean_cis_h2, levels=c("0.01", "0.05", "0.1"))

	df2$method_name = factor(df2$method_name, levels=c("ss_reml_2_step", "ss_reml", "ss_bayes_lmm"))

	p<-ggplot(data=df2, aes(x=mean_cis_h2, y=est_h2, fill=method_name)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="expr cis h2", y="medated h2", fill="") +
  		geom_hline(yintercept=0.1, linetype=2) 
  	return(p)


}



make_se_barplot_for_proportional_mediation_estimates2 <- function(summary_file) {
	df <- read.table(summary_file, header=TRUE, sep="\t")
	df2 = df[as.character(df$heritability_type)=="med",]


	df2$mean_cis_h2 = factor(df2$mean_cis_h2, levels=c("0.01", "0.05", "0.1"))

	df2 = df2[as.character(df2$method_name) != "ss_reml",]

	df2$method_name = factor(df2$method_name, levels=c("ss_reml_2_step", "ss_bayes_lmm"))

	p<-ggplot(data=df2, aes(x=mean_cis_h2, y=est_h2, fill=method_name)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="expr cis h2", y="medated h2", fill="") +
  		geom_hline(yintercept=0.1, linetype=2) 
  	return(p)


}












med_h2_simulation_dir <- args[1]




#####################
# Make SE barplot showing estimates in single causal gene simulatin
#####################
summary_file <- paste0(med_h2_simulation_dir, "mediation_h2_summary.txt")
pp <- make_se_barplot_for_proportional_mediation_estimates(summary_file)
output_file <- paste0(med_h2_simulation_dir, "proportional_med_h2_estimates_se_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=5.5, units="in")


#####################
# Make SE barplot showing estimates in single causal gene simulatin
#####################
summary_file <- paste0(med_h2_simulation_dir, "mediation_h2_summary.txt")
pp <- make_se_barplot_for_proportional_mediation_estimates2(summary_file)
output_file <- paste0(med_h2_simulation_dir, "proportional_med_h2_estimates_se_barplot2.pdf")
ggsave(pp, file=output_file, width=7.2, height=3.8, units="in")