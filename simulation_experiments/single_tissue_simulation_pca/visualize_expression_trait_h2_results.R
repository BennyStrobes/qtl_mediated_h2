args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(RColorBrewer)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}




make_se_barplot_stratefied_by_sample_size_and_h2_version <- function(df, methods) {
	df2 <- df[df$method_name %in% methods,]
	df2$eQTL_SS <- factor(df2$eQTL_SS)
	
	true_mean = mean(df2$true_mean)

	pp<-ggplot(data=df2, aes(x=eQTL_SS, y=est_mean, fill=h2_version)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_mean_lb, ymax=est_mean_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="") +
  		geom_hline(yintercept=true_mean, linetype=2) 

  	return(pp)

}


make_se_barplot_stratefied_by_sample_size_and_priors <- function(df, methods) {
	df2 <- df[df$method_name %in% methods,]
	df2$eQTL_SS <- factor(df2$eQTL_SS)

	df2 <- df2[df2$h2_version=="joint_h2",]

	true_mean = mean(df2$true_mean)

	df2$method_name = factor(df2$method_name, levels=methods)

	pp<-ggplot(data=df2, aes(x=eQTL_SS, y=est_mean, fill=method_name)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_mean_lb, ymax=est_mean_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="") +
  		geom_hline(yintercept=true_mean, linetype=2) 

  	return(pp)

}







#####################
# Command line args
#####################
results_dir = args[1]
viz_dir = args[2]


summary_file <- paste0(viz_dir, "expression_trait_h2_results_summary.txt")
df <- read.table(summary_file, header=TRUE)


methods <- c("gibbs_lmm_multivariate_ig_prior_1e-06")
pp_gibbs <- make_se_barplot_stratefied_by_sample_size_and_h2_version(df, methods)
output_file <- paste0(viz_dir, "se_barplot_h2_for_gibbs_sampler_stratefied_by_sample_size_and_h2_version.pdf")
#ggsave(pp_gibbs, file=output_file, width=7.2, height=4.5, units="in")

methods <- c("VI_lmm_multivariate_ig_prior_1e-06")
pp_vi <- make_se_barplot_stratefied_by_sample_size_and_h2_version(df, methods)
output_file <- paste0(viz_dir, "se_barplot_h2_for_VI_stratefied_by_sample_size_and_h2_version.pdf")
#ggsave(pp_vi, file=output_file, width=7.2, height=4.5, units="in")

methods <- c("gibbs_lmm_univariate_ig_prior_1e-06")
pp_gibbs_uni <- make_se_barplot_stratefied_by_sample_size_and_h2_version(df, methods)
output_file <- paste0(viz_dir, "se_barplot_h2_for_gibbs_sampler_stratefied_by_sample_size_and_h2_version.pdf")
#ggsave(pp_gibbs, file=output_file, width=7.2, height=4.5, units="in")

methods <- c("VI_lmm_univariate_ig_prior_1e-06")
pp_vi_uni <- make_se_barplot_stratefied_by_sample_size_and_h2_version(df, methods)
output_file <- paste0(viz_dir, "se_barplot_h2_for_VI_stratefied_by_sample_size_and_h2_version.pdf")
#ggsave(pp_vi, file=output_file, width=7.2, height=4.5, units="in")



joint_plot <- plot_grid(pp_gibbs + labs(title="Multivariate Gibbs"), pp_vi+ labs(title="Multivariate VI"), pp_gibbs_uni + labs(title="Univariate Gibbs"), pp_vi_uni + labs(title="Univariate VI"), ncol=2)
output_file <- paste0(viz_dir, "se_barplot_h2_for_joint_stratefied_by_sample_size_and_h2_version.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=5.5, units="in")




methods <- c("gibbs_lmm_multivariate_ig_prior_0.001", "gibbs_lmm_multivariate_ig_prior_1e-06", "gibbs_lmm_multivariate_ig_prior_1e-10", "gibbs_lmm_multivariate_ig_prior_0.0")
pp_gibbs <- make_se_barplot_stratefied_by_sample_size_and_priors(df, methods)
output_file <- paste0(viz_dir, "se_barplot_h2_for_gibbs_sampler_stratefied_by_sample_size_and_prior_specification.pdf")
#ggsave(pp_gibbs, file=output_file, width=7.2, height=4.5, units="in")

methods <- c("VI_lmm_multivariate_ig_prior_0.001", "VI_lmm_multivariate_ig_prior_1e-06", "VI_lmm_multivariate_ig_prior_1e-10", "VI_lmm_multivariate_ig_prior_0.0")
pp_vi <- make_se_barplot_stratefied_by_sample_size_and_priors(df, methods)
output_file <- paste0(viz_dir, "se_barplot_h2_for_VI_stratefied_by_sample_size_and_prior_specification.pdf")
#ggsave(pp_vi, file=output_file, width=7.2, height=4.5, units="in")

joint_plot <- plot_grid(pp_gibbs + labs(title="Multivariate Gibbs"), pp_vi+ labs(title="Multivariate VI"), ncol=1)
output_file <- paste0(viz_dir, "se_barplot_multivariate_h2_for_joint_stratefied_by_sample_size_and_prior_specification.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=5.5, units="in")





methods <- c("gibbs_lmm_univariate_ig_prior_0.001", "gibbs_lmm_univariate_ig_prior_1e-06", "gibbs_lmm_univariate_ig_prior_1e-10", "gibbs_lmm_univariate_ig_prior_0.0")
pp_gibbs2 <- make_se_barplot_stratefied_by_sample_size_and_priors(df, methods)
output_file <- paste0(viz_dir, "se_barplot_h2_for_gibbs_sampler_stratefied_by_sample_size_and_prior_specification.pdf")
#ggsave(pp_gibbs, file=output_file, width=7.2, height=4.5, units="in")

methods <- c("VI_lmm_univariate_ig_prior_0.001", "VI_lmm_univariate_ig_prior_1e-06", "VI_lmm_univariate_ig_prior_1e-10", "VI_lmm_univariate_ig_prior_0.0")
pp_vi2 <- make_se_barplot_stratefied_by_sample_size_and_priors(df, methods)
output_file <- paste0(viz_dir, "se_barplot_h2_for_VI_stratefied_by_sample_size_and_prior_specification.pdf")
#ggsave(pp_vi, file=output_file, width=7.2, height=4.5, units="in")

joint_plot <- plot_grid(pp_gibbs2 + labs(title="Univariate Gibbs"), pp_vi2+ labs(title="Univariate VI"), ncol=1)
output_file <- paste0(viz_dir, "se_barplot_univariate_h2_for_joint_stratefied_by_sample_size_and_prior_specification.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=5.5, units="in")


legender <- get_legend(pp_gibbs+guides(fill=guide_legend(nrow=2,byrow=TRUE)))
joint_plot <- plot_grid(pp_gibbs + labs(title="Multivariate Gibbs") + theme(legend.position="none"), pp_vi+ labs(title="Multivariate VI")+ theme(legend.position="none"), pp_gibbs2 + labs(title="Univariate Gibbs")+ theme(legend.position="none"), pp_vi2+ labs(title="Univariate VI")+ theme(legend.position="none"), ncol=2)
joint_plot2 <- plot_grid(joint_plot, legender,ncol=1, rel_heights=c(1,.2))
output_file <- paste0(viz_dir, "se_barplot_h2_for_joint_stratefied_by_sample_size_and_prior_specification.pdf")
ggsave(joint_plot2, file=output_file, width=7.2, height=5.0, units="in")

