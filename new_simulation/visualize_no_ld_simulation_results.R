args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(RColorBrewer)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


make_se_barplot_for_single_causal_gene_simulation <- function(input_stem, eqtl_sss) {
	ss_arr <- c()
	mean_arr <- c()
	lb_arr <- c()
	ub_arr <- c()


	for (ss_iter in 1:length(eqtl_sss)) {
		ss <- eqtl_sss[ss_iter]
		filer <- paste0(input_stem, ss, "_param_summary.txt")
		df <- read.table(filer, header=TRUE, sep="\t")
		df = df[as.character(df$method) == "scipy_odr_effect_est_stringent",]
		ss_arr <- c(ss_arr, ss)
		mean_arr <- c(mean_arr, df$est[1])
		lb_arr <- c(lb_arr, df$est_95_lb[1])
		ub_arr <- c(ub_arr, df$est_95_ub[1])
	}

	df <- data.frame(alpha=mean_arr, eQTL_ss=ss_arr, alpha_lb=lb_arr, alpha_ub=ub_arr)

	df$eQTL_ss <- factor(df$eQTL_ss, levels=c("300","1000"))

	p<-ggplot(data=df, aes(x=eQTL_ss, y=alpha)) +
  		geom_bar(fill="skyblue",stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=alpha_lb, ymax=alpha_ub), width=.4, position=position_dodge(.9), colour="orange")  +
  		figure_theme() +
  		labs(x="eQTL Sample size", y="gene-trait effect") +
  		theme(legend.position="bottom")  +
  		geom_hline(yintercept=-1.41, linetype=2)

  	return(p)


}


make_se_barplot_for_sparse_simulation <- function(input_file) {
	df <- read.table(input_file,header=TRUE,sep="\t")

	df_med = df[as.character(df$effect_type) == "mediated",]
	df_med$eqtl_ss = factor(df_med$eqtl_ss, levels=c("300"))
	df_med$eqtl_ss = factor(df_med$eqtl_ss, levels=c("300"))
	df_med$estimation_method = factor(df_med$estimation_method, levels=c("posterior_distr", "pmces"))

	p_med<-ggplot(data=df_med, aes(x=eqtl_ss, y=bias, fill=estimation_method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=bias_lb, ymax=bias_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c("skyblue", "grey"))+
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Bias (mediated h2)", fill="", title="Mediated h2") +
  		geom_hline(yintercept=0.0) 

	df_nm = df[as.character(df$effect_type) == "non-mediated",]
	df_nm$eqtl_ss = factor(df_nm$eqtl_ss, levels=c("300"))
	df_nm$eqtl_ss = factor(df_nm$eqtl_ss, levels=c("300"))
	df_nm$estimation_method = factor(df_nm$estimation_method, levels=c("posterior_distr", "pmces"))

	p_nm<-ggplot(data=df_nm, aes(x=eqtl_ss, y=bias, fill=estimation_method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=bias_lb, ymax=bias_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c("skyblue", "grey"))+
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Bias (non-mediated h2)", fill="", title="Non-mediated h2") +
  		geom_hline(yintercept=0.0) 

  	p <- plot_grid(p_nm, p_med, ncol=1)

  	return(p)

}




make_se_barplot_for_polygenic_simulation <- function(input_file) {
	df <- read.table(input_file,header=TRUE,sep="\t")

	df_med = df[as.character(df$effect_type) == "mediated",]
	df_med$eqtl_ss = factor(df_med$eqtl_ss, levels=c("100","300"))
	df_med$eqtl_ss = factor(df_med$eqtl_ss, levels=c("100","300"))
	df_med$estimation_method = factor(df_med$estimation_method, levels=c("posterior_distr", "pmces"))

	p_med<-ggplot(data=df_med, aes(x=eqtl_ss, y=bias, fill=estimation_method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=bias_lb, ymax=bias_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c("skyblue", "grey"))+
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Bias (mediated h2)", fill="", title="Mediated h2") +
  		geom_hline(yintercept=0.0) 

	df_nm = df[as.character(df$effect_type) == "non-mediated",]
	df_nm$eqtl_ss = factor(df_nm$eqtl_ss, levels=c("100","300"))
	df_nm$eqtl_ss = factor(df_nm$eqtl_ss, levels=c("100","300"))
	df_nm$estimation_method = factor(df_nm$estimation_method, levels=c("posterior_distr", "pmces"))

	p_nm<-ggplot(data=df_nm, aes(x=eqtl_ss, y=bias, fill=estimation_method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=bias_lb, ymax=bias_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c("skyblue", "grey"))+
  		figure_theme() +
  		labs(x="eQTL Sample size", y="Bias (non-mediated h2)", fill="", title="Non-mediated h2") +
  		geom_hline(yintercept=0.0) 

  	p <- plot_grid(p_nm, p_med, ncol=1)

  	return(p)

}




make_se_barplot_for_three_causal_gene_simulation <- function(input_stem, eqtl_sss) {
	ss_arr <- c()
	mean_arr <- c()
	lb_arr <- c()
	ub_arr <- c()
	gene_arr <- c()


	for (ss_iter in 1:length(eqtl_sss)) {
		ss <- eqtl_sss[ss_iter]
		filer <- paste0(input_stem, ss, "_sim_corr_0.0_param_summary.txt")
		df <- read.table(filer, header=TRUE, sep="\t")
		df = df[as.character(df$method) == "scipy_odr_effect_est_stringent",]
		ss_arr <- c(ss_arr, ss)
		mean_arr <- c(mean_arr, df$est[1])
		lb_arr <- c(lb_arr, df$est_95_lb[1])
		ub_arr <- c(ub_arr, df$est_95_ub[1])
		gene_arr <- c(gene_arr, "gene1")

		ss_arr <- c(ss_arr, ss)
		mean_arr <- c(mean_arr, df$est[2])
		lb_arr <- c(lb_arr, df$est_95_lb[2])
		ub_arr <- c(ub_arr, df$est_95_ub[2])
		gene_arr <- c(gene_arr, "gene2")
		
		ss_arr <- c(ss_arr, ss)
		mean_arr <- c(mean_arr, df$est[3])
		lb_arr <- c(lb_arr, df$est_95_lb[3])
		ub_arr <- c(ub_arr, df$est_95_ub[3])
		gene_arr <- c(gene_arr, "gene3")


	}

	df <- data.frame(alpha=mean_arr*-1, eQTL_ss=ss_arr, alpha_ub=lb_arr*-1, alpha_lb=ub_arr*-1, gene=factor(gene_arr, levels=c("gene1", "gene2", "gene3")))


	df$eQTL_ss <- factor(df$eQTL_ss, levels=c("300","1000"))



	p<-ggplot(data=df, aes(x=eQTL_ss, y=alpha, fill=gene)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=alpha_lb, ymax=alpha_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL Sample size", y="gene-trait effect", fill="") +
  		geom_hline(yintercept=1.41, linetype=2) +
  		geom_hline(yintercept=0.0) +
  		theme(legend.position="none")


  	return(p)


}





no_ld_simulation_res_dir = args[1]


#####################
# Make SE barplot showing estimates in single causal gene simulatin
#####################
if (FALSE) {
input_file = paste0(no_ld_simulation_res_dir, "organized_results_polygenic_simulation.txt")
pp <- make_se_barplot_for_polygenic_simulation(input_file)
output_file <- paste0(no_ld_simulation_res_dir, "polygenic_sim_med_h2_bias.pdf")
ggsave(pp, file=output_file, width=7.2, height=5.5, units="in")
}

#####################
# Make SE barplot showing estimates in single causal gene simulatin
#####################
input_file = paste0(no_ld_simulation_res_dir, "organized_results_sparse_simulation.txt")
pp <- make_se_barplot_for_sparse_simulation(input_file)
output_file <- paste0(no_ld_simulation_res_dir, "sparse_sim_med_h2_bias.pdf")
ggsave(pp, file=output_file, width=7.2, height=5.5, units="in")


if (FALSE) {
#####################
# Make SE barplot showing estimates in single causal gene simulatin
#####################
input_stem = paste0(no_ld_simulation_res_dir, "merged_no_ld_no_nm_variants_1_caus_tiss_eqtlss_")
eqtl_sss = c("300", "1000")
pp <- make_se_barplot_for_single_causal_gene_simulation(input_stem, eqtl_sss)
output_file <- paste0(no_ld_simulation_res_dir, "one_caus_gene_gene_trait_effect_est.pdf")
ggsave(pp, file=output_file, width=7.2, height=3.5, units="in")

#####################
# Make SE barplot showing estimates in 3 causal gene simulatin
#####################
input_stem = paste0(no_ld_simulation_res_dir, "merged_no_ld_no_nm_variants_3_caus_tiss_eqtlss_")
eqtl_sss = c("300", "1000")
pp <- make_se_barplot_for_three_causal_gene_simulation(input_stem, eqtl_sss)
output_file <- paste0(no_ld_simulation_res_dir, "three_caus_gene_gene_trait_effect_est.pdf")
ggsave(pp, file=output_file, width=7.2, height=3.5, units="in")
}
