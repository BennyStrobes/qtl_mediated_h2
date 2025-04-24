args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
options(warn=1)






figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}





make_simulation_bias_barplot <- function(mediated_h2_results, simulation_type, sim_identifier) { 
	# Load in Data
	filer <- paste0(mediated_h2_results, "single_gene_sparse_eqtl_", simulation_type, "_", sim_identifier, "_eqtl_snps.txt")
	df1 <- read.table(filer, header=TRUE,sep="\t")
	df1 <- df1[as.character(df1$sim_identifier) == sim_identifier,]




	method_arr <- c()
	prime_arr <- c()
	est_arr <- c()
	est_lb_arr <- c()
	est_ub_arr <- c()

	# MESC
	tmp <- df1[as.character(df1$method) == "mesc",]
	method_arr <- c(method_arr, "mesc")
	prime_arr <- c(prime_arr, "causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])
	# MESC prime
	tmp <- df1[as.character(df1$method) == "mesc_prime",]
	method_arr <- c(method_arr, "mesc")
	prime_arr <- c(prime_arr, "non_causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])
	# MESC_w_gene_window
	tmp <- df1[as.character(df1$method) == "mesc_w_gene_window",]
	method_arr <- c(method_arr, "mesc_w_gene_window")
	prime_arr <- c(prime_arr, "causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])
	# MESC_prime_w_gene_window
	tmp <- df1[as.character(df1$method) == "mesc_prime_w_gene_window",]
	method_arr <- c(method_arr, "mesc_w_gene_window")
	prime_arr <- c(prime_arr, "non_causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])
	# GECS
	tmp <- df1[as.character(df1$method) == "gecs_2_thresh",]
	method_arr <- c(method_arr, "gecs")
	prime_arr <- c(prime_arr, "causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])
	# MESC_prime_w_gene_window
	tmp <- df1[as.character(df1$method) == "gecs_prime_2_thresh",]
	method_arr <- c(method_arr, "gecs")
	prime_arr <- c(prime_arr, "non_causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])
	


	df <- data.frame(method=factor(method_arr, levels=c("mesc", "mesc_w_gene_window", "gecs")), tissue=factor(prime_arr), alpha_sq=est_arr, alpha_sq_ub=est_ub_arr, alpha_sq_lb=est_lb_arr)

	tmp <- df1[as.character(df1$method) == "corr_marginal_sq_effects",]
	corr_marginal_sq_effects <-tmp$est[1]
	tmp <- df1[as.character(df1$method) == "corr_causal_sq_effects",]
	corr_causal_sq_effects <-tmp$est[1]
	tmp <- df1[as.character(df1$method) == "sq_corr_causal_effects",]
	sq_corr_causal_effects <-tmp$est[1]	

	pp<-ggplot(data=df, aes(x=method, y=alpha_sq, fill=tissue)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=alpha_sq_lb, ymax=alpha_sq_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="Method", y="alpha_sq", fill="", title=paste0(simulation_type, " / ", sim_identifier)) +
  		geom_hline(yintercept=1.0, linetype=2) + 
  		geom_hline(yintercept=corr_marginal_sq_effects, linetype=2, color="orange") + 
  		geom_hline(yintercept=corr_causal_sq_effects, linetype=2, color='green') + 
  		geom_hline(yintercept=sq_corr_causal_effects, linetype=2, color='purple') +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  	return(pp)

}



make_simulation_bias_barplot_over_gecs_range <- function(mediated_h2_results, simulation_type, sim_identifier) { 
	# Load in Data
	filer <- paste0(mediated_h2_results, "single_gene_sparse_eqtl_", simulation_type, "_", sim_identifier, "_eqtl_snps.txt")
	df1 <- read.table(filer, header=TRUE,sep="\t")
	df1 <- df1[as.character(df1$sim_identifier) == sim_identifier,]




	method_arr <- c()
	prime_arr <- c()
	est_arr <- c()
	est_lb_arr <- c()
	est_ub_arr <- c()

	if (FALSE) {
	# MESC
	tmp <- df1[as.character(df1$method) == "mesc",]
	method_arr <- c(method_arr, "mesc")
	prime_arr <- c(prime_arr, "causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])
	# MESC prime
	tmp <- df1[as.character(df1$method) == "mesc_prime",]
	method_arr <- c(method_arr, "mesc")
	prime_arr <- c(prime_arr, "non_causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])
	# MESC_w_gene_window
	tmp <- df1[as.character(df1$method) == "mesc_w_gene_window",]
	method_arr <- c(method_arr, "mesc_w_gene_window")
	prime_arr <- c(prime_arr, "causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])
	# MESC_prime_w_gene_window
	tmp <- df1[as.character(df1$method) == "mesc_prime_w_gene_window",]
	method_arr <- c(method_arr, "mesc_w_gene_window")
	prime_arr <- c(prime_arr, "non_causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])
	}

	# GECS
	tmp <- df1[as.character(df1$method) == "gecs_all_snp_pairs",]
	method_arr <- c(method_arr, "gecs_all")
	prime_arr <- c(prime_arr, "causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])
	# MESC_prime_w_gene_window
	tmp <- df1[as.character(df1$method) == "gecs_prime_all_snp_pairs",]
	method_arr <- c(method_arr, "gecs_all")
	prime_arr <- c(prime_arr, "non_causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])
	# GECS
	tmp <- df1[as.character(df1$method) == "gecs_10_thresh",]
	method_arr <- c(method_arr, "gecs_10")
	prime_arr <- c(prime_arr, "causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])
	# MESC_prime_w_gene_window
	tmp <- df1[as.character(df1$method) == "gecs_prime_10_thresh",]
	method_arr <- c(method_arr, "gecs_10")
	prime_arr <- c(prime_arr, "non_causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])

	# GECS
	tmp <- df1[as.character(df1$method) == "gecs_2_thresh",]
	method_arr <- c(method_arr, "gecs_2")
	prime_arr <- c(prime_arr, "causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])
	# MESC_prime_w_gene_window
	tmp <- df1[as.character(df1$method) == "gecs_prime_2_thresh",]
	method_arr <- c(method_arr, "gecs_2")
	prime_arr <- c(prime_arr, "non_causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])

	# GECS
	tmp <- df1[as.character(df1$method) == "gecs_.5_thres_thresh",]
	method_arr <- c(method_arr, "gecs_.5")
	prime_arr <- c(prime_arr, "causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])
	# MESC_prime_w_gene_window
	tmp <- df1[as.character(df1$method) == "gecs_prime_.5_thresh",]
	method_arr <- c(method_arr, "gecs_.5")
	prime_arr <- c(prime_arr, "non_causal_tissue")
	est_arr <- c(est_arr, tmp$est[1])
	est_lb_arr <- c(est_lb_arr, tmp$est_lb[1])
	est_ub_arr <- c(est_ub_arr, tmp$est_ub[1])

	


	df <- data.frame(method=factor(method_arr, levels=c("gecs_all", "gecs_10", "gecs_2", "gecs_.5")), tissue=factor(prime_arr), alpha_sq=est_arr, alpha_sq_ub=est_ub_arr, alpha_sq_lb=est_lb_arr)

	tmp <- df1[as.character(df1$method) == "corr_marginal_sq_effects",]
	corr_marginal_sq_effects <-tmp$est[1]
	tmp <- df1[as.character(df1$method) == "corr_causal_sq_effects",]
	corr_causal_sq_effects <-tmp$est[1]
	tmp <- df1[as.character(df1$method) == "sq_corr_causal_effects",]
	sq_corr_causal_effects <-tmp$est[1]	

	pp<-ggplot(data=df, aes(x=method, y=alpha_sq, fill=tissue)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=alpha_sq_lb, ymax=alpha_sq_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="Method", y="alpha_sq", fill="", title=paste0(simulation_type, " / ", sim_identifier)) +
  		geom_hline(yintercept=1.0, linetype=2) + 
  		geom_hline(yintercept=corr_marginal_sq_effects, linetype=2, color="orange") + 
  		geom_hline(yintercept=corr_causal_sq_effects, linetype=2, color='green') + 
  		geom_hline(yintercept=sq_corr_causal_effects, linetype=2, color='purple') +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  	return(pp)

}




mediated_h2_results = args[1]
viz_dir = args[2]



print(mediated_h2_results)
print(viz_dir)



# Make plot
simulation_type <- "unshared_causal_snps_extra_nm_var_in_gene"
sim_identifier <- "25"
pp <- make_simulation_bias_barplot_over_gecs_range(mediated_h2_results, simulation_type, sim_identifier)
output_file <- paste0(viz_dir, "simulation_unshared_causal_snps_extra_nm_var_in_gene_", sim_identifier, "_est_se_barplot_range_of_gecs_thresholds.pdf")
ggsave(pp, file=output_file, width=7.2, height=3.5, units="in")
print(output_file)


if (FALSE) {
# Make plot
simulation_type <- "unshared_causal_snps_extra_nm_var_in_gene"
sim_identifier <- "25"
pp <- make_simulation_bias_barplot(mediated_h2_results, simulation_type, sim_identifier)
output_file <- paste0(viz_dir, "simulation_unshared_causal_snps_extra_nm_var_in_gene_", sim_identifier, "_est_se_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=3.5, units="in")
print(output_file)

sim_identifier <- "5"
pp <- make_simulation_bias_barplot(mediated_h2_results, simulation_type, sim_identifier)
output_file <- paste0(viz_dir, "simulation_unshared_causal_snps_extra_nm_var_in_gene_", sim_identifier, "_est_se_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=3.5, units="in")
print(output_file)

sim_identifier <- "2"
pp <- make_simulation_bias_barplot(mediated_h2_results, simulation_type, sim_identifier)
output_file <- paste0(viz_dir, "simulation_unshared_causal_snps_extra_nm_var_in_gene_", sim_identifier, "_est_se_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=3.5, units="in")
print(output_file)
}











if (FALSE) {
# Make plot
simulation_type <- "shared_causal_snps"
sim_identifier <- "0.2"
pp_2 <- make_simulation_bias_barplot(simulation_type, sim_identifier)
sim_identifier <- "0.4"
pp_4 <- make_simulation_bias_barplot(simulation_type, sim_identifier)
sim_identifier <- "0.6"
pp_6 <- make_simulation_bias_barplot(simulation_type, sim_identifier)
sim_identifier <- "0.8"
pp_8 <- make_simulation_bias_barplot(simulation_type, sim_identifier)
joint_plot <- plot_grid(pp_2, pp_4, pp_6, ncol=1)

output_file <- paste0("simulation_shared_causal_snps_est_se_barplot.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=6.5, units="in")


# Make plot
simulation_type <- "unshared_causal_snps"
sim_identifier <- "25"
pp_2 <- make_simulation_bias_barplot(simulation_type, sim_identifier)
sim_identifier <- "100"
pp_4 <- make_simulation_bias_barplot(simulation_type, sim_identifier)
sim_identifier <- "1000"
pp_6 <- make_simulation_bias_barplot(simulation_type, sim_identifier)
joint_plot <- plot_grid(pp_2, pp_4, pp_6, ncol=1)

output_file <- paste0("simulation_unshared_causal_snps_est_se_barplot.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=6.5, units="in")

}



if (FALSE) {
# Make plot
simulation_type <- "unshared_causal_snps"

sim_identifier <- "50"
pp <- make_simulation_bias_barplot(simulation_type, sim_identifier)

output_file <- paste0("viz/simulation_unshared_causal_snps_est_se_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=3.5, units="in")
print(output_file)


# Make plot
simulation_type <- "unshared_causal_snps_less_nm_var_in_gene"

sim_identifier <- "50"
pp <- make_simulation_bias_barplot(simulation_type, sim_identifier)

output_file <- paste0("viz/simulation_unshared_causal_snps_less_nm_var_in_gene_est_se_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=3.5, units="in")
print(output_file)

# Make plot
simulation_type <- "unshared_causal_snps_extra_nm_var_in_gene"
sim_identifier <- "50"
pp <- make_simulation_bias_barplot(simulation_type, sim_identifier)
output_file <- paste0("viz/simulation_unshared_causal_snps_extra_nm_var_in_gene_est_se_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=3.5, units="in")
print(output_file)


# Make plot
simulation_type <- "unshared_causal_snps_extra_nm_var_in_middle_of_gene"
sim_identifier <- "50"
pp <- make_simulation_bias_barplot(simulation_type, sim_identifier)
output_file <- paste0("viz/simulation_unshared_causal_snps_extra_nm_var_in_middle_of_gene_est_se_barplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=3.5, units="in")
print(output_file)
}

