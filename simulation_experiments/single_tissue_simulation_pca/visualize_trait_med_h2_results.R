args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(RColorBrewer)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}



poster_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=30), text = element_text(size=25),axis.text=element_text(size=30), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=30), legend.title = element_text(size=30)))
}




make_se_barplot_for_nm_h2_old <- function(df, method) {
	df2 <- df[as.character(df$method)==method,]
	df2 <- df2[as.character(df2$heritability_type)=="nm_h2"| as.character(df2$heritability_type)=="alt_nm_h2",]
	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,300,1000,10000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("nm_h2", "alt_nm_h2"))

	sim_h2 = df2$sim_h2[1]

	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2, fill=heritability_type)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0(method, " / non-mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)


}

make_se_barplot_for_med_h2_old <- function(df, method) {
	df2 <- df[as.character(df$method)==method,]
	df2 <- df2[as.character(df2$heritability_type)=="med_h2"| as.character(df2$heritability_type)=="alt_med_h2",]
	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,300,1000,10000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("med_h2", "alt_med_h2"))

	sim_h2 = df2$sim_h2[1]

	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2, fill=heritability_type)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0(method, " / mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)


}





make_se_barplot_for_nm_h2_cross_priors <- function(df, method) {
	df2 <- df[as.character(df$method)==method,]
	df2 <- df2[as.character(df2$heritability_type)=="nm_h2"| as.character(df2$heritability_type)=="alt_nm_h2",]
	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,300,1000,10000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("nm_h2", "alt_nm_h2"))

	sim_h2 = df2$sim_h2[1]

	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2, fill=heritability_type)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0(method, " / non-mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)


}




make_se_barplot_for_nm_h2 <- function(df, method) {
	df2 <- df[as.character(df$method)==method,]
	df2 <- df2[as.character(df2$heritability_type)=="alt_nm_h2",]
	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,300,1000,10000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("alt_nm_h2"))

	sim_h2 = df2$sim_h2[1]

	blue_colors=brewer.pal(n = 9, name = "Blues")

	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2)) +
  		geom_bar(stat="identity", position=position_dodge(), fill=blue_colors[7]) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="QTL sample size", y="heritability", fill="", title=paste0("non-mediated heritability")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)


}

make_se_barplot_for_nm_h2_joint_vs_marginal <- function(df, method) {
	df2 <- df[as.character(df$method)==method,]
	df2 <- df2[(as.character(df2$heritability_type)=="alt_nm_h2") | (as.character(df2$heritability_type)=="nm_h2"),]
	
	print(df2)
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,300,1000,10000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("nm_h2","alt_nm_h2"))

	sim_h2 = df2$sim_h2[1]

	red_colors=brewer.pal(n = 9, name = "Reds")

	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2, fill=heritability_type)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("non-mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)


}

make_se_barplot_for_med_h2_joint_vs_marginal <- function(df, method) {
	df2 <- df[as.character(df$method)==method,]
	df2 <- df2[(as.character(df2$heritability_type)=="alt_med_h2") | (as.character(df2$heritability_type)=="med_h2"),]
	
	print(df2)
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,300,1000,10000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("med_h2","alt_med_h2"))

	sim_h2 = df2$sim_h2[1]

	red_colors=brewer.pal(n = 9, name = "Reds")

	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2, fill=heritability_type)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)


}



make_se_barplot_for_nm_h2_joint_vs_mesc <- function(df, method, df_tglr) {
	df2 <- df[as.character(df$method)==method,]
	df2 <- df2[(as.character(df2$heritability_type)=="alt_nm_h2"),]

	df3 <- df_tglr[as.character(df_tglr$method) =="tglr_unscaled",]
	df3 <- df3[as.character(df3$heritability_type) =="nm_h2",]

	
	sample_size_arr <- c()
	method_arr <- c()
	h2_est_arr <- c()
	h2_lb_arr <- c()
	h2_ub_arr <- c()

	sample_size_arr <- c(sample_size_arr, df2$eqtl_sample_size)
	method_arr <- c(method_arr, rep("MPMH", length(df2$eqtl_sample_size)))
	h2_est_arr <- c(h2_est_arr, df2$est_h2)
	h2_lb_arr <- c(h2_lb_arr, df2$est_h2_lb)
	h2_ub_arr <- c(h2_ub_arr, df2$est_h2_ub)

	sample_size_arr <- c(sample_size_arr, df3$eqtl_sample_size)
	method_arr <- c(method_arr, rep("MESC", length(df3$eqtl_sample_size)))
	h2_est_arr <- c(h2_est_arr, df3$est_h2+.005)
	h2_lb_arr <- c(h2_lb_arr, df3$est_h2_lb+.005)
	h2_ub_arr <- c(h2_ub_arr, df3$est_h2_ub+.005)

	final_df <- data.frame(eqtl_sample_size=sample_size_arr, est_h2=h2_est_arr, method=method_arr, est_h2_lb=h2_lb_arr, est_h2_ub=h2_ub_arr)


	final_df$eqtl_sample_size <- factor(final_df$eqtl_sample_size, levels=c(100,300,1000,10000))

	final_df$method <- factor(final_df$method, levels=c("MESC", "MPMH"))

	sim_h2 = df2$sim_h2[1]

	red_colors=brewer.pal(n = 9, name = "Reds")

	pp<-ggplot(data=final_df, aes(x=eqtl_sample_size, y=est_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("non-mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)
}



make_se_barplot_for_med_h2_joint_vs_mesc <- function(df, method, df_tglr) {
	df2 <- df[as.character(df$method)==method,]
	df2 <- df2[(as.character(df2$heritability_type)=="alt_med_h2"),]

	df3 <- df_tglr[as.character(df_tglr$method) =="tglr_unscaled",]
	df3 <- df3[as.character(df3$heritability_type) =="med_h2",]

	
	sample_size_arr <- c()
	method_arr <- c()
	h2_est_arr <- c()
	h2_lb_arr <- c()
	h2_ub_arr <- c()

	sample_size_arr <- c(sample_size_arr, df2$eqtl_sample_size)
	method_arr <- c(method_arr, rep("MPMH", length(df2$eqtl_sample_size)))
	h2_est_arr <- c(h2_est_arr, df2$est_h2)
	h2_lb_arr <- c(h2_lb_arr, df2$est_h2_lb)
	h2_ub_arr <- c(h2_ub_arr, df2$est_h2_ub)

	sample_size_arr <- c(sample_size_arr, df3$eqtl_sample_size)
	method_arr <- c(method_arr, rep("MESC", length(df3$eqtl_sample_size)))
	h2_est_arr <- c(h2_est_arr, df3$est_h2-.005)
	h2_lb_arr <- c(h2_lb_arr, df3$est_h2_lb-.005)
	h2_ub_arr <- c(h2_ub_arr, df3$est_h2_ub-.005)

	final_df <- data.frame(eqtl_sample_size=sample_size_arr, est_h2=h2_est_arr, method=method_arr, est_h2_lb=h2_lb_arr, est_h2_ub=h2_ub_arr)


	final_df$eqtl_sample_size <- factor(final_df$eqtl_sample_size, levels=c(100,300,1000,10000))

	final_df$method <- factor(final_df$method, levels=c("MESC", "MPMH"))

	sim_h2 = df2$sim_h2[1]

	red_colors=brewer.pal(n = 9, name = "Reds")

	pp<-ggplot(data=final_df, aes(x=eqtl_sample_size, y=est_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c("grey", red_colors[7]))+
  		poster_figure_theme() +
  		labs(x="QTL sample size", y="Heritability", fill="", title=paste0("Mediated heritability")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)
}



make_se_barplot_for_med_h2 <- function(df, method) {
	df2 <- df[as.character(df$method)==method,]
	df2 <- df2[as.character(df2$heritability_type)=="alt_med_h2",]
	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,300,1000,10000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("alt_med_h2"))

	sim_h2 = df2$sim_h2[1]

	red_colors=brewer.pal(n = 9, name = "Reds")

	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2)) +
  		geom_bar(stat="identity", position=position_dodge(), fill=red_colors[7]) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="QTL sample size", y="heritability", fill="", title=paste0("mediated heritability")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)


}


make_se_barplot_for_total_h2 <- function(df, method) {
	df2 <- df[as.character(df$method)==method,]
	df2 <- df2[as.character(df2$heritability_type)=="total_h2",]
	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,300,1000,10000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("total_h2"))

	sim_h2 = df2$sim_h2[1]



	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2)) +
  		geom_bar(stat="identity", position=position_dodge(), fill='grey16') +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="QTL sample size", y="heritability", fill="", title=paste0("Total heritability")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)
}

make_se_barplot_for_total_h2_over_iterations <- function(visualize_trait_med_h2_dir, iter_names, method) {

	eqtl_sample_size_arr <- c()
	est_h2_arr <- c()
	est_h2_lb_arr <- c()
	est_h2_ub_arr <- c()
	iter_name_arr <- c()

	for (itera in 1:length(iter_names)) {
		iter_name <- iter_names[itera]

		input_file <- paste0(visualize_trait_med_h2_dir, "med_h2_results_summary_averaged", iter_name,".txt")
		tmp_df <-read.table(input_file, header=TRUE, sep="\t")

		tmp_df <- tmp_df[as.character(tmp_df$method)==method,]
		tmp_df <- tmp_df[as.character(tmp_df$heritability_type)=="total_h2",]

		eqtl_sample_size_arr <- c(eqtl_sample_size_arr, tmp_df$eqtl_sample_size)
		est_h2_arr <- c(est_h2_arr, tmp_df$est_h2)
		est_h2_lb_arr <- c(est_h2_lb_arr, tmp_df$est_h2_lb)
		est_h2_ub_arr <- c(est_h2_ub_arr, tmp_df$est_h2_ub)

		iter_name_arr <- c(iter_name_arr, rep(iter_name, length(tmp_df$est_h2_ub)))

		sim_h2 = tmp_df$sim_h2[1]
	}

	df <- data.frame(iteration=iter_name_arr,eqtl_sample_size=eqtl_sample_size_arr, est_h2=est_h2_arr, est_h2_lb=est_h2_lb_arr, est_h2_ub=est_h2_ub_arr)

	df$iteration = factor(df$iteration, levels=iter_names)
	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,300,1000,10000))

	red_colors=brewer.pal(n = 9, name = "Greys")


	pp<-ggplot(data=df, aes(x=iteration, fill=eqtl_sample_size, y=est_h2)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c(red_colors[3], red_colors[5], red_colors[7], red_colors[9]))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("total h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) +
  		theme(legend.position="top") +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


  	return(pp)

}

make_se_barplot_for_nm_h2_over_iterations <- function(visualize_trait_med_h2_dir, iter_names, method) {

	eqtl_sample_size_arr <- c()
	est_h2_arr <- c()
	est_h2_lb_arr <- c()
	est_h2_ub_arr <- c()
	iter_name_arr <- c()

	for (itera in 1:length(iter_names)) {
		iter_name <- iter_names[itera]

		input_file <- paste0(visualize_trait_med_h2_dir, "med_h2_results_summary_averaged", iter_name,".txt")
		tmp_df <-read.table(input_file, header=TRUE, sep="\t")

		tmp_df <- tmp_df[as.character(tmp_df$method)==method,]
		tmp_df <- tmp_df[as.character(tmp_df$heritability_type)=="alt_nm_h2",]

		eqtl_sample_size_arr <- c(eqtl_sample_size_arr, tmp_df$eqtl_sample_size)
		est_h2_arr <- c(est_h2_arr, tmp_df$est_h2)
		est_h2_lb_arr <- c(est_h2_lb_arr, tmp_df$est_h2_lb)
		est_h2_ub_arr <- c(est_h2_ub_arr, tmp_df$est_h2_ub)

		iter_name_arr <- c(iter_name_arr, rep(iter_name, length(tmp_df$est_h2_ub)))

		sim_h2 = tmp_df$sim_h2[1]
	}

	df <- data.frame(iteration=iter_name_arr,eqtl_sample_size=eqtl_sample_size_arr, est_h2=est_h2_arr, est_h2_lb=est_h2_lb_arr, est_h2_ub=est_h2_ub_arr)

	df$iteration = factor(df$iteration, levels=iter_names)
	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,300,1000,10000))

	red_colors=brewer.pal(n = 9, name = "Blues")


	pp<-ggplot(data=df, aes(x=iteration, fill=eqtl_sample_size, y=est_h2)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c(red_colors[3], red_colors[5], red_colors[7], red_colors[9]))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("non-mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) +
  		theme(legend.position="top") +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


  	return(pp)

}

make_se_barplot_for_med_h2_over_iterations <- function(visualize_trait_med_h2_dir, iter_names, method) {

	eqtl_sample_size_arr <- c()
	est_h2_arr <- c()
	est_h2_lb_arr <- c()
	est_h2_ub_arr <- c()
	iter_name_arr <- c()

	for (itera in 1:length(iter_names)) {
		iter_name <- iter_names[itera]

		input_file <- paste0(visualize_trait_med_h2_dir, "med_h2_results_summary_averaged", iter_name,".txt")
		tmp_df <-read.table(input_file, header=TRUE, sep="\t")

		tmp_df <- tmp_df[as.character(tmp_df$method)==method,]
		tmp_df <- tmp_df[as.character(tmp_df$heritability_type)=="alt_med_h2",]

		eqtl_sample_size_arr <- c(eqtl_sample_size_arr, tmp_df$eqtl_sample_size)
		est_h2_arr <- c(est_h2_arr, tmp_df$est_h2)
		est_h2_lb_arr <- c(est_h2_lb_arr, tmp_df$est_h2_lb)
		est_h2_ub_arr <- c(est_h2_ub_arr, tmp_df$est_h2_ub)

		iter_name_arr <- c(iter_name_arr, rep(iter_name, length(tmp_df$est_h2_ub)))

		sim_h2 = tmp_df$sim_h2[1]
	}

	df <- data.frame(iteration=iter_name_arr,eqtl_sample_size=eqtl_sample_size_arr, est_h2=est_h2_arr, est_h2_lb=est_h2_lb_arr, est_h2_ub=est_h2_ub_arr)

	df$iteration = factor(df$iteration, levels=iter_names)
	df$eqtl_sample_size <- factor(df$eqtl_sample_size, levels=c(100,300,1000,10000))

	red_colors=brewer.pal(n = 9, name = "Reds")


	pp<-ggplot(data=df, aes(x=iteration, fill=eqtl_sample_size, y=est_h2)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c(red_colors[3], red_colors[5], red_colors[7], red_colors[9]))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) +
  		theme(legend.position="top") +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


  	return(pp)

}




########################
# Command line args
########################
trait_med_h2_inference_dir = args[1]
visualize_trait_med_h2_dir = args[2]


# Load in data
summary_file <- paste0(visualize_trait_med_h2_dir, "med_h2_results_summary_averaged_wrong_eqtls.txt")
df <- read.table(summary_file, header=TRUE, sep="\t")

summary_file <- paste0(visualize_trait_med_h2_dir, "tglr_med_h2_results_summary_averaged.txt")
df_tglr <- read.table(summary_file, header=TRUE, sep="\t")


print(summary_file)

med_plot <- make_se_barplot_for_med_h2(df, "marginal_pca_sumstat_gibbs_0.0")
output_file <- paste0(visualize_trait_med_h2_dir, "marginal_pca_sumstat_gibbs_h2_est_se_barplot_wrong_etls.pdf")
ggsave(med_plot, file=output_file, width=7.2, height=4.5, units="in")

print(output_file)


##################
## Marginal PCA
##################
# Non-mediated plot
if (FALSE) {
nm_plot <- make_se_barplot_for_nm_h2(df, "marginal_pca_sumstat_gibbs_0.0")
# Mediated plot
med_plot <- make_se_barplot_for_med_h2(df, "marginal_pca_sumstat_gibbs_0.0")
# Total h2 plot
tot_plot <- make_se_barplot_for_total_h2(df, "marginal_pca_sumstat_gibbs_0.0")
# Combine together with cowplot
joint <- plot_grid(med_plot + poster_figure_theme(), nm_plot+ poster_figure_theme(), tot_plot+ poster_figure_theme(), ncol=1)
# Save
output_file <- paste0(visualize_trait_med_h2_dir, "marginal_pca_sumstat_gibbs_h2_est_se_barplot.pdf")
ggsave(joint, file=output_file, width=16.2, height=12.5, units="in")
}

if (FALSE) {
###############################################
# Plot showing joint vs MESC plot
# Mediated plot
med_plot <- make_se_barplot_for_med_h2_joint_vs_mesc(df, "marginal_pca_sumstat_gibbs_0.0",df_tglr) + theme(legend.position="bottom")
output_file <- paste0(visualize_trait_med_h2_dir, "marginal_pca_sumstat_gibbs_med_est_se_barplot_joint_vs_mesc.pdf")
ggsave(med_plot, file=output_file, width=16.2, height=9.7, units="in")

# Non-Mediated plot
med_plot <- make_se_barplot_for_nm_h2_joint_vs_mesc(df, "marginal_pca_sumstat_gibbs_0.0",df_tglr)
output_file <- paste0(visualize_trait_med_h2_dir, "marginal_pca_sumstat_gibbs_nm_est_se_barplot_joint_vs_mesc.pdf")
ggsave(med_plot, file=output_file, width=7.2, height=3.7, units="in")
}


###############################################
# Plot showing joint vs marginal effects
# Mediated plot
if (FALSE) {
med_plot <- make_se_barplot_for_med_h2_joint_vs_marginal(df, "marginal_pca_sumstat_gibbs_0.0")
output_file <- paste0(visualize_trait_med_h2_dir, "marginal_pca_sumstat_gibbs_med_est_se_barplot_joint_vs_marginal.pdf")
ggsave(med_plot, file=output_file, width=7.2, height=3.7, units="in")
print(output_file)
# Non-Mediated plot
med_plot <- make_se_barplot_for_nm_h2_joint_vs_marginal(df, "marginal_pca_sumstat_gibbs_0.0")
output_file <- paste0(visualize_trait_med_h2_dir, "marginal_pca_sumstat_gibbs_nm_est_se_barplot_joint_vs_marginal.pdf")
ggsave(med_plot, file=output_file, width=7.2, height=3.7, units="in")
}

###############################################
# Plot showing various values across iterations
if (FALSE) {
options("scipen"=10)
iter_values <- seq(0, 600000,25000)
iter_names <- c()
for (ii in 2:length(iter_values)) {
	iter_names <- c(iter_names, paste0(as.character(iter_values[(ii-1)]), ":", as.character(iter_values[ii])))
}

# Mediated plot
med_plot <- make_se_barplot_for_med_h2_over_iterations(visualize_trait_med_h2_dir, iter_names, "marginal_pca_sumstat_gibbs_0.0")
output_file <- paste0(visualize_trait_med_h2_dir, "marginal_pca_sumstat_gibbs_med_h2_est_se_barplot_across_iterations.pdf")
ggsave(med_plot, file=output_file, width=7.2, height=3.7, units="in")
print(output_file)

nm_plot <- make_se_barplot_for_nm_h2_over_iterations(visualize_trait_med_h2_dir, iter_names, "marginal_pca_sumstat_gibbs_0.0")
output_file <- paste0(visualize_trait_med_h2_dir, "marginal_pca_sumstat_gibbs_nm_h2_est_se_barplot_across_iterations.pdf")
ggsave(nm_plot, file=output_file, width=7.2, height=3.7, units="in")
print(output_file)

total_plot <- make_se_barplot_for_total_h2_over_iterations(visualize_trait_med_h2_dir, iter_names, "marginal_pca_sumstat_gibbs_0.0")
output_file <- paste0(visualize_trait_med_h2_dir, "marginal_pca_sumstat_gibbs_total_h2_est_se_barplot_across_iterations.pdf")
ggsave(total_plot, file=output_file, width=7.2, height=3.7, units="in")
print(output_file)
}
