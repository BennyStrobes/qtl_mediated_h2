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

make_multimethod_se_barplot_for_nm_h2 <- function(df, methods) {
	df2 <- df[as.character(df$method)%in%methods,]
	df2 <- df2[as.character(df2$heritability_type)=="nm_h2",]
	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,200,300,1000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("nm_h2"))

	sim_h2 = df2$sim_h2[1]

	df2$method = factor(df2$method, levels=methods)

	#blue_colors=brewer.pal(n = 9, name = "Blues")

	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("non-mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)


}


make_se_barplot_for_nm_h2 <- function(df, method) {
	df2 <- df[as.character(df$method)==method,]
	df2 <- df2[as.character(df2$heritability_type)=="nm_h2",]
	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,200,300,1000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("nm_h2"))

	sim_h2 = df2$sim_h2[1]

	blue_colors=brewer.pal(n = 9, name = "Blues")

	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2)) +
  		geom_bar(stat="identity", position=position_dodge(), fill=blue_colors[7]) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("non-mediated h2")) +
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


make_multimethod_se_barplot_for_causal_tissue_med_h2_v1 <- function(df, methods) {
	df2 <- df[as.character(df$method)%in%methods,]
	df2 <- df2[as.character(df2$heritability_type)=="causal_tissue_med_h2",]
	df2$eqtl_sample_size = as.character(df2$eqtl_sample_size)
	valid_ss=c("100","200","300","1000")
	df2 <- df2[as.character(df2$eqtl_sample_size)%in%valid_ss,]

	print(df2)

	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c("100","200","300","1000"))
	df2$heritability_type = factor(df2$heritability_type, levels=c("causal_tissue_med_h2"))

	sim_h2 = df2$sim_h2[1]

	df2$method = factor(df2$method,levels=methods)


	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2, fill=snp_representation)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("Causal tissue mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)
}



make_multimethod_se_barplot_for_causal_tissue_med_h2 <- function(df, methods) {
	df2 <- df[as.character(df$method)%in%methods,]
	df2 <- df2[as.character(df2$heritability_type)=="causal_tissue_med_h2",]
	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c("100", "100-1000","200","300","1000"))
	df2$heritability_type = factor(df2$heritability_type, levels=c("causal_tissue_med_h2"))

	sim_h2 = df2$sim_h2[1]

	df2$method = factor(df2$method,levels=methods)


	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("Causal tissue mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)
}

make_multimethod_se_barplot_for_noncausal_tissue_med_h2_v1 <- function(df, methods) {
	df2 <- df[as.character(df$method)%in%methods,]
	df2 <- df2[as.character(df2$heritability_type)=="non_causal_tissue_med_h2",]
	
	df2$eqtl_sample_size = as.character(df2$eqtl_sample_size)
	df2 <- df2[as.character(df2$eqtl_sample_size)%in%c("100","200","300","1000"),]



	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c("100","200","300","1000"))
	df2$heritability_type = factor(df2$heritability_type, levels=c("non_causal_tissue_med_h2"))

	sim_h2 = df2$sim_h2[1]

	df2$method = factor(df2$method,levels=methods)


	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2, fill=snp_representation)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("Aggregate non-Causal tissue mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)
}


make_multimethod_se_barplot_for_noncausal_tissue_med_h2 <- function(df, methods) {
	df2 <- df[as.character(df$method)%in%methods,]
	df2 <- df2[as.character(df2$heritability_type)=="non_causal_tissue_med_h2",]
	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c("100", "100-1000","200","300","1000"))
	df2$heritability_type = factor(df2$heritability_type, levels=c("non_causal_tissue_med_h2"))

	sim_h2 = df2$sim_h2[1]

	df2$method = factor(df2$method,levels=methods)


	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("Aggregate non-Causal tissue mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)
}


make_single_method_mediated_h2_se_barplot_for_permuted_and_unpermuted <- function(df, method) {
	df2 <- df[as.character(df$method)==method,]
	df2 <- df2[as.character(df2$heritability_type)=="med_h2",]
	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,300,1000,10000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("med_h2"))

	sim_h2 = df2$sim_h2[1]

	df2$permuted <- factor(df2$permuted, levels=c(TRUE, FALSE))

	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2, fill=permuted)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c("grey", "skyblue3"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="Permuted data", title=paste0(method, " mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)


}

make_multimethod_se_barplot_for_med_h2 <- function(df, methods) {
	df2 <- df[as.character(df$method)%in%methods,]
	df2 <- df2[as.character(df2$heritability_type)=="med_h2",]
	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,200,300,1000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("med_h2"))

	sim_h2 = df2$sim_h2[1]

	df2$method = factor(df2$method,levels=methods)


	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)
}

make_se_barplot_for_med_h2 <- function(df, method) {
	df2 <- df[as.character(df$method)==method,]
	df2 <- df2[as.character(df2$heritability_type)=="med_h2",]
	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,300,1000,10000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("med_h2"))

	sim_h2 = df2$sim_h2[1]

	red_colors=brewer.pal(n = 9, name = "Reds")

	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2)) +
  		geom_bar(stat="identity", position=position_dodge(), fill=red_colors[7]) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("mediated h2")) +
  		geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)
}

make_se_barplot_for_per_tissue_category_med_h2 <- function(df, method) {
	df2 <- df[as.character(df$method)==method,]
	df2 <- df2[(as.character(df2$heritability_type)!="med_h2_tissue2") &(as.character(df2$heritability_type)!="nm_h2")&(as.character(df2$heritability_type)!="total_h2")&(as.character(df2$heritability_type)!="med_h2_tissue0")&(as.character(df2$heritability_type)!="med_h2_tissue1")&(as.character(df2$heritability_type)!="med_h2_tissue3")&(as.character(df2$heritability_type)!="med_h2_tissue4") ,]


	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,300,1000,10000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("med_h2", "causal_tissue_med_h2", "non_causal_tissue_med_h2"))

	sim_h2 = df2$sim_h2[1]

	print(df2)


	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, fill=heritability_type, y=est_h2)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c(red_colors[6], "grey", "grey", "grey", "grey"))+
  		figure_theme() +
  		labs(x="QTL sample size", y="heritability", fill="", title=paste0("mediated heritability")) +
  		geom_hline(yintercept=sim_h2, linetype=2)  +
  		theme(legend.position="bottom")+guides(fill=guide_legend(nrow=2,byrow=TRUE))

  	return(pp)
}


make_se_barplot_for_per_tissue_med_h2 <- function(df, method) {
	df2 <- df[as.character(df$method)==method,]
	df2 <- df2[(as.character(df2$heritability_type)!="med_h2") &(as.character(df2$heritability_type)!="nm_h2")&(as.character(df2$heritability_type)!="total_h2")&(as.character(df2$heritability_type)!="causal_tissue_med_h2")&(as.character(df2$heritability_type)!="non_causal_tissue_med_h2") ,]


	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,300,1000,10000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("med_h2_tissue0", "med_h2_tissue1", "med_h2_tissue2", "med_h2_tissue3", "med_h2_tissue4"))

	sim_h2 = df2$sim_h2[1]

	red_colors=brewer.pal(n = 9, name = "Reds")

	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, fill=heritability_type, y=est_h2)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		scale_fill_manual(values=c(red_colors[6], "grey", "grey", "grey", "grey"))+
  		figure_theme() +
  		labs(x="QTL sample size", y="heritability", fill="", title=paste0("mediated heritability")) +
  		geom_hline(yintercept=sim_h2, linetype=2)  +
  		theme(legend.position="bottom")+guides(fill=guide_legend(nrow=2,byrow=TRUE))

  	return(pp)
}


make_multi_method_mediated_h2_se_barplot_for_permuted_and_unpermuted_ratio <- function(df, methods) {
	df2 <- df[as.character(df$method)%in%methods,]
	df2 <- df2[as.character(df2$heritability_type)=="med_h2",]

	med_h2_ratio_arr <- c()
	method_arr <- c()
	eqtl_ss_arr <- c()

	eqtl_sample_sizes <- c(100, 300, 1000, 10000)

	for (ss_iter in 1:length(eqtl_sample_sizes)) {
		eqtl_ss <- eqtl_sample_sizes[ss_iter]
		for (method_iter in 1:length(methods)) {
			method <- as.character(methods[method_iter])

			df3 <- df2[as.character(df2$method) == method,]
			df3 <- df3[df3$eqtl_sample_size == eqtl_ss,]

			med_h2_ratio <- df3$est_h2[2]/df3$est_h2[1]

			med_h2_ratio_arr <- c(med_h2_ratio_arr, med_h2_ratio)
			method_arr <- c(method_arr, method)
			eqtl_ss_arr <- c(eqtl_ss_arr, eqtl_ss)
		}
	}

	df_final <- data.frame(method=factor(method_arr, levels=methods), eqtl_ss=factor(eqtl_ss_arr, levels=eqtl_sample_sizes), med_h2_ratio=med_h2_ratio_arr)

	pp<-ggplot(data=df_final, aes(x=eqtl_ss, fill=method, y=med_h2_ratio)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		figure_theme() +
  		labs(x="QTL sample size", y="permuted med h2 / un-permuted med h2", fill="",) +
  		theme(legend.position="bottom")

  	return(pp)

}


make_multimethod_se_barplot_for_total_h2 <- function(df, methods) {
	df2 <- df[as.character(df$method)%in%methods,]
	df2 <- df2[as.character(df2$heritability_type)=="total_h2",]
	
	df2$eqtl_sample_size <- factor(df2$eqtl_sample_size, levels=c(100,200,300,1000))
	df2$heritability_type = factor(df2$heritability_type, levels=c("total_h2"))

	sim_h2 = df2$sim_h2[1]

	df2$method = factor(df2$method,levels=methods)

	pp<-ggplot(data=df2, aes(x=eqtl_sample_size, y=est_h2, fill=method)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("Total h2")) +
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
  		geom_bar(stat="identity", position=position_dodge(), fill='grey') +
  		geom_errorbar(aes(ymin=est_h2_lb, ymax=est_h2_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="eQTL SS", y="h2", fill="", title=paste0("Total h2")) +
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

make_histogram_showing_distribution_of_tissue_0_estimates <- function(per_sim_df, eqtl_sample_size) {
	df <- per_sim_df[as.character(per_sim_df$eqtl_sample_size) == eqtl_sample_size,]

	pp <-ggplot(data = df, aes(x = est_tissue0_joint)) +
  		geom_histogram() +
  		figure_theme() + 
  		labs(x="Tissue 0 mediated h2 estimate", title=paste0("eQTL Sample Size: ", eqtl_sample_size))
  	return(pp)
}

make_coverage_standard_error_barplot <- function(df_cal, genetic_element_class) {




	df <- df_cal[as.character(df_cal$genetic_element) == genetic_element_class,]

	df$eQTL_ss <- factor(df$eQTL_ss, levels=c("100", "100-1000", "200", "300", "1000"))
	df$expected_coverage <- factor(df$expected_coverage, levels=c(0.50,0.70, 0.90, 0.95))



df_lines <- data.frame(
  expected_coverage = factor(c(0.50,0.70,0.90,0.95),
                             levels = levels(df$expected_coverage)),
  hline_y          = c(0.50, 0.70, 0.90, 0.95)
)

	pp<-ggplot(data=df, aes(x=expected_coverage, y=observed_coverage, fill=eQTL_ss)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=observed_coverage_lb, ymax=observed_coverage_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="Expected coverage", y="Coverage", fill="eQTL SS", title=genetic_element_class)


pp=pp +
  geom_errorbar(
    data       = df_lines,
    aes(x       = expected_coverage,
        ymin    = hline_y,
        ymax    = hline_y),
    width      = 0.9,        # matches the full width of a bar group
    inherit.aes = FALSE,     # don’t carry over your original aes()
    color      = "grey40",    # pick whatever color you like
    size       = .5
  )

  	return(pp)

}


make_power_standard_error_barplot <- function(df_cal, genetic_element_class) {




	df <- df_cal[as.character(df_cal$genetic_element) == genetic_element_class,]

	df <- df[as.character(df$method_name) == "calibrated_mesc",]

	df$eQTL_ss <- factor(df$eQTL_ss, levels=c("100", "100-1000", "200", "300", "1000"))


	df$pvalue_threshold <- factor(df$pvalue_threshold, levels=c(0.20,0.10, 0.05, 0.01))

df_lines <- data.frame(
  pvalue_threshold = factor(c(0.20,0.10, 0.05, 0.01),
                             levels = levels(df$pvalue_threshold)),
  hline_y          = c(0.20,0.10, 0.05, 0.01)
)



	pp<-ggplot(data=df, aes(x=pvalue_threshold, y=power, fill=eQTL_ss)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		labs(x="p-value threshold", y="Power", fill="eQTL SS", title=genetic_element_class)

pp=pp +
  geom_errorbar(
    data       = df_lines,
    aes(x       = pvalue_threshold,
        ymin    = hline_y,
        ymax    = hline_y),
    width      = 0.9,        # matches the full width of a bar group
    inherit.aes = FALSE,     # don’t carry over your original aes()
    color      = "grey40",    # pick whatever color you like
    size       = .5
  )

  	return(pp)

}

make_fstat_violinplot_w_stdExpr_data <- function(df, remove_1000=FALSE) {
	#df <- fstat_df[as.character(fstat_df$method_name) == "category0",]
	if (remove_1000) {
		df <- df[as.character(df$eQTL_ss) != "1000",]
		df$eQTL_ss <- factor(df$eQTL_ss, levels=c("100", "100-1000", "200", "300"))
	} else {
		df$eQTL_ss <- factor(df$eQTL_ss, levels=c("100", "100-1000", "200", "300", "1000"))
	}


	pp<-ggplot(df, aes(x=eQTL_ss, y=fstat, color=method_name)) +
  		geom_violin(trim=FALSE) +
  		figure_theme() + 
  		labs(x="eQTL SS", y="F statistic", color="cis-snp h2 bin") 

  	return(pp)

}


make_fstat_violinplot <- function(fstat_df, remove_1000=FALSE) {
	df <- fstat_df[as.character(fstat_df$method_name) == "category0",]
	if (remove_1000) {
		df <- df[as.character(df$eQTL_ss) != "1000",]
		df$eQTL_ss <- factor(df$eQTL_ss, levels=c("100", "100-1000", "200", "300"))
	} else {
		df$eQTL_ss <- factor(df$eQTL_ss, levels=c("100", "100-1000", "200", "300", "1000"))
	}

	pp<-ggplot(df, aes(x=eQTL_ss, y=fstat, color=eQTL_ss)) +
  		geom_violin(trim=FALSE) +
  		figure_theme() + 
  		labs(x="eQTL SS", y="F statistic", color="") +
  		theme(legend.position="none")

  	return(pp)

}


make_fstat_method_comparison_violinplot <- function(fstat_df,fstat_random_bins_df, remove_1000=FALSE) {
	df <- fstat_df[as.character(fstat_df$method_name) == "category0",]
	df2 <- fstat_random_bins_df[as.character(fstat_random_bins_df$method_name) == "category0",]
	if (remove_1000) {
		df <- df[as.character(df$eQTL_ss) != "1000",]
		df$eQTL_ss <- factor(df$eQTL_ss, levels=c("100", "100-1000", "200", "300"))
		df2 <- df2[as.character(df2$eQTL_ss) != "1000",]
		df2$eQTL_ss <- factor(df2$eQTL_ss, levels=c("100", "100-1000", "200", "300"))


	} else {
		df$eQTL_ss <- factor(df$eQTL_ss, levels=c("100", "100-1000", "200", "300", "1000"))
		d2f$eQTL_ss <- factor(df$eQTL_ss, levels=c("100", "100-1000", "200", "300", "1000"))
	}


	df$source <- "single_gene_bin"
	df2$source <- "four_gene_bins"

	# stack them
	new_df <- rbind(df, df2)

	# make the new column a factor if you like
	new_df$source <- factor(new_df$source, levels = c("single_gene_bin","four_gene_bins"))


	pp<-ggplot(new_df, aes(x=eQTL_ss, y=fstat, color=source)) +
  		geom_violin(trim=FALSE) +
  		figure_theme() + 
  		labs(x="eQTL SS", y="F statistic", color="")

  	return(pp)

}




########################
# Command line args
########################
trait_med_h2_inference_dir = args[1]
organized_trait_med_h2_results_dir = args[2]
visualize_trait_med_h2_dir = args[3]

if (FALSE) {
###################################################
# Joint plot showing mediated heritability in causal tissue and mediated heritability in non-causal tissues
# Make plot for causal tissue
eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="linear"
sq_sumstat_threshold="100.0"
weighting = "unweighted"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture,"_squared_marginal_sumstats_", sq_sumstat_threshold,"_", weighting)


joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_calibrated_ldsc_summary_averaged.txt")
joint_ldsc_df <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")

causal_med_plot <- make_multimethod_se_barplot_for_causal_tissue_med_h2(joint_ldsc_df, c("uncalibrated_mesc","calibrated_mesc"))
# Make plot for non-causal tissue
non_causal_med_plot <- make_multimethod_se_barplot_for_noncausal_tissue_med_h2(joint_ldsc_df, c("uncalibrated_mesc","calibrated_mesc"))
# Get legend
legender <- get_legend(non_causal_med_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(causal_med_plot+ theme(legend.position="none"), non_causal_med_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_calibrated_mesc_tissue_simulation_causal_tissue_med_vs_non_causal_tissue_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=5.5, units="in")


###################################################
# Joint plot showing mediated heritability in causal tissue and mediated heritability in non-causal tissues
# Make plot for causal tissue
eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="linear"
sq_sumstat_threshold="100.0"
weighting = "weighted"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture,"_squared_marginal_sumstats_", sq_sumstat_threshold,"_", weighting)


joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_calibrated_ldsc_summary_averaged.txt")
joint_ldsc_df <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")

causal_med_plot <- make_multimethod_se_barplot_for_causal_tissue_med_h2(joint_ldsc_df, c("uncalibrated_mesc","calibrated_mesc"))
# Make plot for non-causal tissue
non_causal_med_plot <- make_multimethod_se_barplot_for_noncausal_tissue_med_h2(joint_ldsc_df, c("uncalibrated_mesc","calibrated_mesc"))
# Get legend
legender <- get_legend(non_causal_med_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(causal_med_plot+ theme(legend.position="none"), non_causal_med_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_calibrated_mesc_tissue_simulation_causal_tissue_med_vs_non_causal_tissue_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=5.5, units="in")


###################################################
# Joint plot showing mediated heritability in causal tissue and mediated heritability in non-causal tissues
# Make plot for causal tissue
eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="random_bins"
sq_sumstat_threshold="100.0"
weighting = "unweighted"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture,"_squared_marginal_sumstats_", sq_sumstat_threshold,"_", weighting)


joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_calibrated_ldsc_summary_averaged.txt")
joint_ldsc_df <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")

causal_med_plot <- make_multimethod_se_barplot_for_causal_tissue_med_h2(joint_ldsc_df, c("uncalibrated_mesc","calibrated_mesc"))
# Make plot for non-causal tissue
non_causal_med_plot <- make_multimethod_se_barplot_for_noncausal_tissue_med_h2(joint_ldsc_df, c("uncalibrated_mesc","calibrated_mesc"))
# Get legend
legender <- get_legend(non_causal_med_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(causal_med_plot+ theme(legend.position="none"), non_causal_med_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_calibrated_mesc_tissue_simulation_causal_tissue_med_vs_non_causal_tissue_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=5.5, units="in")



###################################################
# Joint plot showing mediated heritability in causal tissue and mediated heritability in non-causal tissues
# Make plot for causal tissue
eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="linear"
sq_sumstat_threshold="100.0"
weighting = "unweighted"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture,"_squared_marginal_sumstats_", sq_sumstat_threshold,"_", weighting)


joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_calibrated_ldsc_summary_variance_weighted_averaged.txt")
joint_ldsc_df <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")

causal_med_plot <- make_multimethod_se_barplot_for_causal_tissue_med_h2(joint_ldsc_df, c("uncalibrated_mesc","calibrated_mesc"))
# Make plot for non-causal tissue
non_causal_med_plot <- make_multimethod_se_barplot_for_noncausal_tissue_med_h2(joint_ldsc_df, c("uncalibrated_mesc","calibrated_mesc"))
# Get legend
legender <- get_legend(non_causal_med_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(causal_med_plot+ theme(legend.position="none"), non_causal_med_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_calibrated_mesc_tissue_simulation_causal_tissue_med_vs_non_causal_tissue_variance_weighted_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=5.5, units="in")

}


###################################################
# Make bookstrapped standard error coverage plots
if (FALSE) {
eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="linear"
sq_sumstat_threshold="100.0"
weighting = "unweighted"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture,"_squared_marginal_sumstats_", sq_sumstat_threshold,"_", weighting)


calibration_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_calibration_summary.txt")
df_cal <- read.table(calibration_summary_file, header=TRUE, sep="\t")

genetic_element_class <- "total_nm_h2"
pp1 <- make_coverage_standard_error_barplot(df_cal, genetic_element_class)

genetic_element_class <- "total_med_h2"
pp2 <- make_coverage_standard_error_barplot(df_cal, genetic_element_class)

genetic_element_class <- "causal_tissue_med_h2"
pp3 <- make_coverage_standard_error_barplot(df_cal, genetic_element_class)

joint_pp <- plot_grid(pp1, pp2, pp3, ncol=1)

# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_coverage_calibrated_mesc_se_barplot.pdf")
ggsave(joint_pp, file=output_file, width=7.2, height=5.5, units="in")
}




###################################################
# Make bookstrapped power plots
if (FALSE) {
eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="linear"
sq_sumstat_threshold="100.0"
weighting = "unweighted"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture,"_squared_marginal_sumstats_", sq_sumstat_threshold,"_", weighting)


power_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_power_summary.txt")
df_power <- read.table(power_summary_file, header=TRUE, sep="\t")

genetic_element_class <- "total_nm_h2"
pp1 <- make_power_standard_error_barplot(df_power, genetic_element_class)

genetic_element_class <- "total_med_h2"
pp2 <- make_power_standard_error_barplot(df_power, genetic_element_class)

genetic_element_class <- "causal_tissue_med_h2"
pp3 <- make_power_standard_error_barplot(df_power, genetic_element_class)

joint_pp <- plot_grid(pp1, pp2, pp3, ncol=1)

# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_power_calibrated_mesc_se_barplot.pdf")
ggsave(joint_pp, file=output_file, width=7.2, height=5.5, units="in")
}





###################################################
# Make F-statistic plots
eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="linear"
sq_sumstat_threshold="100.0"
weighting = "unweighted"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture,"_squared_marginal_sumstats_", sq_sumstat_threshold,"_", weighting)
if (FALSE) {
fstat_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_fstat_summary.txt")
fstat_df <- read.table(fstat_summary_file, header=TRUE, sep="\t")

pp <- make_fstat_violinplot(fstat_df)
# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_fstat_violinplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")


pp <- make_fstat_violinplot(fstat_df, remove_1000=TRUE)
# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_fstat_remove_1000_violinplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.5, units="in")
}

###################################################
# Make F-statistic plots
eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="stdExpr"
inference_gt_architecture="stdExpr"
sq_sumstat_threshold="100.0"
weighting = "unweighted"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture,"_squared_marginal_sumstats_", sq_sumstat_threshold,"_", weighting)
fstat_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_fstat_summary.txt")
fstat_df <- read.table(fstat_summary_file, header=TRUE, sep="\t")
pp <- make_fstat_violinplot_w_stdExpr_data(fstat_df, remove_1000=TRUE)
# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_fstat_stdExpr_remove_1000_violinplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=3.79, units="in")


###################################################
# Make F-statistic plots comparing linear model + random bins
if (FALSE) {
eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="linear"
sq_sumstat_threshold="100.0"
weighting = "unweighted"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture,"_squared_marginal_sumstats_", sq_sumstat_threshold,"_", weighting)
fstat_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_fstat_summary.txt")
fstat_df <- read.table(fstat_summary_file, header=TRUE, sep="\t")
simulated_gt_architecture="linear"
inference_gt_architecture="random_bins"
sq_sumstat_threshold="100.0"
weighting = "unweighted"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture,"_squared_marginal_sumstats_", sq_sumstat_threshold,"_", weighting)
fstat_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_fstat_summary.txt")
fstat_random_bins_df <- read.table(fstat_summary_file, header=TRUE, sep="\t")

pp <- make_fstat_method_comparison_violinplot(fstat_df,fstat_random_bins_df, remove_1000=TRUE)

output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_fstat_method_comparison_remove_1000_violinplot.pdf")
ggsave(pp, file=output_file, width=7.2, height=4.0, units="in")

}



















if (FALSE) {
###################################################
# Joint plot showing mediated heritability in causal tissue and mediated heritability in non-causal tissues
# Make plot for causal tissue

eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="stdExpr"
inference_gt_architecture="stdExpr"
sq_sumstat_threshold="100.0"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture,"_squared_marginal_sumstats_", sq_sumstat_threshold)



joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_calibrated_ldsc_summary_averaged.txt")
joint_ldsc_df <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")

causal_med_plot <- make_multimethod_se_barplot_for_causal_tissue_med_h2(joint_ldsc_df, c("two_step_ldsc","calibrated_two_step_ldsc"))
# Make plot for non-causal tissue
non_causal_med_plot <- make_multimethod_se_barplot_for_noncausal_tissue_med_h2(joint_ldsc_df, c("two_step_ldsc","calibrated_two_step_ldsc"))
# Get legend
legender <- get_legend(non_causal_med_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(causal_med_plot+ theme(legend.position="none"), non_causal_med_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_calibrated_mesc_tissue_simulation_causal_tissue_med_vs_non_causal_tissue_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=5.5, units="in")



###################################################
# Joint plot showing mediated heritability in causal tissue and mediated heritability in non-causal tissues
# Make plot for causal tissue

eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="stdExpr"
inference_gt_architecture="linear"
sq_sumstat_threshold="100.0"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture,"_squared_marginal_sumstats_", sq_sumstat_threshold)



joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_calibrated_ldsc_summary_averaged.txt")
joint_ldsc_df <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")

causal_med_plot <- make_multimethod_se_barplot_for_causal_tissue_med_h2(joint_ldsc_df, c("two_step_ldsc","calibrated_two_step_ldsc"))
# Make plot for non-causal tissue
non_causal_med_plot <- make_multimethod_se_barplot_for_noncausal_tissue_med_h2(joint_ldsc_df, c("two_step_ldsc","calibrated_two_step_ldsc"))
# Get legend
legender <- get_legend(non_causal_med_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(causal_med_plot+ theme(legend.position="none"), non_causal_med_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_calibrated_mesc_tissue_simulation_causal_tissue_med_vs_non_causal_tissue_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=5.5, units="in")



###################################################
# Joint plot showing mediated heritability in causal tissue and mediated heritability in non-causal tissues
# Make plot for causal tissue

eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="random_bins"
sq_sumstat_threshold="100"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture,"_squared_marginal_sumstats_", sq_sumstat_threshold)



joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_calibrated_ldsc_summary_averaged.txt")
joint_ldsc_df <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")

causal_med_plot <- make_multimethod_se_barplot_for_causal_tissue_med_h2(joint_ldsc_df, c("two_step_ldsc","calibrated_two_step_ldsc"))
# Make plot for non-causal tissue
non_causal_med_plot <- make_multimethod_se_barplot_for_noncausal_tissue_med_h2(joint_ldsc_df, c("two_step_ldsc","calibrated_two_step_ldsc"))
# Get legend
legender <- get_legend(non_causal_med_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(causal_med_plot+ theme(legend.position="none"), non_causal_med_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_calibrated_mesc_tissue_simulation_causal_tissue_med_vs_non_causal_tissue_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=5.5, units="in")



###################################################
# Joint plot showing mediated heritability in causal tissue and mediated heritability in non-causal tissues
# Make plot for causal tissue

eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="stdExpr"
sq_sumstat_threshold="100"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture,"_squared_marginal_sumstats_", sq_sumstat_threshold)



joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_calibrated_ldsc_summary_averaged.txt")
joint_ldsc_df <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")

causal_med_plot <- make_multimethod_se_barplot_for_causal_tissue_med_h2(joint_ldsc_df, c("two_step_ldsc","calibrated_two_step_ldsc"))
# Make plot for non-causal tissue
non_causal_med_plot <- make_multimethod_se_barplot_for_noncausal_tissue_med_h2(joint_ldsc_df, c("two_step_ldsc","calibrated_two_step_ldsc"))
# Get legend
legender <- get_legend(non_causal_med_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(causal_med_plot+ theme(legend.position="none"), non_causal_med_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_calibrated_mesc_tissue_simulation_causal_tissue_med_vs_non_causal_tissue_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=5.5, units="in")

}




if (FALSE) {
eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="linear"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture)

joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_joint_ldsc_binned_summary_averaged.txt")
joint_ldsc_df <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")
print(joint_ldsc_df)

###################################################
# Joint plot showing mediated heritability in causal tissue and mediated heritability in non-causal tissues
# Make plot for causal tissue
causal_med_plot <- make_multimethod_se_barplot_for_causal_tissue_med_h2(joint_ldsc_df, c("two_step_ldsc","joint_ldsc"))
# Make plot for non-causal tissue
non_causal_med_plot <- make_multimethod_se_barplot_for_noncausal_tissue_med_h2(joint_ldsc_df, c("two_step_ldsc","joint_ldsc"))
# Get legend
legender <- get_legend(non_causal_med_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(causal_med_plot+ theme(legend.position="none"), non_causal_med_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_tissue_simulation_causal_tissue_med_vs_non_causal_tissue_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.5, units="in")
}


if (FALSE) {
eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="linear"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture)

joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_joint_ldsc_binned_summary_averaged.txt")
joint_ldsc_df_bins20 <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")
joint_ldsc_df_bins20['snp_representation'] = rep(eqtl_snp_representation, dim(joint_ldsc_df_bins20)[1])
print(joint_ldsc_df_bins20)

eqtl_snp_representation="bins_10"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="linear"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture)

joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_joint_ldsc_binned_summary_averaged.txt")
joint_ldsc_df_bins10 <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")
joint_ldsc_df_bins10['snp_representation'] = rep(eqtl_snp_representation, dim(joint_ldsc_df_bins10)[1])

eqtl_snp_representation="pca_90"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="linear"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture)

joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_joint_ldsc_binned_summary_averaged.txt")
joint_ldsc_df_pca90 <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")
joint_ldsc_df_pca90['snp_representation'] = rep(eqtl_snp_representation, dim(joint_ldsc_df_pca90)[1])


joint_ldsc_df <- rbind(joint_ldsc_df_bins10, joint_ldsc_df_bins20, joint_ldsc_df_pca90)

###################################################
# Joint plot showing mediated heritability in causal tissue and mediated heritability in non-causal tissues
# Make plot for causal tissue
causal_med_plot <- make_multimethod_se_barplot_for_causal_tissue_med_h2_v1(joint_ldsc_df, c("joint_ldsc"))
# Make plot for non-causal tissue
non_causal_med_plot <- make_multimethod_se_barplot_for_noncausal_tissue_med_h2_v1(joint_ldsc_df, c("joint_ldsc"))
# Get legend
legender <- get_legend(non_causal_med_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(causal_med_plot+ theme(legend.position="none"), non_causal_med_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, "cmp_eqtl_snp_reps_tissue_simulation_causal_tissue_med_vs_non_causal_tissue_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.5, units="in")
}



if (FALSE) {
eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="stdExpr"
inference_gt_architecture="stdExpr"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture)

joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_joint_ldsc_binned_summary_averaged.txt")
joint_ldsc_df <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")
print(joint_ldsc_df)

###################################################
# Joint plot showing mediated heritability in causal tissue and mediated heritability in non-causal tissues
# Make plot for causal tissue
causal_med_plot <- make_multimethod_se_barplot_for_causal_tissue_med_h2(joint_ldsc_df, c("two_step_ldsc","joint_ldsc"))
# Make plot for non-causal tissue
non_causal_med_plot <- make_multimethod_se_barplot_for_noncausal_tissue_med_h2(joint_ldsc_df, c("two_step_ldsc","joint_ldsc"))
# Get legend
legender <- get_legend(non_causal_med_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(causal_med_plot+ theme(legend.position="none"), non_causal_med_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_tissue_simulation_causal_tissue_med_vs_non_causal_tissue_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.5, units="in")
}


if (FALSE) {
eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="stdExpr"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture)

joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_joint_ldsc_binned_summary_averaged.txt")
joint_ldsc_df <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")
print(joint_ldsc_df)

###################################################
# Joint plot showing mediated heritability in causal tissue and mediated heritability in non-causal tissues
# Make plot for causal tissue
causal_med_plot <- make_multimethod_se_barplot_for_causal_tissue_med_h2(joint_ldsc_df, c("two_step_ldsc","joint_ldsc"))
# Make plot for non-causal tissue
non_causal_med_plot <- make_multimethod_se_barplot_for_noncausal_tissue_med_h2(joint_ldsc_df, c("two_step_ldsc","joint_ldsc"))
# Get legend
legender <- get_legend(non_causal_med_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(causal_med_plot+ theme(legend.position="none"), non_causal_med_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, run_string, "_tissue_simulation_causal_tissue_med_vs_non_causal_tissue_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.5, units="in")

}






if (FALSE) {
eqtl_snp_representation="bins_20"
non_med_anno="genotype_intercept"
simulated_gt_architecture="linear"
inference_gt_architecture="linear"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture)

joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_joint_ldsc_binned_summary_averaged.txt")
joint_ldsc_df_bins20 <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")
joint_ldsc_df_bins20['snp_representation'] = rep(non_med_anno, dim(joint_ldsc_df_bins20)[1])

eqtl_snp_representation="bins_20"
non_med_anno="full_anno"
simulated_gt_architecture="linear"
inference_gt_architecture="linear"
run_string = paste0(eqtl_snp_representation, "_", non_med_anno, "_", simulated_gt_architecture, "_",inference_gt_architecture)

joint_ldsc_summary_file <- paste0(organized_trait_med_h2_results_dir, "med_h2_5_causal_tissue_", run_string,"_sim_results_joint_ldsc_binned_summary_averaged.txt")
joint_ldsc_df_bins10 <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")
joint_ldsc_df_bins10['snp_representation'] = rep(non_med_anno, dim(joint_ldsc_df_bins10)[1])


joint_ldsc_df <- rbind(joint_ldsc_df_bins10, joint_ldsc_df_bins20)

joint_ldsc_df$snp_representation = factor(joint_ldsc_df$snp_representation, levels=c("genotype_intercept", "full_anno"))

print(joint_ldsc_df)


###################################################
# Joint plot showing mediated heritability in causal tissue and mediated heritability in non-causal tissues
# Make plot for causal tissue
causal_med_plot <- make_multimethod_se_barplot_for_causal_tissue_med_h2_v1(joint_ldsc_df, c("joint_ldsc"))
# Make plot for non-causal tissue
non_causal_med_plot <- make_multimethod_se_barplot_for_noncausal_tissue_med_h2_v1(joint_ldsc_df, c("joint_ldsc"))
# Get legend
legender <- get_legend(non_causal_med_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(causal_med_plot+ theme(legend.position="none"), non_causal_med_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, "cmp_non_med_anno_reps_tissue_simulation_causal_tissue_med_vs_non_causal_tissue_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.5, units="in")

}



























if (FALSE) {
###################################################
# Joint plot showing mediated heritability, non-mediated heritability, and total heritablity for multiple methods
# Non-mediated plot
nm_plot <- make_multimethod_se_barplot_for_nm_h2(joint_ldsc_df, c("two_step_ldsc","joint_ldsc"))
# Mediated plot
med_plot <- make_multimethod_se_barplot_for_med_h2(joint_ldsc_df, c("two_step_ldsc","joint_ldsc"))
# Total h2 plot
tot_plot <- make_multimethod_se_barplot_for_total_h2(joint_ldsc_df, c("two_step_ldsc","joint_ldsc"))

# Get legend
legender <- get_legend(nm_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(med_plot+ theme(legend.position="none"), nm_plot+ theme(legend.position="none"), tot_plot+ theme(legend.position="none"),legender, ncol=1, rel_heights=c(1,1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, "multimethod_5_tissue_simulation_med_nm_total_h2_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.5, units="in")
print(output_file)


}

if (FALSE) {
# Joint plot showing mediated heritability, non-mediated heritability, and total heritablity for multiple methods
# Non-mediated plot
nm_plot <- make_multimethod_se_barplot_for_nm_h2(joint_ldsc_df, c("joint_ldsc", "joint_ldsc_cis_var"))
# Mediated plot
med_plot <- make_multimethod_se_barplot_for_med_h2(joint_ldsc_df, c("joint_ldsc", "joint_ldsc_cis_var"))
# Total h2 plot
tot_plot <- make_multimethod_se_barplot_for_total_h2(joint_ldsc_df, c("joint_ldsc", "joint_ldsc_cis_var"))

# Get legend
legender <- get_legend(nm_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(med_plot+ theme(legend.position="none"), nm_plot+ theme(legend.position="none"), tot_plot+ theme(legend.position="none"),legender, ncol=1, rel_heights=c(1,1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, "multimethod_cis_var_5_tissue_simulation_med_nm_total_h2_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.5, units="in")
print(output_file)
}


###################################################
if (FALSE) {
# Joint plot showing mediated heritability in causal tissue and mediated heritability in non-causal tissues
# Make plot for causal tissue
causal_med_plot <- make_multimethod_se_barplot_for_causal_tissue_med_h2(joint_ldsc_df, c("joint_ldsc", "joint_ldsc_cis_var"))
# Make plot for non-causal tissue
non_causal_med_plot <- make_multimethod_se_barplot_for_noncausal_tissue_med_h2(joint_ldsc_df, c("joint_ldsc", "joint_ldsc_cis_var"))
# Get legend
legender <- get_legend(non_causal_med_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(causal_med_plot+ theme(legend.position="none"), non_causal_med_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, "multimethod_cisvar_5_tissue_simulation_causal_tissue_med_vs_non_causal_tissue_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.5, units="in")
print(output_file)

}













##############
# OLD
##############

if (FALSE) {
###################################################
# Load in data
# Joint ldsc
joint_ldsc_summary_file <- paste0(visualize_trait_med_h2_dir, "med_h2_5_causal_tissue_sim_results_joint_ldsc_binned_summary_averaged.txt")
joint_ldsc_df <- read.table(joint_ldsc_summary_file, header=TRUE, sep="\t")
# MESC
mesc_summary_file <- paste0(visualize_trait_med_h2_dir, "med_h2_5_causal_tissue_sim_results_mesc_summary_averaged.txt")
mesc_df <- read.table(mesc_summary_file, header=TRUE, sep="\t")
# Joint data
df <- rbind(joint_ldsc_df, mesc_df)
df$permuted = rep(FALSE, dim(df)[1])

# Joint ldsc permuted
joint_ldsc_perm_summary_file <- paste0(visualize_trait_med_h2_dir, "med_h2_5_causal_tissue_permuted_eqtls_sim_results_joint_ldsc_binned_summary_averaged.txt")
joint_ldsc_perm_df <- read.table(joint_ldsc_perm_summary_file, header=TRUE, sep="\t")
# MESC permuted
mesc_perm_summary_file <- paste0(visualize_trait_med_h2_dir, "med_h2_5_causal_tissue_permuted_sim_results_mesc_summary_averaged.txt")
mesc_perm_df <- read.table(mesc_perm_summary_file, header=TRUE, sep="\t")
# Joint data
perm_df <- rbind(joint_ldsc_perm_df, mesc_perm_df)
perm_df$permuted = rep(TRUE, dim(df)[1])

# Compbine permuted and unpermuted data frames
big_df <- rbind(df, perm_df)
}


###################################################
# Joint plot showing mediated heritability, non-mediated heritability, and total heritablity for multiple methods
if (FALSE) {
# Non-mediated plot
nm_plot <- make_multimethod_se_barplot_for_nm_h2(df, c("mesc", "two_step_ldsc","joint_ldsc"))
# Mediated plot
med_plot <- make_multimethod_se_barplot_for_med_h2(df, c("mesc", "two_step_ldsc","joint_ldsc"))
# Total h2 plot
tot_plot <- make_multimethod_se_barplot_for_total_h2(df, c("mesc", "two_step_ldsc","joint_ldsc"))
# Get legend
legender <- get_legend(nm_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(med_plot+ theme(legend.position="none"), nm_plot+ theme(legend.position="none"), tot_plot+ theme(legend.position="none"),legender, ncol=1, rel_heights=c(1,1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, "multimethod_5_tissue_simulation_med_nm_total_h2_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.5, units="in")
print(output_file)
}

###################################################
# Joint plot showing mediated heritability in causal tissue and mediated heritability in non-causal tissues
if (FALSE) {
# Make plot for causal tissue
causal_med_plot <- make_multimethod_se_barplot_for_causal_tissue_med_h2(df, c("mesc", "two_step_ldsc","joint_ldsc"))
# Make plot for non-causal tissue
non_causal_med_plot <- make_multimethod_se_barplot_for_noncausal_tissue_med_h2(df, c("mesc", "two_step_ldsc","joint_ldsc"))
# Get legend
legender <- get_legend(non_causal_med_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(causal_med_plot+ theme(legend.position="none"), non_causal_med_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.2))

# Save
output_file <- paste0(visualize_trait_med_h2_dir, "multimethod_5_tissue_simulation_causal_tissue_med_vs_non_causal_tissue_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.5, units="in")
print(output_file)
}


###################################################
# Joint plot showing total mediated heritability in permuted and unpermuted data for each method
if (FALSE) {
# MESC
mesc_plot <- make_single_method_mediated_h2_se_barplot_for_permuted_and_unpermuted(big_df, "mesc")
# Two step ldsc
two_step_ldsc_plot <- make_single_method_mediated_h2_se_barplot_for_permuted_and_unpermuted(big_df, "two_step_ldsc")
# joint ldsc
joint_ldsc_plot <- make_single_method_mediated_h2_se_barplot_for_permuted_and_unpermuted(big_df, "joint_ldsc")

# Get legend
legender <- get_legend(joint_ldsc_plot + theme(legend.position="bottom"))
# Combine together with cowplot
joint <- plot_grid(mesc_plot+ theme(legend.position="none"), two_step_ldsc_plot+ theme(legend.position="none"), joint_ldsc_plot+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1, 1,.2))


# Save
output_file <- paste0(visualize_trait_med_h2_dir, "multimethod_5_tissue_simulation_permuted_med_vs_unpermuted_med_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.5, units="in")
print(output_file)
}



###################################################
if (FALSE) {
# Joint plot showing total permuted to unpermuted mediated heritability fraction
perm_med_fraction_plot <- make_multi_method_mediated_h2_se_barplot_for_permuted_and_unpermuted_ratio(big_df, c("mesc", "two_step_ldsc", "joint_ldsc"))
output_file <- paste0(visualize_trait_med_h2_dir, "multimethod_5_tissue_simulation_permuted_med_vs_unpermuted_med_ratio_barplot.pdf")
ggsave(perm_med_fraction_plot, file=output_file, width=7.2, height=4.5, units="in")
print(output_file)
}

















################
# OLD
###############


###############
# Single causal tissue
###############
# Load in data
if (FALSE) {
summary_file <- paste0(visualize_trait_med_h2_dir, "med_h2_single_causal_tissue_sim_results_summary_averaged.txt")
df <- read.table(summary_file, header=TRUE, sep="\t")



# Non-mediated plot
nm_plot <- make_multimethod_se_barplot_for_nm_h2(df, c("two_step_inference","joint_inference"))
# Mediated plot
med_plot <- make_multimethod_se_barplot_for_med_h2(df, c("two_step_inference","joint_inference"))
# Total h2 plot
tot_plot <- make_multimethod_se_barplot_for_total_h2(df, c("two_step_inference","joint_inference"))
# Combine together with cowplot
joint <- plot_grid(med_plot, nm_plot, tot_plot, ncol=1)
# Save
output_file <- paste0(visualize_trait_med_h2_dir, "multimethod_ldsc_single_causal_tissue_est_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.5, units="in")
}

if (FALSE) {
# Non-mediated plot
nm_plot <- make_se_barplot_for_nm_h2(df, "joint_inference")
# Mediated plot
med_plot <- make_se_barplot_for_med_h2(df, "joint_inference")
# Total h2 plot
tot_plot <- make_se_barplot_for_total_h2(df, "joint_inference")
# Combine together with cowplot
joint <- plot_grid(med_plot, nm_plot, tot_plot, ncol=1)
# Save
output_file <- paste0(visualize_trait_med_h2_dir, "joint_ldsc_single_causal_tissue_est_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.5, units="in")


# Non-mediated plot
nm_plot <- make_se_barplot_for_nm_h2(df, "two_step_inference")
# Mediated plot
med_plot <- make_se_barplot_for_med_h2(df, "two_step_inference")
# Total h2 plot
tot_plot <- make_se_barplot_for_total_h2(df, "two_step_inference")
# Combine together with cowplot
joint <- plot_grid(med_plot, nm_plot, tot_plot, ncol=1)
# Save
output_file <- paste0(visualize_trait_med_h2_dir, "two_step_ldsc_single_causal_tissue_est_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.5, units="in")
}

###############
# 5 tissues
###############
if (FALSE) {
# Load in data
summary_file <- paste0(visualize_trait_med_h2_dir, "med_h2_5_causal_tissue_sim_results_binned_summary_averaged.txt")
df <- read.table(summary_file, header=TRUE, sep="\t")


# Plot showing per tissue category mediated heritability
# Mediated plot
med_plot <- make_se_barplot_for_per_tissue_category_med_h2(df, "joint_inference_binned")
output_file <- paste0(visualize_trait_med_h2_dir, "joint_binned_ldsc_per_tissue_category_med_h2_est_se_barplot.pdf")
ggsave(med_plot, file=output_file, width=7.2, height=4.5, units="in")

med_plot <- make_se_barplot_for_per_tissue_category_med_h2(df, "two_step_inference")
output_file <- paste0(visualize_trait_med_h2_dir, "two_binned_step_ldsc_per_tissue_category_med_h2_est_se_barplot.pdf")
ggsave(med_plot, file=output_file, width=7.2, height=4.5, units="in")



# Plot showing per tissue mediated heritability
# Mediated plot
med_plot <- make_se_barplot_for_per_tissue_med_h2(df, "joint_inference_binned")
output_file <- paste0(visualize_trait_med_h2_dir, "joint_binned_ldsc_per_tissue_med_h2_est_se_barplot.pdf")
ggsave(med_plot, file=output_file, width=7.2, height=4.5, units="in")

med_plot <- make_se_barplot_for_per_tissue_med_h2(df, "two_step_inference")
output_file <- paste0(visualize_trait_med_h2_dir, "two_step_binned_ldsc_per_tissue_med_h2_est_se_barplot.pdf")
ggsave(med_plot, file=output_file, width=7.2, height=4.5, units="in")
}


if (FALSE) {
# Load in per simulation df
per_sim_summary_file <- paste0(visualize_trait_med_h2_dir, "med_h2_5_causal_tissue_sim_results_binned_summary_concatenated.txt")
per_sim_df <- read.table(per_sim_summary_file, sep="\t", header=TRUE)

# Compare distribution of tissue 0 scores
eqtl_sample_size="100"
est_distr_100 = make_histogram_showing_distribution_of_tissue_0_estimates(per_sim_df, eqtl_sample_size)

eqtl_sample_size="1000"
est_distr_1000 = make_histogram_showing_distribution_of_tissue_0_estimates(per_sim_df, eqtl_sample_size)

joint_histo <- plot_grid(est_distr_100, est_distr_1000, ncol=1)

output_file <- paste0(visualize_trait_med_h2_dir, "joint_binned_ldsc_tissue0_est_distribution_histo.pdf")
ggsave(joint_histo, file=output_file, width=7.2, height=4.5, units="in")
}





if (FALSE) {
# Load in data
summary_file <- paste0(visualize_trait_med_h2_dir, "med_h2_5_causal_tissue_sim_results_summary_averaged.txt")
df <- read.table(summary_file, header=TRUE, sep="\t")


# Plot showing per tissue category mediated heritability
# Mediated plot
med_plot <- make_se_barplot_for_per_tissue_category_med_h2(df, "joint_inference")
output_file <- paste0(visualize_trait_med_h2_dir, "joint_ldsc_per_tissue_category_med_h2_est_se_barplot.pdf")
ggsave(med_plot, file=output_file, width=7.2, height=4.5, units="in")

med_plot <- make_se_barplot_for_per_tissue_category_med_h2(df, "two_step_inference")
output_file <- paste0(visualize_trait_med_h2_dir, "two_step_ldsc_per_tissue_category_med_h2_est_se_barplot.pdf")
ggsave(med_plot, file=output_file, width=7.2, height=4.5, units="in")




# Plot showing per tissue mediated heritability
# Mediated plot
med_plot <- make_se_barplot_for_per_tissue_med_h2(df, "joint_inference")
output_file <- paste0(visualize_trait_med_h2_dir, "joint_ldsc_per_tissue_med_h2_est_se_barplot.pdf")
ggsave(med_plot, file=output_file, width=7.2, height=4.5, units="in")

med_plot <- make_se_barplot_for_per_tissue_med_h2(df, "two_step_inference")
output_file <- paste0(visualize_trait_med_h2_dir, "two_step_ldsc_per_tissue_med_h2_est_se_barplot.pdf")
ggsave(med_plot, file=output_file, width=7.2, height=4.5, units="in")
}

































##################
## Marginal PCA
##################
# Joint plot showing aggregate mediated, non-mediated, and total heritability
if (FALSE) {
# Non-mediated plot
nm_plot <- make_se_barplot_for_nm_h2(df, "marginal_pca_sumstat_gibbs_0.0")
# Mediated plot
med_plot <- make_se_barplot_for_med_h2(df, "marginal_pca_sumstat_gibbs_0.0")
# Total h2 plot
tot_plot <- make_se_barplot_for_total_h2(df, "marginal_pca_sumstat_gibbs_0.0")
# Combine together with cowplot
joint <- plot_grid(med_plot, nm_plot, tot_plot, ncol=1)
# Save
output_file <- paste0(visualize_trait_med_h2_dir, "marginal_pca_sumstat_gibbs_h2_est_se_barplot.pdf")
ggsave(joint, file=output_file, width=7.2, height=6.5, units="in")
}


if (FALSE) {
# Plot showing per tissue mediated heritability
# Mediated plot
med_plot <- make_se_barplot_for_per_tissue_med_h2(df, "marginal_pca_sumstat_gibbs_0.0")
output_file <- paste0(visualize_trait_med_h2_dir, "marginal_pca_sumstat_gibbs_per_tissue_med_h2_est_se_barplot.pdf")
ggsave(med_plot, file=output_file, width=7.2, height=4.5, units="in")


med_plot2 <- make_se_barplot_for_per_tissue_med_h2(df_eqtl_cov, "marginal_pca_sumstat_gibbs_0.0")
output_file <- paste0(visualize_trait_med_h2_dir, "marginal_pca_eqtl_cov_sumstat_gibbs_per_tissue_med_h2_est_se_barplot.pdf")
ggsave(med_plot2, file=output_file, width=7.2, height=4.5, units="in")


legender <- get_legend(med_plot + theme(legend.position="right") +guides(fill=guide_legend(nrow=5,byrow=TRUE)))
joint_plot <- plot_grid(plot_grid(med_plot + labs(title="No eQTL covariance") + theme(legend.position="none"), med_plot2 + labs(title="eQTL covariance")+ theme(legend.position="none"), ncol=1), legender, ncol=2, rel_widths=c(8,2.6))
output_file <- paste0(visualize_trait_med_h2_dir, "marginal_pca_eqtl_cov_sumstat_gibbs_per_tissue_med_h2_est_se_barplot_joint.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=4.5, units="in")
}


if (FALSE) {
# Plot showing per tissue mediated heritability
# Mediated plot
med_plot <- make_se_barplot_for_per_tissue_med_h2(df, "marginal_pca_sumstat_gibbs_0.0") + poster_figure_theme()
output_file <- paste0(visualize_trait_med_h2_dir, "marginal_pca_sumstat_gibbs_per_tissue_med_h2_est_se_barplot_for_poster.pdf")
ggsave(med_plot, file=output_file, width=16.2, height=7.7, units="in")
print(output_file)
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
