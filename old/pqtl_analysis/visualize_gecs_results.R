args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(dplyr)
library(reshape)
library(stringr)
library(reshape2)
library(ggbeeswarm)
library(RColorBrewer)
options(warn=1)
options(bitmapType='cairo')


figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}



make_bar_plot_showiing_average_mediated_h2_across_traits <- function(df) {
	# Initialize vectors to keep track of things
	method_vec <- c()
	anno_vec <- c()
	full_name_vec <- c()
	mean_vec <- c()
	mean_lb_vec <- c()
	mean_ub_vec <- c()

	method = "tglr"
	full_name = "mesc"
	anno = "baselineLD"
	tmp_df = df[as.character(df$method)==method,]
	tmp_df = tmp_df[as.character(tmp_df$non_mediated_annotations)==anno,]
	arr <- tmp_df$frac_h2_med
	meaner <- mean(arr)
	se_mean <- sd(arr)/sqrt(length(arr))
	meaner_ub = meaner + 1.96*se_mean
	meaner_lb = meaner - 1.96*se_mean
	method_vec <- c(method_vec, method)
	anno_vec <- c(anno_vec, anno)
	mean_vec <- c(mean_vec, meaner)
	mean_lb_vec <- c(mean_lb_vec, meaner_lb)
	mean_ub_vec <- c(mean_ub_vec, meaner_ub)
	full_name_vec <- c(full_name_vec, paste0(full_name, " : ", anno))

	method = "tglr"
	full_name = "mesc"
	anno = "genotype_intercept"
	tmp_df = df[as.character(df$method)==method,]
	tmp_df = tmp_df[as.character(tmp_df$non_mediated_annotations)==anno,]
	arr <- tmp_df$frac_h2_med
	meaner <- mean(arr)
	se_mean <- sd(arr)/sqrt(length(arr))
	meaner_ub = meaner + 1.96*se_mean
	meaner_lb = meaner - 1.96*se_mean
	method_vec <- c(method_vec, method)
	anno_vec <- c(anno_vec, anno)
	mean_vec <- c(mean_vec, meaner)
	mean_lb_vec <- c(mean_lb_vec, meaner_lb)
	mean_ub_vec <- c(mean_ub_vec, meaner_ub)
	full_name_vec <- c(full_name_vec, paste0(full_name, " : ", anno))

	method = "gecs"
	anno = "baselineLD"
	tmp_df = df[as.character(df$method)==method,]
	tmp_df = tmp_df[as.character(tmp_df$non_mediated_annotations)==anno,]
	arr <- tmp_df$frac_h2_med
	meaner <- mean(arr)
	se_mean <- sd(arr)/sqrt(length(arr))
	meaner_ub = meaner + 1.96*se_mean
	meaner_lb = meaner - 1.96*se_mean
	method_vec <- c(method_vec, method)
	anno_vec <- c(anno_vec, anno)
	mean_vec <- c(mean_vec, meaner)
	mean_lb_vec <- c(mean_lb_vec, meaner_lb)
	mean_ub_vec <- c(mean_ub_vec, meaner_ub)
	full_name_vec <- c(full_name_vec, paste0(method, " : ", anno))

	method = "gecs"
	anno = "genotype_intercept"
	tmp_df = df[as.character(df$method)==method,]
	tmp_df = tmp_df[as.character(tmp_df$non_mediated_annotations)==anno,]
	arr <- tmp_df$frac_h2_med
	meaner <- mean(arr)
	se_mean <- sd(arr)/sqrt(length(arr))
	meaner_ub = meaner + 1.96*se_mean
	meaner_lb = meaner - 1.96*se_mean
	method_vec <- c(method_vec, method)
	anno_vec <- c(anno_vec, anno)
	mean_vec <- c(mean_vec, meaner)
	mean_lb_vec <- c(mean_lb_vec, meaner_lb)
	mean_ub_vec <- c(mean_ub_vec, meaner_ub)
	full_name_vec <- c(full_name_vec, paste0(method, " : ", anno))


	final_df <- data.frame(method=method_vec, nm_anno=anno_vec, frac_h2_med=mean_vec, lb=mean_lb_vec, ub=mean_ub_vec, full_name=full_name_vec)

	final_df$full_name = factor(final_df$full_name, levels=c("mesc : genotype_intercept", "mesc : baselineLD", "gecs : genotype_intercept", "gecs : baselineLD"))

	p<-ggplot(final_df) +
  		geom_bar(aes(x=full_name, y=frac_h2_med),stat="identity",position=position_dodge(), fill="cyan4")+figure_theme() +
  		theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
  		theme(legend.position="bottom") +
  		geom_errorbar( aes(x=full_name, ymin=lb, ymax=ub), width=.3,colour="grey30", position=position_dodge()) +
  		labs(x="", y="Fraction h2 mediated")

  	return(p)
}


make_scatterplot_showing_med_h2_for_each_trait <- function(df, method1_name, method1_anno, method1_formal_name, method2_name, method2_anno,method2_formal_name, outlier_removal=FALSE) {
	tmp1_df = df[as.character(df$method)==method1_name,]
	tmp1_df = tmp1_df[as.character(tmp1_df$non_mediated_annotations)==method1_anno,]

	tmp2_df = df[as.character(df$method)==method2_name,]
	tmp2_df = tmp2_df[as.character(tmp2_df$non_mediated_annotations)==method2_anno,]


	final_df <- data.frame(method1_val=tmp1_df$frac_h2_med, method2_val=tmp2_df$frac_h2_med, method1_ub=tmp1_df$frac_h2_med_ub, method1_lb=tmp1_df$frac_h2_med_lb, method2_lb=tmp2_df$frac_h2_med_lb, method2_ub=tmp2_df$frac_h2_med_ub)

	if (outlier_removal) {
		final_df <- final_df[final_df$method1_val < .2,]
		final_df <- final_df[final_df$method2_val < .2,]

	}


	pp <- ggplot(final_df, aes(x=method1_val, y=method2_val)) + geom_point() +
	figure_theme() +
	labs(x=paste0(method1_formal_name, " : ", method1_anno), y=paste0(method2_formal_name, " : ", method2_anno)) +
	geom_abline(slope=1)

	return(pp)

}



gecs_results_dir <- args[1]
visualize_gecs_results_dir <- args[2]

summary_file <- paste0(gecs_results_dir, "gecs_results_summary.txt")

df <- read.table(summary_file,header=TRUE, sep='\t')

if (FALSE) {
# Make bar plot showing average mediated h2 across traits
bar_plot <- make_bar_plot_showiing_average_mediated_h2_across_traits(df)
output_file <- paste0(visualize_gecs_results_dir, "frac_h2_med_average_barplot.pdf")
ggsave(bar_plot, file=output_file, width=7.2, height=4.5, units="in")
print(output_file)
}

# Make scatter plot showing similarity of mediated h2 across traits
method1_name <- "tglr"
method1_formal_name <- "mesc"
method1_anno <- "genotype_intercept"
method2_name <- "tglr"
method2_formal_name <- "mesc"
method2_anno <- "baselineLD"
scatterplot <- make_scatterplot_showing_med_h2_for_each_trait(df, method1_name, method1_anno, method1_formal_name, method2_name, method2_anno,method2_formal_name)
output_file <- paste0(visualize_gecs_results_dir, "scatterplot_showing_per_trait_med_h2_",method1_name,"_", method1_anno,"_", method2_name,"_", method2_anno,".pdf")
ggsave(scatterplot, file=output_file, width=7.2, height=4.5, units="in")
print(output_file)


# Make scatter plot showing similarity of mediated h2 across traits
method1_name <- "gecs"
method1_formal_name <- "gecs"
method1_anno <- "genotype_intercept"
method2_name <- "gecs"
method2_formal_name <- "gecs"
method2_anno <- "baselineLD"
scatterplot <- make_scatterplot_showing_med_h2_for_each_trait(df, method1_name, method1_anno, method1_formal_name, method2_name, method2_anno,method2_formal_name)
output_file <- paste0(visualize_gecs_results_dir, "scatterplot_showing_per_trait_med_h2_",method1_name,"_", method1_anno,"_", method2_name,"_", method2_anno,".pdf")
ggsave(scatterplot, file=output_file, width=7.2, height=4.5, units="in")
print(output_file)

# Make scatter plot showing similarity of mediated h2 across traits
method1_name <- "tglr"
method1_formal_name <- "mesc"
method1_anno <- "baselineLD"
method2_name <- "gecs"
method2_formal_name <- "gecs"
method2_anno <- "baselineLD"
scatterplot <- make_scatterplot_showing_med_h2_for_each_trait(df, method1_name, method1_anno, method1_formal_name, method2_name, method2_anno,method2_formal_name)
output_file <- paste0(visualize_gecs_results_dir, "scatterplot_showing_per_trait_med_h2_",method1_name,"_", method1_anno,"_", method2_name,"_", method2_anno,".pdf")
ggsave(scatterplot, file=output_file, width=7.2, height=4.5, units="in")
print(output_file)



# Make scatter plot showing similarity of mediated h2 across traits
method1_name <- "tglr"
method1_formal_name <- "mesc"
method1_anno <- "baselineLD"
method2_name <- "gecs"
method2_formal_name <- "gecs"
method2_anno <- "baselineLD"
scatterplot <- make_scatterplot_showing_med_h2_for_each_trait(df, method1_name, method1_anno, method1_formal_name, method2_name, method2_anno,method2_formal_name, outlier_removal=TRUE)
output_file <- paste0(visualize_gecs_results_dir, "scatterplot_showing_per_trait_med_h2_",method1_name,"_", method1_anno,"_", method2_name,"_", method2_anno,"_outlier_removal.pdf")
ggsave(scatterplot, file=output_file, width=7.2, height=4.5, units="in")
print(output_file)


