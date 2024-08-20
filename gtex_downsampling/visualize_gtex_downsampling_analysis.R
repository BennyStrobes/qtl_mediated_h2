args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}








scatter_showing_med_h2_with_sample_size <- function(df, title="") {
	pp <- ggplot(df, aes(y=med_h2, x=sample_size)) +
	      geom_errorbar(aes(ymin = med_h2_lb, ymax=med_h2_ub), width = 0.2, color='grey45') +
  		geom_point(size=2.0, color="dodgerblue") +
  		figure_theme() +
  		labs(y="Fraction h2 mediated by predicted expression", x="Aggregate GTEx sample size", title=title) +
  		theme(legend.position="bottom")

  	return(pp)
}




scatter_showing_med_h2_with_sample_size_unweighted <- function(df, title="") {
	df$lb = df$unweighted_mean - (1.96*df$med_h2_se)
	df$ub = df$unweighted_mean + (1.96*df$med_h2_se)

	print(cor.test(df$sample_size, df$unweighted_mean))
	pp <- ggplot(df, aes(y=unweighted_mean, x=sample_size)) +
  		geom_point(size=2.0, color="dodgerblue") +
  		geom_errorbar(aes(ymin = lb, ymax=ub), width = 0.2, color='grey45') +
  		figure_theme() +
  		labs(y="Fraction h2 mediated by predicted expression", x="Aggregate GTEx sample size", title=title) +
  		theme(legend.position="bottom")

  	return(pp)
}



scatter_showing_mesc_med_h2_vs_tglr_med_h2 <- function(mesc_df, tglr_df, title="") {
	df <- data.frame(mesc=mesc_df$med_h2, mesc_lb=mesc_df$med_h2-(1.96*mesc_df$med_h2_se), mesc_ub=mesc_df$med_h2+(1.96*mesc_df$med_h2_se), tglr=tglr_df$med_h2, tglr_lb=tglr_df$med_h2-(1.96*tglr_df$med_h2_se), tglr_ub=tglr_df$med_h2+(1.96*tglr_df$med_h2_se))

	pp <- ggplot(df, aes(y=mesc, x=tglr)) +
	      geom_errorbar(aes(ymin = mesc_lb, ymax=mesc_ub), color='grey55') +
	      geom_errorbar(aes(xmin=tglr_lb, xmax=tglr_ub), color='grey55') +
  		geom_point(size=2.3, color="dodgerblue") +
  		figure_theme() +
  		labs(y="MESC Fraction h2 mediated by predicted expression", x="TGLR Fraction h2 mediated by predicted expression", title=title) +
  		theme(legend.position="bottom") +
  		geom_abline(slope=1)


  	return(pp)

}

scatter_showing_mesc_med_h2_vs_mesc_with_other_sample_size <- function(mesc_df, tglr_df, title="") {
	df <- data.frame(mesc=mesc_df$med_h2, mesc_lb=mesc_df$med_h2-(1.96*mesc_df$med_h2_se), mesc_ub=mesc_df$med_h2+(1.96*mesc_df$med_h2_se), tglr=tglr_df$med_h2, tglr_lb=tglr_df$med_h2-(1.96*tglr_df$med_h2_se), tglr_ub=tglr_df$med_h2+(1.96*tglr_df$med_h2_se))

	pp <- ggplot(df, aes(y=mesc, x=tglr)) +
	      geom_errorbar(aes(ymin = mesc_lb, ymax=mesc_ub), color='grey55') +
	      geom_errorbar(aes(xmin=tglr_lb, xmax=tglr_ub), color='grey55') +
  		geom_point(size=2.3, color="dodgerblue") +
  		figure_theme() +
  		labs(y="MESC (full data) Fraction h2 mediated", x="MESC (half data) Fraction h2 mediated", title=title) +
  		theme(legend.position="bottom") +
  		geom_abline(slope=1)


  	return(pp)

}


make_se_barplot_of_mesc_results_across_traits <- function(mesc_df, ordered_traits,title="") {

	df <- data.frame(trait_name=mesc_df$trait_name, mesc=mesc_df$med_h2, mesc_lb=mesc_df$med_h2-(1.96*mesc_df$med_h2_se), mesc_ub=mesc_df$med_h2+(1.96*mesc_df$med_h2_se))


	reordered_traits = ordered_traits[order(df$mesc)]

	df$trait_name = factor(ordered_traits, levels=reordered_traits)






	p<-ggplot(data=df, aes(x=trait_name, y=mesc)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=mesc_lb, ymax=mesc_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9.5)) +
  		labs(x="", y="Fraction h2 med.", title=title)   	
  	return(p)

}


make_se_barplot_of_correlations_between_gene_ld_scores_and_variant_annotations <- function(df, threshold=FALSE, font_size=8.5) {

	df <- data.frame(variant_anno=df$variant_anno, correlation=df$correlation, correlation_lb=df$correlation-(1.96*df$correlation_se), correlation_ub=df$correlation+(1.96*df$correlation_se))

	if (threshold==TRUE) {
		df <- df[df$correlation > 0.4,]
	}

	ordered_anno = as.character(df$variant_anno)
	reordered_anno = ordered_anno[order(df$correlation)]


	df$variant_anno = factor(ordered_anno, levels=reordered_anno)

	p<-ggplot(data=df, aes(x=variant_anno, y=correlation)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=correlation_lb, ymax=correlation_ub), width=.4, position=position_dodge(.9))  +
  		figure_theme() +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=font_size)) +
  		labs(x="", y="Correlation")   	
  	return(p)



}







########################
# Command line args
########################
mesc_results_dir= args[1]
tglr_results_dir=args[2]
gene_ldscore_variant_ld_score_corr_dir=args[3]
visualize_gtex_downsampling_dir=args[4]


# Create mapping from trait identifier to trait_name readable
trait_name_readable <- "/n/scratch/users/b/bes710/tmp_storage/ldsc/sumstats_formatted_2024/non_redundent_traits_EUR_2024_readable.txt"
trait_name_readable_df <- read.table(trait_name_readable, header=TRUE, sep="\t")
trait_name_readable <- as.character(trait_name_readable_df$trait_identifier)



###########################
# Scatter plot comparing mesc med-h2 as a function of sample size
###########################
if (FALSE) {
# Load in data for genotype intercept
anno_version="genotype_intercept"
mesc_res_file = paste0(mesc_results_dir, "downsampled_med_h2_meta_analyzed_summary_", anno_version, ".txt")
df <- read.table(mesc_res_file, header=TRUE)

scatter <- scatter_showing_med_h2_with_sample_size(df, title=anno_version)
output_file <- paste0(visualize_gtex_downsampling_dir, "mesc_med_h2_vs_sample_size_", anno_version, "_scatter.pdf")
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")

scatter <- scatter_showing_med_h2_with_sample_size_unweighted(df, title=anno_version)
output_file <- paste0(visualize_gtex_downsampling_dir, "mesc_med_h2_vs_sample_size_", anno_version, "_unweighted_scatter.pdf")
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")
}

print("HELLO")
# Load in data for BaselineLD
anno_version="baselineLD_no_qtl"
mesc_res_file = paste0(mesc_results_dir, "downsampled_med_h2_meta_analyzed_summary_", anno_version, ".txt")
df <- read.table(mesc_res_file, header=TRUE)

scatter <- scatter_showing_med_h2_with_sample_size(df, title=anno_version)
output_file <- paste0(visualize_gtex_downsampling_dir, "mesc_med_h2_vs_sample_size_", anno_version, "_scatter.pdf")
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")


scatter <- scatter_showing_med_h2_with_sample_size_unweighted(df, title=anno_version)
output_file <- paste0(visualize_gtex_downsampling_dir, "mesc_med_h2_vs_sample_size_", anno_version, "_unweighted_scatter.pdf")
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")


###########################
# Bar plot showing full data analysis MESC med-h2 results per trait
###########################
if (FALSE) {
# Load in data for baselineLD intercept
anno_version="baselineLD_no_qtl"
mesc_file <- paste0(mesc_results_dir, "cross_trait_med_h2_summary_full_data_", anno_version, ".txt")
mesc_df <- read.table(mesc_file, header=TRUE, sep="\t")

se_barplot <- make_se_barplot_of_mesc_results_across_traits(mesc_df, trait_name_readable, title=anno_version)

output_file <- paste0(visualize_gtex_downsampling_dir, "mesc_med_h2_across_traits_", anno_version, "_se_barplot.pdf")
ggsave(se_barplot, file=output_file, width=7.2, height=5.5, units="in")
}

###########################
# Scatter plot (across traits) comparing mesc (full sample size) and MESC (half sample size)
###########################
if (FALSE) {
# Load in data for baselineLD intercept
anno_version="baselineLD_no_qtl"
mesc_file <- paste0(mesc_results_dir, "cross_trait_med_h2_summary_full_data_", anno_version, ".txt")
mesc_full_df <- read.table(mesc_file, header=TRUE, sep="\t")
mesc_file <- paste0(mesc_results_dir, "cross_trait_med_h2_summary_downsample_6797_", anno_version, ".txt")
mesc_half_df <- read.table(mesc_file, header=TRUE, sep="\t")

output_file <- paste0(visualize_gtex_downsampling_dir, "mesc_full_data_mesc_half_data_scatter_comparison_", anno_version, ".pdf")
scatter <- scatter_showing_mesc_med_h2_vs_mesc_with_other_sample_size(mesc_full_df, mesc_half_df, title=anno_version)
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")
}



###########################
# Scatter plot (across traits) comparing mesc and TGLR at full sample size
###########################
if (FALSE) {
# Load in data for genotype intercept
anno_version="genotype_intercept"
mesc_file <- paste0(mesc_results_dir, "cross_trait_med_h2_summary_full_data_", anno_version, ".txt")
mesc_df <- read.table(mesc_file, header=TRUE, sep="\t")

tglr_file <- paste0(tglr_results_dir, "cross_trait_med_h2_summary_full_data_", anno_version, ".txt")
tglr_df <- read.table(tglr_file, header=TRUE, sep="\t")

output_file <- paste0(visualize_gtex_downsampling_dir, "mesc_tglr_scatter_comparison_", anno_version, ".pdf")
scatter <- scatter_showing_mesc_med_h2_vs_tglr_med_h2(mesc_df, tglr_df, title=anno_version)
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")

# Load in data for baselineLD intercept
anno_version="baselineLD_no_qtl"
mesc_file <- paste0(mesc_results_dir, "cross_trait_med_h2_summary_full_data_", anno_version, ".txt")
mesc_df <- read.table(mesc_file, header=TRUE, sep="\t")

tglr_file <- paste0(tglr_results_dir, "cross_trait_med_h2_summary_full_data_", anno_version, ".txt")
tglr_df <- read.table(tglr_file, header=TRUE, sep="\t")

output_file <- paste0(visualize_gtex_downsampling_dir, "mesc_tglr_scatter_comparison_", anno_version, ".pdf")
scatter <- scatter_showing_mesc_med_h2_vs_tglr_med_h2(mesc_df, tglr_df, title=anno_version)
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")
}




###########################
# Standard error bar plot showing correlation between gene ld scores and variant ld scores
###########################
if (FALSE) {
anno_version="baselineLD_no_qtl"
anno_corr_file <- paste0(gene_ldscore_variant_ld_score_corr_dir, "full_datagene_ld_score_variant_ld_score_correlations.txt")
anno_corr_df <- read.table(anno_corr_file, header=TRUE,sep="\t")
anno_corr_df = anno_corr_df[as.character(anno_corr_df$gene_anno) == "aggregate",]

se_barplot <- make_se_barplot_of_correlations_between_gene_ld_scores_and_variant_annotations(anno_corr_df)
output_file <- paste0(visualize_gtex_downsampling_dir, "corr_gene_ld_score_variant_ld_score_se_barplot.pdf")
ggsave(se_barplot, file=output_file, width=9.2, height=5.5, units="in")

se_barplot <- make_se_barplot_of_correlations_between_gene_ld_scores_and_variant_annotations(anno_corr_df, threshold=TRUE,font_size=11)
output_file <- paste0(visualize_gtex_downsampling_dir, "corr_gene_ld_score_variant_ld_score_se_barplot2.pdf")
ggsave(se_barplot, file=output_file, width=7.2, height=5.5, units="in")
}
