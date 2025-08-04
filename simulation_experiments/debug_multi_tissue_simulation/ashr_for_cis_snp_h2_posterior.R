args = commandArgs(trailingOnly=TRUE)
library(ashr)







input_file <- args[1]
output_file <- args[2]



df <- read.table(input_file, header=TRUE, sep="\t")

meany = df$cis_snp_h2
mean_se = df$cis_snp_h2_se

ash_obj = ash(meany, mean_se, mixcompdist="halfuniform")

post_means = get_pm(ash_obj)

write.table(post_means, file=output_file, quote=FALSE, col.names=FALSE, row.names=FALSE)