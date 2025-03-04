
library(clusterProfiler)


# MAKE BARPLOTS: ---------------------------------------------------------------
filename <- file.path(output_dir_figures, 
                      paste0(control_group,"_GO_results_barplot.png"))

png(filename, width = 10, height = 10, units = "in", res = 300)

# barplot for GO term enrichment
control_group_barplot <- barplot(get(paste0(control_group, "_GO_results")), showCategory = 20)
print(control_group_barplot)

dev.off()

filename <- file.path(output_dir_figures, 
                      paste0(mutant_group,"_GO_results_barplot.png"))

png(filename, width = 10, height = 10, units = "in", res = 300)

# barplot for GO term enrichment
mutant_group_barplot <- barplot(get(paste0(mutant_group, "_GO_results")), showCategory = 20)
print(mutant_group_barplot)

dev.off()

print(paste0("GO results saved in ", output_dir_figures))