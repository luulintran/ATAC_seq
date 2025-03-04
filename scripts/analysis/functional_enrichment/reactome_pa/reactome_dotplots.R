library(ReactomePA)

# MAKE DOTPLOTS: ---------------------------------------------------------------

# Control_group:
filename <- file.path(output_dir_figures, 
                      paste0(control_group,"_reactome_dotplot.png"))

png(filename, width = 10, height = 10, units = "in", res = 300)

control_group_dotplot <- dotplot(get(paste0(control_group, "_reactome_pa")))
print(control_group_dotplot)

dev.off()


# Mutant_group: 
filename <- file.path(output_dir_figures, 
                      paste0(mutant_group,"_reactome_dotplot.png"))

png(filename, width = 10, height = 10, units = "in", res = 300)

mutant_group_dotplot <- dotplot(get(paste0(mutant_group, "_reactome_pa")))
print(mutant_group_dotplot)

dev.off()

print(paste0("Reactome PA results saved in ", output_dir_figures))

