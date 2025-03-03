# EXPLORATORY ANALYSIS: --------------------------------------------------------

# Venn Diagram
filename <- file.path(output_dir_figures,"venn_diagram.png")

png(filename, width = 5, height = 5, units = "in", res = 300)

dba.plotVenn(dbObj, contrast = 1, bDB = TRUE,
             bGain = TRUE, bLoss = TRUE, bAll = FALSE)

dev.off()

# PCA Plot ******************************************************

# PCA on all sites
filename <- file.path(output_dir_figures,"pca_plot_allsites.png")

png(filename, width = 5, height = 5, units = "in", res = 300)

dba.plotPCA(dbObj, DBA_CONDITION, label = DBA_ID)

dev.off()

# PCA on differential sites
filename <- file.path(output_dir_figures,"pca_plot_diffsites.png")
png(filename, width = 5, height = 5, units = "in", res = 300)
dba.plotPCA(dbObj, contrast = 1, label = DBA_CONDITION)

dev.off()

# MA plot ******************************************************
filename <- file.path(output_dir_figures,"ma_plot.png")

png(filename, width = 5, height = 5, units = "in", res = 300)

dba.plotMA(dbObj)

dev.off()

# Volcano plot ******************************************************
filename <- file.path(output_dir_figures,"volcanoplot.png")

png(filename, width = 5, height = 5, units = "in", res = 300)

dba.plotVolcano(dbObj)

dev.off()

# Box plot ******************************************************
filename <- file.path(output_dir_figures,"boxplot.png")
png(filename, width = 5, height = 5, units = "in", res = 300)
dba.plotBox(dbObj)

dev.off()

# Heatmap  ******************************************************

# Heatmap colors
filename <- file.path(output_dir_figures,"heatmap_diffsites.png")
png(filename, width = 5, height = 5, units = "in", res = 300)

hmap <- colorRampPalette(c("white", "#938fc2", "#40027e")) (n =100)
readscores <- dba.plotHeatmap(dbObj, contrast = 1, correlations = FALSE,
                              scale = "row", colScheme = hmap)
print(readscores)

dev.off()

# Profile plot ******************************************************
filename <- file.path(output_dir_figures,"profile_plot.png")
png(filename, width = 5, height = 5, units = "in", res = 300)

profiles <- dba.plotProfile(dbObj, merge = DBA_REPLICATE)
dba.plotProfile(profiles)

dev.off()

print(paste0("Exploratory analysis figures saved in ", output_dir_figures))