library(ggplot2)

# MUTANT GROUP BARPLOT: --------------------------------------------------------
# Create the bar plot for mutant group
mutant_go_barplot <- ggplot(get(paste0(mutant_group, "_GO_terms")), 
                            aes(x = reorder(Description, log_p.adjust), 
                                y = log_p.adjust)) + 
  geom_bar(stat = "identity", fill = "#f48c67") + 
  coord_flip() +  # Flip coordinates to make it horizontal
  labs(title = "Upregulated GO Terms", 
       x = "GO Term", 
       y = "-log10(p.adjust)") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# SAVE PLOT: -------------------------------------------------------------------
filename <- file.path(output_dir_figures, paste0(mutant_group, "_upregulated_go_terms.pdf"))
pdf(filename, width = 5, height = 3)
print(mutant_go_barplot)

dev.off()

# CONTROL GROUP BARPLOT: --------------------------------------------------------
# Create the bar plot
control_go_barplot <- ggplot(get(paste0(control_group, "_GO_terms")), 
                            aes(x = reorder(Description, log_p.adjust), 
                                y = log_p.adjust)) + 
  geom_bar(stat = "identity", fill = "#ca8bc9") + 
  coord_flip() +  # make horizontal barplot
  labs(title = "Downregulated GO terms", 
       x = "GO Term", 
       y = "-log10(p.adjust)") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# SAVE PLOT: -------------------------------------------------------------------
filename <- file.path(output_dir_figures, paste0(control_group, "_upregulated_go_terms.pdf"))
pdf(filename, width = 5, height = 3)
print(control_go_barplot)

dev.off()

print(paste0("Bar plots for selected GO terms saved in ", output_dir_figures))