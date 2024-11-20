##### Plots to journal #####

# Join all figures
Figure_1 <- arrangeGrob(figure_1a,
                        figure_2a,
                        figure_1b,
                        figure_2b,
                        figure_1c,
                        figure_2c,
                        layout_matrix = rbind(c(1, 2),
                                              c(3, 4),
                                              c(3, 4),
                                              c(5, 6),
                                              c(5,6)))
# Adding letters
figure_one <- ggpubr::as_ggplot(Figure_1) +
  draw_plot_label(label = LETTERS[1:6],
                  x = c(0, 0.5, 0, 0.5, 0, 0.5),
                  y = c(.99, .99, .790, .790, .350, .350))
# (*.pdf) file
ggsave(filename = "PCA_Journal_Plot/Figure_1_PCA_Journal_new.pdf",
       plot = figure_one, width = 150, height = 180,
       units = "mm", dpi = 300, scale = 2.5)
# (*.png) file
ggsave(filename = "PCA_Journal_Plot/Figure_1_PCA_Journal_new.png",
       plot = figure_one, width = 150, height = 180,
       units = "mm", dpi = 300, scale = 2.5)
# (*.jpg) file
ggsave(filename = "PCA_Journal_Plot/Figure_1_PCA_Journal_new.jpg",
       plot = figure_one, width = 150, height = 180,
       units = "mm", dpi = 300, scale = 2.5)
