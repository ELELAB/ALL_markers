### This script contains functions used for visualizing DEA results of TARGET-ALL-P2 gene expression data 


### ---------------------------------------------------------------------------------------------------------

### This function, upset_plot, creates an UpSet plot of intersections between
### different sets. This function takes five inputs:
### file_name:      File name of UpSet plot to be saved as a pdf
### dir_output:     Directory where UpSet plot should be saved
### sets_list:      A list containing those sets that are to be compared 
### names_sets:     A character vector containing names of each element in the list, sets_list
### title_plot:     The title to be included in the UpSet plot
### The function returns an UpSet plot. 
upset_plot <- function(file_name,
                       dir_output,
                       sets_list, 
                       names_sets, 
                       title_plot) {
  
  # Save UpSet plot as pdf
  pdf(file = paste0(dir_output, 
                    "/", 
                    file_name),
      width = 18,
      height = 8)
  
  # Add names to elements in list of sets
  names(sets_list) <- names_sets
  
  # Make combination matrix of sets to be used in UpSet plot
  comb_mat_sets <- make_comb_mat(sets_list)
  
  # Generate color palette for plot using viridis package
  n_col <- max(comb_degree(comb_mat_sets))
  palette_col <- viridis_pal(option = "viridis")(n_col)
  
  # Create UpSet plot
  upset_p <- UpSet(comb_mat_sets, 
                   set_order = names(sets_list),
                   pt_size = unit(5, 
                                  "mm"), 
                   lwd = 3, 
                   height = unit(4, 
                                 "cm"),
                   comb_col = palette_col[comb_degree(comb_mat_sets)],
                   top_annotation = upset_top_annotation(comb_mat_sets, 
                                                         height = unit(12, 
                                                                       "cm"),
                                                         ylim = c(0, 
                                                                  max(comb_size(comb_mat_sets))),
                                                         bar_width = 0.7, 
                                                         axis_param = list(side = "left", 
                                                                           at = seq(from = 0,
                                                                                    to = max(comb_size(comb_mat_sets)),
                                                                                    by = 500)),
                                                         annotation_name_side = "left", 
                                                         annotation_name_gp = gpar(cex = 1), 
                                                         annotation_name_offset = unit(1.5,
                                                                                       "cm")),
                   right_annotation = upset_right_annotation(comb_mat_sets, 
                                                             width = unit(3, 
                                                                          "cm"), 
                                                             gp = gpar(fill = "darkseagreen"),
                                                             axis_param = list(at = seq(from = 0,
                                                                                        to = max(set_size(comb_mat_sets)),
                                                                                        by = 2000)), 
                                                             annotation_name_offset = unit(1.5, 
                                                                                           "cm")),
                   row_names_gp = gpar(fontsize = 12))
  
  # Add number of elements in each set on top of bars in plot
  draw_upset <- draw(upset_p)
  col_ord <- column_order(draw_upset)
  c_s <- comb_size(comb_mat_sets)
  decorate_annotation("intersection_size", {
    grid.text(c_s[col_ord], 
              x = seq(c_s), 
              y = unit(c_s[col_ord], 
                       "native") + 
                unit(2, "pt"), 
              gp = gpar(fontsize = 12, 
                        fontface = "bold"),
              just = "bottom",
              default.units = "native")
  })
  
  # Add title to plot
  grid.text(label = title_plot, 
            x = unit(20, 
                     "cm"), 
            y = unit(18, 
                     "cm"), 
            gp = gpar(fontsize = 18),
            just = "centre")
  
  # End with dev.off() to save plot
  dev.off()
  
}

### ---------------------------------------------------------------------------------------------------------


### ---------------------------------------------------------------------------------------------------------

### This function, venn_diagram, creates a Venn diagram of intersections between
### different sets. This function takes five inputs:
### sets_list:        A list containing those sets that are to be compared 
### names_sets:       A character vector containing names of each element in the list, sets_list
### title_plot:       The title to be included in the Venn diagram
### dir_output:       The path to the directory where the plot should be saved
### file_name:        File name of plot that will be saved; best if saved as .png
### The function returns a Venn diagram. 
venn_diagram <- function(sets_list,
                         names_sets,
                         title_plot,
                         dir_output,
                         file_name) {
  
  # Assign names to list of sets
  names(sets_list) <- names_sets
  
  # Define colors to use in Venn diagram
  n_col <- length(sets_list)
  palette_col <- viridis_pal(alpha = 0.8,
                             option = "viridis")(n_col)
  
  #Create Venn diagram and save it as pdf in plots directory of TCGA project
  venn.diagram(sets_list, 
               filename = paste0(dir_output, 
                                 "/", 
                                 file_name), 
               category.names = names_sets, 
               output = TRUE, 
               main = title_plot, 
               main.cex = 1.5, 
               main.fontface = "bold", 
               main.fontfamily = "sans", 
               lwd = 3, 
               lty = "blank", 
               fill = palette_col, 
               cex = 2, 
               fontface = "bold", 
               fontfamily = "sans",
               cat.cex = 1.5, 
               cat.fontface = "bold", 
               cat.fontfamily = "sans")
  
}

### ---------------------------------------------------------------------------------------------------------


### ---------------------------------------------------------------------------------------------------------

### This function, volcano_plot, creates a volcano plot of differentially expressed
### genes. This function takes ten inputs:
### highlight_genes:    A boolean value either TRUE or FALSE indicating if certain
###                     genes should be highlight in the plot
### highlight_up:       A character vector of upregulated genes which should be 
###                     highlighted in the plot
### highlight_down:     A character vector of downregulated which should be
###                     highlighted in the plot
### DEG:                A vector of DEGs
### logFC:              A vector of logFC values of DEGs
### p_value:            A vector of adjusted FDR values of DEGs
### p_value_cut:        An integer indicating cut-off of FDR values
### title_plot:         The title to be included in the Venn diagram
### dir_output:         The path to the directory where the plot should be saved
### file_name:          File name of plot that will be saved; best if saved as .png
### The function returns an Volcano plot. 
volcano_plot <- function(highlight_genes,
                         highlight_up,
                         highlight_down,
                         DEG,
                         logFC,
                         p_value,
                         p_value_cut,
                         title_plot,
                         dir_output,
                         file_name) {
  
  if (highlight_genes == TRUE) {
    
    # Create volcano plot and save plot in specified directory
    TCGAVisualize_volcano(x = logFC, 
                          y = p_value, 
                          filename = paste0(dir_output, 
                                            "/", 
                                            file_name), 
                          ylab = "-Log[10](FDR corrected -P values)", 
                          xlab = "Log2 Fold Change", 
                          y.cut = p_value_cut,
                          title = title_plot, 
                          legend =  NULL, 
                          label = c("Not significant", 
                                    "Upregulated genes", 
                                    "Downregulated genes"), 
                          color = c("black", 
                                    "red", 
                                    "green"),
                          names = DEG, 
                          names.fill = TRUE, 
                          show.names = "highlighted", 
                          highlight = c(highlight_up, 
                                        highlight_down), 
                          highlight.color = c("orange"))
    
  } else {
    
    # Create volcano plot and save plot in specified directory
    TCGAVisualize_volcano(x = logFC, 
                          y = p_value, 
                          filename = paste0(dir_output, 
                                            "/", 
                                            file_name), 
                          ylab = "-Log[10](FDR corrected -P values)", 
                          xlab = "Log2 Fold Change", 
                          y.cut = p_value_cut,
                          title = title_plot, 
                          legend =  NULL, 
                          label = c("Not significant", 
                                    "Upregulated genes", 
                                    "Downregulated genes"), 
                          color = c("black", 
                                    "red", 
                                    "green"),
                          names.fill = FALSE)
    
  }
  
}

### ---------------------------------------------------------------------------------------------------------


