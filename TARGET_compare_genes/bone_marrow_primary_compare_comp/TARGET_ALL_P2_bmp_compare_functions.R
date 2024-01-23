### This script contains functions used to compare results of various gene discovery analyses
### performed on TARGET-ALL-P2 gene expression data 


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
  comb_mat_sets <- make_comb_mat(sets_list,
                                 mode = "intersect")
  
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








