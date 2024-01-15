### Functions for enrichment analysis


### ---------------------------------------------------------------------------------------------------------

# This function, goplot, visualizes results of an enrichment analysis.
# This function takes three arguments: 
# data:           table containing results from an enrichment analysis.
# title:          string specifying the title of the plot.
# top:            integer specifying number of terms to plot on the y-axis.
# The function returns a ggplot visualizing the enrichment analysis. 
goplot <- function(data, 
                   title, 
                   top) {
  myplot <- data %>% 
    separate(col = Overlap, 
             into = c("Count",
                      "Total"), 
             sep ="/", 
             remove =FALSE) %>% 
    mutate(Count = as.numeric(Count),
           Gene.Ratio = Count/as.numeric(Total)) %>% 
    filter(Adjusted.P.value < 0.05) %>%
    slice_min(order_by = Adjusted.P.value,
              n = top,
              with_ties = FALSE) %>% 
    ggplot(data = .,
           mapping = aes(x = Gene.Ratio,
                         y = reorder(Term, 
                                     Gene.Ratio), 
                         color = Adjusted.P.value, 
                         size = Count)) +
    geom_point() +
    scale_y_discrete(labels = function(Term) str_wrap(Term, width = 20)) +
    #scale_x_continuous(breaks = seq(from = 0,
    #                                to = 1,
    #                                by = 0.1)) +
    scale_color_gradient(low = "red", 
                         high = "blue") +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 18,
                                     color = "black"),
          axis.title.x = element_text(size = 18,
                                      color = "black"),
          axis.text.y = element_text(size = 18,
                                     color = "black"),
          axis.ticks = element_line(size = 2),
          title = element_text(size = 12,
                               color = "black"),
          panel.border = element_rect(size = 2),
          legend.text = element_text(size = 18),
          legend.key.size = unit(1.3, "cm")) +
    labs(x = "Gene ratio",
         y = "",
         title = title)
  return(myplot)
}

### ---------------------------------------------------------------------------------------------------------






