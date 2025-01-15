library(tidyverse)
library(CEMiTool)
library(enrichR)
library(annotables)

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
    #filter(Adjusted.P.value < 0.05) %>%
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
    scale_y_discrete(labels = function(Term) str_wrap(Term, width = 50)) +
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


# Create vector of databases to use for enrichment analysis 
enrichr_dbs <- c("GO_Molecular_Function_2021",
                 "GO_Biological_Process_2021",
                 "MSigDB_Hallmark_2020")

# load data and metadata
data <- get(load('../processing/iCOPE_ALL_filtered_data.rda'))
meta <- get(load("../get_data/iCOPE_ALL_meta_data.rda"))

# first we turn ENSEMBL ID in data to HUGO names - this seems to be a requirement for CEMiTool
# and we need to do it anyway to understand what genes are present

# discard rows that are duplicate for ENSEMBL ID and HUGO but have different Entrez ID
unique_hugo <- annotables::grch38 %>% 
  select(c('ensgene', 'symbol')) %>% 
  distinct()

# merge to find corresponding HUGO names
data <- as_tibble(data, rownames='ensembl') %>%
  left_join(unique_hugo, by=join_by(ensembl==ensgene), relationship="one-to-one") %>%
  filter(symbol != "") %>%
  group_by(symbol) %>%
  filter(n() == 1) %>%
  ungroup()

# process data, dividing T and B subtypes
meta <- meta %>% select(c('ids', 'subtype'))

B_samples <- meta %>% filter(subtype=='B')
B_data <- data %>% as_tibble() %>% select(all_of(c(B_samples$ids, 'symbol'))) %>%
  column_to_rownames('symbol')
B_cem <- cemitool(B_data)
B_modules <- module_genes(B_cem)
write_csv(B_modules, file='modules_ALL_B.csv')

B_modules <- B_modules %>% filter(modules != "Not.Correlated")
for (module in unique(B_modules$modules)) {
  this_module <- B_modules |> filter(modules == module)
  enrich_analysis <- enrichr(genes = this_module$genes, 
                             databases = enrichr_dbs)
  for (db in enrichr_dbs) {
    plot <- goplot(enrich_analysis[[db]], paste0("B-ALL, ", module, ', ', db), 15)
    ggsave(filename = paste0("B_", module, "_", db, ".pdf"),
         plot = plot,
         width = 30,
         height = 15)
    ggsave(filename = paste0("B_", module, "_", db, ".png"),
           plot = plot,
           width = 30,
           height = 15)
    write_csv(enrich_analysis[[db]], file=paste0("B_", module, "_", db, ".csv"))
  }
}


T_samples <- meta %>% filter(subtype=='T')
T_data <- data %>% select(all_of(c(T_samples$ids, 'symbol'))) %>%
  column_to_rownames('symbol')
T_cem <- cemitool(T_data)
T_modules <- module_genes(T_cem)
write_csv(T_modules, file='modules_ALL_T.csv')

T_modules <- T_modules %>% filter(modules != "Not.Correlated")
for (module in unique(T_modules$modules)) {
  this_module <- T_modules |> filter(modules == module)
  enrich_analysis <- enrichr(genes = this_module$genes, 
                             databases = enrichr_dbs)
  for (db in enrichr_dbs) {
    plot <- goplot(enrich_analysis[[db]], paste0("T-ALL, ", module, ', ', db), 15)
    ggsave(filename = paste0("T_", module, "_", db, ".pdf"),
        plot = plot,
        width = 30,
        height = 15)
    ggsave(filename = paste0("T_", module, "_", db, ".png"),
           plot = plot,
           width = 30,
           height = 15)
    
    write_csv(enrich_analysis[[db]], file=paste0("T_", module, "_", db, ".csv"))
  }
}




