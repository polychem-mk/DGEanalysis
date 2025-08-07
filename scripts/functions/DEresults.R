# function DEresults(), returns a list of 3:

#  [[1]] Table with top n differentially expressed genes.
#  [[2]] Summary and keywords for these genes. 
#  [[3]] ggplot for normalized counts filtered for the specified batch (or all
#        batches if ‘batch=NULL’).

# The DEresults() function shows the results for all batches. It does not rerun
# the DESeq() function, but only retrieves the results. However, the data for
# the normalized counts plot can be filtered using the ‘batch’ argument.

# Arguments:

# de        output of DESeq() function.
# md        metadata.
# level.1 and level.2 are the base level and the comparison level in the factor 
#           variable 'cancer_type', they are provided as arguments to 'contrast'
#           in the results() function.
# alpha     is cutoff for padj to filter results.

# sort.by   can be "stat", "padj" or "log2FoldChange". The results will be 
#           sorted by 'sort.by'.

# batch     is the batch number, if NULL, then the results of all batches are
#           displayed on the counts plot.
# n         number of top results in the results table and on the counts plot.


DEresults <- function(de, 
                      md,
                      level.1 = "HC",
                      level.2 = NULL,
                      alpha = 0.01,
                      sort.by = "stat", 
                      batch = NULL,
                      n = 20){
  
  # set batch to character if it is not NULL
  if(!is.null(batch)){
    batch_i  = as.character(batch)
  }else{
    batch_i  = NULL
  }
  
  # Extract results for level.1 and level.2 contrasts
  res = results(de,             
                contrast = c("cancer_type", level.2, level.1),
                alpha = alpha)
  
  # Turn 'res' into data frame and filter results by alpha
  res = res %>%                
    as.data.frame() %>%
    filter(padj <= alpha)  
  
  # Sort results by 'sort.by' argument
  if(sort.by == "stat"){
    res = res %>%
      arrange(desc(abs(stat))) 
  }else if(sort.by == "log2FoldChange"){
    res = res %>%
      arrange(desc(abs(log2FoldChange))) 
  }else{
    res = res %>%
      arrange(padj) 
  }
  
   # Add annotation 
  res$ENSEMBL = rownames(res)
  
  suppressMessages({
    ann_data = AnnotationDbi::select(org.Hs.eg.db, 
                                     keys = res$ENSEMBL, 
                                     columns = c("SYMBOL", "GENENAME", "GENETYPE"), 
                                     keytype = "ENSEMBL") %>%
      # org.Hs.eg.db may have more than one row per ensemble ID, so keep only
      # the first match.
      distinct(ENSEMBL, .keep_all = TRUE)
  })
  
  # table 'res' contains results for 'padj' greater than 'alpha', and is sorted 
  # by the 'sort.by' argument; we add the gene name, symbol, and type; if these 
  # values are missing, we replace the gene symbol with the ensemble ID.
  # The variable 'dir' indicates whether the gene is up or down regulated.
  res = res %>%
    left_join(ann_data, by = "ENSEMBL") %>%
    mutate(SYMBOL = ifelse(is.na(SYMBOL),ENSEMBL , SYMBOL),
           dir = ifelse(log2FoldChange < 0, "down", "up"))
  
  # Keywords _______________________________
  
  # The DEresults() function shows 20 keywords, 10 of which are from the top
  # n gene names, and the remaining 10 are from all results.
  
  n_res = nrow(res)  # number of results in 'res'
  
  # percentage of down regulated genes in 'res':
  dr = round(100*sum(res$dir == "down", na.rm = TRUE)/n_res, 2)
  
  # Summary showing the number of differentially expressed genes with 'padj'
  # less than 'alpha' and the percentages of down- and up-regulated genes:
  genes_summary = paste("There are ", n_res,
                        " differentially expressed genes, ",
                        dr, "% are down-regulated and ", 
                        round(100-dr, 2), "% are up-regulated.",
                         sep = "")
  
  # Table with the number of up and down regulated genes by gene type
  gene_types = table(res$GENETYPE, res$dir) %>%
    as.data.frame( ) %>%
    pivot_wider(names_from = Var2, values_from = Freq)
  
  colnames(gene_types)[1] = "gene.type"
  
  # words to be excluded from keywords
  stop_words = c("protein", "and", "rich", "like", "subunit", "repeat",
                 "motif", "domain", "associated", "interacting", "family",
                 "alpha", "beta", "A", "B", "component", "small", "solute",
                 "variable", "molecule", "containing", "member", "of" ,
                 "homolog", "related", "constant", "acid") 
  
  # get all keywords from table 'res', add word count and sort by count in
  # descending order
  keywords_all = res$GENENAME %>%
    str_split(" ") %>% unlist() %>% table() %>% 
    sort(decreasing = TRUE)  %>% names()
  
  # the same for top n results
  keywords_top = res$GENENAME[1:n] %>%
    str_split(" ") %>% unlist() %>% table() %>% 
    sort(decreasing = TRUE)  %>% names()
  
  # filter out 'stop words'
  keywords_top = setdiff(keywords_top[str_detect(keywords_top, "\\d",
                                                 negate = TRUE)], 
                     stop_words)
  keywords_all = setdiff(keywords_all[str_detect(keywords_all, "\\d",
                                                 negate = TRUE)], 
                     c(stop_words, keywords_top[1:10]))
  
  # 'kw' list contains:
  #  1 genes_summary    string, the number of differentially expressed genes  
  #                     and the percentages of down- and up-regulated genes;
  #  2 gene_types       data frame with the number of up and down regulated 
  #                     genes by gene type;
  #  3  string with 10 keywords from top n results and 10 keywords from all
  #                     results
  kw = list(genes_summary, gene_types, c(keywords_top[1:10], keywords_all[1:10]))
  
   # Top n ____________________________
  
  # Filter for the top n results
  res = res %>% 
    slice_head(n = n)
  
  # Data for the normalized counts plot:
  res_plot = lapply(res$ENSEMBL, function(i){
    
    # get data with plotCounts() function, 
    # add butch and mutational_subclass from metadata,
    # filter data for control samples and cancer_type specified in level.2 argument
    # add Enseml ID:
    
   if(!is.null(batch_i)){
     # If a batch is specified, filter the data by batch.
      res_i = plotCounts(de,
                         i, 
                         intgroup = "cancer_type",
                         returnData = TRUE)  %>%
        bind_cols(batch = md$batch,
                  mutational_subclass = md$mutational_subclass) %>%
        filter(cancer_type %in% c(level.2, level.1) ) %>%
        filter(batch == batch_i) %>%
        mutate(ENSEMBL = i)
      
    }else{
      res_i = plotCounts(de,
                         i, 
                         intgroup = "cancer_type",
                         returnData = TRUE)  %>%
        bind_cols(batch = md$batch,
                  mutational_subclass = md$mutational_subclass) %>%
        filter(cancer_type %in% c(level.2, level.1) ) %>%
        mutate(ENSEMBL = i)
    }
    
    res_i
    
  }) %>% bind_rows()
  
  # instead of the ensemble IDs we will use gene symbols
  res_plot = res_plot %>%
    left_join(res[ , c("ENSEMBL", "SYMBOL")], by = "ENSEMBL")
  
  # Set the plot title
  if(!is.null(batch)){
    plot_title = paste("Top ", n, 
                       " differentially expressed genes, batch ", 
                       batch_i, sep = "")
  }else{
    plot_title = paste("Top ", n, 
                       " differentially expressed genes", sep = "")
  }
  
  # Create normalized counts plot for top n differentially expressed genes
  set.seed(1)
  plot = res_plot %>%
    ggplot(aes(x = SYMBOL , y = count, colour = cancer_type)) +
    geom_jitter(width = 0.3, 
                size = 2, 
                shape = 16,
                alpha = 0.5) +
    stat_summary(fun.data = mean_sdl,
                 position = position_dodge(width = 0.2),
                 linewidth = 1.1) +
    # The plotCounts() function displays the results in logarithmic scale,  
    # but for the ggplot we need to add a logarithmic y-axis scale.
    scale_y_log10() +
    scale_color_manual(values = c("#5d6d7e" , "#F08080" ),
                       labels = c(paste(level.1, " (control)", sep = ""), 
                                  level.2)) +
    labs(y = "normalized counts",
         title = plot_title) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          plot.title = element_text(size = rel(1)), 
          axis.title = element_text( size = rel(0.9)),
          legend.title = element_blank() )
  
  # return:
  list(res, kw, plot) 
  
}





