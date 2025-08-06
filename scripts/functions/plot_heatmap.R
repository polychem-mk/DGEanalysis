# plot_heatmap() is a custom function that filters cancer type and creates 
# heatmaps; it adds column annotations ("mutational_subclass") and row
# annotations ("batch").

plot_heatmap <- function(cm,                  # correlation matrix
                         md,                  # metadata 
                         level_to_compare,    # one of cancer_type levels
                         method = "complete", # clustering_method for pheatmap()
    # values for annotation_col and annotation_row arguments for pheatmap():
                         ann_cols = "mutational_subclass",
                         ann_rows = "batch"){
  
  # indexes for cancer_type i in metadata
  ind = which(md$cancer_type == level_to_compare)
  
  # dataframe with values for column annotations
  ann_cols_df = dplyr::select(md, all_of(ann_cols))
  
  # dataframe with values for row annotations
  ann_rows_df = dplyr::select(md, all_of(ann_rows))
  
  # plot name 
  plot_name = paste("Correlation heatmap,", level_to_compare, "samples")
  
  # unique values for column annotations
  ann_cols_lvls = as.character(unique(ann_cols_df[ind , 1]))
  
  # unique values for row annotations
  ann_rows_lvls = as.character(unique(ann_rows_df[ind , 1]))
  
  plot_color  = list(
    batch = c("1" = "#6c8db1",  "2" =  "#b16c8a" ,"3" =  "#3CB371", 
              "4" ="#87CEEB"  , "5" = "#FF8674" , "6" = "#ADFF2F"),
    
    mutational_subclass = c(HER2 = "#f7dc6f", wt = "#85929e",
                            Triple_Negative = "#6ca05d", HER2_PIK3CA = "#b7950b",
                            PIK3CA ="#db804f", KRAS = "#f941d7", MET = "#fbb7fb",
                            EGFR_MET ="#6c408c", EGFR =  "#c0f1ea",
                            KRAS_MET =   "#c7b9ff" ) ,
    
    cancer_type = c(breast = "#C11B17", hepatobiliary = "#5CB3FF",
                    CRC = "#a629f5", GBM = "#a7d170", lung = "#f78909", 
                    pancreas = "#0711fe", HC = "#474747"))
  
  # In 'ann_color', only keep values that are present in the metadata subset
  # with cancer_type i, otherwise all levels will be displayed in the heatmap
  # plot:
  ann_color = plot_color
  ann_color[[ann_cols]] = plot_color[[ann_cols]][ann_cols_lvls]
  ann_color[[ann_rows]] = plot_color[[ann_rows]][ann_rows_lvls]
  
  # Filter the correlation matrix and make the subset for cancer_type i,
  cm[ind, ind] %>%          # create a heat map:
    pheatmap( annotation_col =  ann_cols_df, 
              annotation_row =  ann_rows_df,
              annotation_colors = ann_color,
              show_rownames = FALSE,
              show_colnames = FALSE,
              clustering_method = method,
              main = plot_name)
}



