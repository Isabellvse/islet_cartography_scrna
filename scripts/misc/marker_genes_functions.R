meta_analysis <- function(df){
  df_filt <- tidyr::drop_na(df, log2FoldChange)
  
  res <- meta::metagen(TE = log2FoldChange, seTE = lfcSE, studlab = dataset,
                       data = df_filt, n.e = n_donors, method.tau = "DL", prediction = TRUE,
                       common = FALSE, random = TRUE)
  
  res2 <- metabind(list(res))$data 
  return(res2)
}

process_and_save <- function(path, name, save_path) {
  result <- vroom::vroom(path) |> 
    dplyr::select(gene_symbol, dataset, target, log2FoldChange, lfcSE, n_donors) |> 
    tidyr::nest(.by = 'gene_symbol') |> 
    dplyr::filter(purrr::map_int(data, ~sum(!is.na(.x$log2FoldChange))) > 1) |> 
    dplyr::mutate(meta = furrr::future_map(data, meta_analysis)) |> 
    dplyr::select(-data) |> 
    tidyr::unnest('meta')
  
  vroom::vroom_write(result, paste0(save_path, "/", name, ".csv"))
}

