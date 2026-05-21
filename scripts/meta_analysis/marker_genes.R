# Description -------------------------------------------------------------
# Meta analysis of marker gene identification

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)


# Load --------------------------------------------------------------------
df <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/deseq_onevsother/deseq2_one_vs_all_per_study.csv"))


# INS ---------------------------------------------------------------------
df_ins <- df |> 
  dplyr::select()
  
res <- metagen(TE = logFC,
                seTE = lfcSE,
                studlab = study,
                data = ins_data,
                sm = "SMD",  # Summary measure
                fixed = FALSE,
                random = TRUE,
                method.tau = "REML",
                method.random.ci = "HK",
                title = "INS Gene Expression in Beta Cells")

