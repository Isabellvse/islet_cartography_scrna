# Description -------------------------------------------------------------
# Create plots for quality control metrics for each study

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load --------------------------------------------------------------------
## Quality metrics ----
paths <- base::list.files(path = here::here("islet_cartography_scrna/data/quality_control/first_pass/quality_metrics"),
                          pattern = "quality_metrics.csv",
                          full.names = TRUE, 
                          recursive = TRUE)

base::names(paths) <- purrr::map_chr(paths, ~ stringr::str_extract(.x, pattern = "(?<=quality_metrics/)[^/]+(?=_quality_metrics\\.csv)"))

quality_met <- purrr::map(paths, ~vroom::vroom(.x, delim = ",", col_names = TRUE))

## Thresholds ----
thresholds <- vroom::vroom(here::here("islet_cartography_scrna/data/quality_control/first_pass/thresholds.csv"), 
                           delim = ",", 
                           col_names = TRUE) 

# Preprocess --------------------------------------------------------------
## Add threshold to quality matrix
# quality_met_thres <- purrr::map(quality_met, function(df){
#   if(base::unique(df$platform) %in% c("droplet", "plate_barcode")){
#     thres <- thresholds |> dplyr::select(name, dplyr::all_of(qc_thres_droplet))
#     df_thres <- dplyr::left_join(x = df, y = thres, by = "name")
#     
#   } else if (base::unique(df$platform) %in% c("plate")) {
#     thres <- thresholds |> dplyr::select(name, dplyr::all_of(qc_thres_plate))
#     df_thres <- dplyr::left_join(x = df, y = thres, by = "name")
#     
#   } else {
#     message("platform could not be found")
#   }
#   return(df_thres)
# })

## Add unmaped reads % to dataframe
## Star quality
quality_met <- purrr::map(quality_met, function(df) {
  df <- df |> 
    dplyr::rowwise() %>%
    dplyr::mutate("Unmapped_reads_%" = sum(dplyr::c_across(tidyselect::starts_with("%_of_reads_unmapped_")), na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  return(df)
})

qs2::qs_save(quality_met, here::here("islet_cartography_scrna/data/quality_control/first_pass/quality_met.qs2"))
quality_met <- qs2::qs_read(here::here("islet_cartography_scrna/data/quality_control/first_pass/quality_met.qs2"))

# Quality metrics ---------------------------------------------------------
purrr::imap(quality_met, function(df, name) {
  pdf(
    file = base::paste0(
      here::here(
        "islet_cartography_scrna/data/quality_control/first_pass/plots/"
      ),
      name,
      "_threshold.pdf"
    ),
    height = 2,
    width = 2
  )
  
  thresholds_filtered <- dplyr::filter(thresholds, name == !!name)
  
  if (base::unique(df$platform) %in% c("droplet", "plate_barcode")) {
    df |>
      dplyr::select(dplyr::all_of(qc_met_thres_droplet)) |>
      purrr::imap(~ {
        lower_col <- paste0("threshold_", .y, "_lower")
        upper_col <- paste0("threshold_", .y, "_upper")
        lower <- thresholds_filtered[[lower_col]][1]
        upper <- thresholds_filtered[[upper_col]][1]
        
        # Check if lower or upper is NA or NULL and replace accordingly
        if (is.na(lower) || is.null(lower)) {
          lower <- NULL
        }
        
        if (is.na(upper) || is.null(upper)) {
          upper <- NULL
        }
        
        print(plot_hist_qc_thres(.x, .y, lower, upper))
      })
  } else if (base::unique(df$platform) %in% c("plate")) {
    df |>
      dplyr::select(dplyr::all_of(qc_met_thres_plate)) |>
      purrr::imap(~ {
        lower_col <- paste0("threshold_", .y, "_lower")
        upper_col <- paste0("threshold_", .y, "_upper")
        lower <- thresholds_filtered[[lower_col]][1]
        upper <- thresholds_filtered[[upper_col]][1]
        
        # Check if lower or upper is NA or NULL and replace accordingly
        if (is.na(lower) || is.null(lower)) {
          lower <- NULL
        }
        
        if (is.na(upper) || is.null(upper)) {
          upper <- NULL
        }
        
        print(plot_hist_qc_thres(.x, .y, lower, upper))
      })
  }
  
  dev.off()
})

pdf(
  file = here::here("islet_cartography_scrna/data/quality_control/first_pass/plots/threshold.pdf"),
  height = 2,
  width = 2
)
purrr::imap(quality_met, function(df, name) {

  thresholds_filtered <- dplyr::filter(thresholds, name == !!name)
  subtitle <- name
  
  if (base::unique(df$platform) %in% c("droplet", "plate_barcode")) {
    df |>
      dplyr::select(dplyr::all_of(qc_met_thres_droplet)) |>
      purrr::imap(~ {
        lower_col <- paste0("threshold_", .y, "_lower")
        upper_col <- paste0("threshold_", .y, "_upper")
        lower <- thresholds_filtered[[lower_col]][1]
        upper <- thresholds_filtered[[upper_col]][1]
        
        # Check if lower or upper is NA or NULL and replace accordingly
        if (is.na(lower) || is.null(lower)) {
          lower <- NULL
        }
        
        if (is.na(upper) || is.null(upper)) {
          upper <- NULL
        }
        
        print(plot_hist_qc_thres(.x, .y, lower, upper, subtitle = subtitle))
      })
  } else if (base::unique(df$platform) %in% c("plate")) {
    df |>
      dplyr::select(dplyr::all_of(qc_met_thres_plate)) |>
      purrr::imap(~ {
        lower_col <- paste0("threshold_", .y, "_lower")
        upper_col <- paste0("threshold_", .y, "_upper")
        lower <- thresholds_filtered[[lower_col]][1]
        upper <- thresholds_filtered[[upper_col]][1]
        
        # Check if lower or upper is NA or NULL and replace accordingly
        if (is.na(lower) || is.null(lower)) {
          lower <- NULL
        }
        
        if (is.na(upper) || is.null(upper)) {
          upper <- NULL
        }
        
        print(plot_hist_qc_thres(.x, .y, lower, upper, subtitle = subtitle))
      })
  }

})
dev.off()


# QC thresholds per donor -------------------------------------------------
qc_donor <- quality_met |> 
  purrr::modify_depth(1, ~ base::split(.x, .x$ic_id_donor)) |> 
  purrr::list_flatten()

thresholds_donor <- qc_donor |> 
  purrr::modify_depth(1, ~ dplyr::select(.x, ic_id_donor, name) |> 
                        dplyr::distinct()) |> 
  dplyr::bind_rows(.id = "ic_id_donor") |> 
  dplyr::left_join(y = thresholds, by = "name")


pdf(
  file = here::here(
      "islet_cartography_scrna/data/quality_control/first_pass/plots/threshold_per_donor.pdf"),
  height = 2,
  width = 2
)
purrr::imap(qc_donor, function(df, name) {
    thresholds_filtered <- dplyr::filter(thresholds_donor, ic_id_donor == !!name)
    subtitle <- name
    if (base::unique(df$platform) %in% c("droplet", "plate_barcode")) {
      df |>
        dplyr::select(dplyr::all_of(qc_met_thres_droplet)) |>
        purrr::imap(~ {
          lower_col <- paste0("threshold_", .y, "_lower")
          upper_col <- paste0("threshold_", .y, "_upper")
          lower <- thresholds_filtered[[lower_col]][1]
          upper <- thresholds_filtered[[upper_col]][1]
          
          # Check if lower or upper is NA or NULL and replace accordingly
          if (is.na(lower) || is.null(lower)) {
            lower <- NULL
          }
          
          if (is.na(upper) || is.null(upper)) {
            upper <- NULL
          }
          
          print(plot_hist_qc_thres(.x, .y, lower, upper, subtitle = subtitle))
        })
    } else if (base::unique(df$platform) %in% c("plate")) {
      df |>
        dplyr::select(dplyr::all_of(qc_met_thres_plate)) |>
        purrr::imap(~ {
          lower_col <- paste0("threshold_", .y, "_lower")
          upper_col <- paste0("threshold_", .y, "_upper")
          lower <- thresholds_filtered[[lower_col]][1]
          upper <- thresholds_filtered[[upper_col]][1]
          
          # Check if lower or upper is NA or NULL and replace accordingly
          if (is.na(lower) || is.null(lower)) {
            lower <- NULL
          }
          
          if (is.na(upper) || is.null(upper)) {
            upper <- NULL
          }
          
          print(plot_hist_qc_thres(.x, .y, lower, upper, subtitle = subtitle))
        })
    }
  })
dev.off()



