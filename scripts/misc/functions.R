#' Create Directories if They Do Not Exist
#'
#' This function checks if each directory in the provided list exists. If any directory
#' does not exist, it creates the directory (along with any necessary parent directories).
#' It prints a message for each directory indicating whether it was created or already exists.
#'
#' @param output_dir A character vector or list of directory paths to be checked and created.
#'
#' @return NULL This function performs a side-effect (creating directories) and does not return anything.
#'
#' @export
#'
#' @examples
#' # Define a list of directories
#' output_dirs <- list("data/quality_control/rna", "data/plots", "data/other_output")
#'
#' # Call the function to create the directories
#' create_directories(output_dirs)
create_directories <- function(output_dir) {
  purrr::walk(output_dir, ~{
    if (!dir.exists(.x)) {
      dir.create(.x, recursive = TRUE)  # create the directory if it doesn't exist
      print(paste0(.x, " has been created!"))
    } else {
      print(paste0(.x, " already exists!"))
    }
  })
}

load_data_with_classes <- function(data_file = "data.csv", class_file = "column_classes.csv") {
  # Read the column classes
  df_classes <- vroom::vroom(class_file, delim = ",", col_types = cols(column = vroom::col_character(),
                                                                       class = vroom::col_character()))

  # Create a named list of column types based on the class information
  col_types <- purrr::set_names(
    purrr::map(df_classes$class, function(x) {
      base::switch(x,
             "logical" = vroom::col_logical(),
             "integer" = vroom::col_integer(),
             "numeric" = vroom::col_double(),
             "character" = vroom::col_character(),
             "date" = vroom::col_date(format = ""),
             "time" = vroom::col_time(format = ""),
             "datetime" = vroom::col_datetime(format = ""),
             "factor" = vroom::col_factor(levels = NULL, ordered = FALSE),
             vroom::col_guess()  # Default to guessing for unknown types
      )
    }),
    df_classes$column
  )

  # Use vroom to load the data with the correct column types
  df_loaded <- vroom::vroom(data_file, col_types = col_types, delim = ",", n_max = 10000)

  return(df_loaded)
}

#' Load and Process Star Quality Data
#'
#' This function loads a data file, processes any columns containing percentage 
#' signs (`%`) by removing the `%` and converting the values to numeric. It then 
#' adds a `name` column, which is assigned the value of the `i` parameter, and 
#' relocates this column to the first position in the resulting data frame.
#'
#' @param path A character string specifying the path to the .csv file.
#' @param i A value (numeric, character, etc.) that will be used to populate the 
#'          `name` column in the data frame. This could be an identifier or a 
#'          label for the data set.
#'
#' @return A data frame where:
#'   - Columns containing percentage signs (`%`) are converted to numeric by 
#'     removing the `%` and converting the remaining values to numeric.
#'   - A new column named `name` is added, containing the value passed in the `i`
#'     parameter.
#'   - The `name` column is relocated to the first position.
#'
#' @export
#'
#' @examples
#' # Example usage:
#' # df <- load_star_quality("path/to/your/file.csv", "SampleName")
#' # head(df)
#'
load_star_quality <- function(path, i){
  df <- path |> vroom::vroom()  |> 
    dplyr::mutate(dplyr::across(.cols = dplyr::contains("%"),
                                .fns = ~stringr::str_remove(.x, pattern = "%") |> 
                                  base::as.numeric()),
                  name = i) |> 
    dplyr::relocate(name)
  return(df)
  }

