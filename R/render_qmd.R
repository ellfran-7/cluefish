#' Render and Move Quarto Documents
#'
#' This function renders a Quarto `.qmd` document to a specified output format and then moves
#' the resulting file to a designated output path. It uses `quarto::quarto_render` to
#' perform the rendering and `fs::file_move` to handle the file movement.
#'
#' @param input_file A character string specifying the path to the input Quarto `.qmd` file.
#' @param output_path A character string specifying the directory where the rendered file should be moved.
#' @param ... Additional arguments passed to `quarto::quarto_render`.
#'
#' @return This function does not return a value. It performs rendering and file moving operations, and displays a message indicating the file movement.
#' @export
#'
#' @examples

render_qmd <- function(input_file, output_file, output_path, file_ext, ...) {
  
  # render the input document and output file will be in the
  # current working directory.
  quarto::quarto_render(input = input_file, output_format = file_ext, ...)
  
  # name of the rendered output file
  output_name <- paste0(output_file, ".", file_ext)
  
  # Ensure the output directory exists
  if (!fs::dir_exists(output_path)) {
    fs::dir_create(output_path)
  }
  
  # move the file to the output path
  fs::file_move(paste0(output_name), output_path)
  
  message(paste0(paste0(output_name, collapse = " and "), " moved to ", output_path))
  
  message(paste0(output_name, collapse = " and ", "\n"))
  message(paste0("--> moved to -->"))
  message(output_path)
}
