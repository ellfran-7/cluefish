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
  
  # Define the name of the rendered output file
  output_name <- paste0(output_file, ".", file_ext)
  
  # Render the input document with the specified output format
  # The rendered file will initially be saved in the current working directory
  quarto::quarto_render(input = input_file, output_file = output_name, output_format = file_ext, ...)
  
  # Check if the output directory exists; create it if it doesn't
  if (!fs::dir_exists(output_path)) {
    fs::dir_create(output_path)
    message(paste("Directory", output_path, "did not exist. Creating it now."))
  }
  
  # Move the rendered output file to the specified output directory
  fs::file_move(output_name, file.path(output_path, output_name))
  
  # Notify that the file has been successfully moved
  message(paste("File", output_name, "has been successfully moved to", output_path))
}
