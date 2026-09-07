pdac_paths <- function() {
  repo_root <- Sys.getenv("PDAC_CODE_REPO_ROOT", unset = getwd())
  project_root <- Sys.getenv("PDAC_PROJECT_ROOT", unset = repo_root)
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = FALSE)
  project_root <- normalizePath(project_root, winslash = "/", mustWork = FALSE)
  list(
    repo_root = repo_root,
    project_root = project_root,
    raw_data_root = file.path(project_root, "data"),
    external_data_root = file.path(project_root, "data_external"),
    analysis_root = file.path(project_root, "analysis"),
    results_root = file.path(project_root, "results"),
    figure_root = file.path(project_root, "analysis", "15_priority_figures")
  )
}

pdac_script_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)
    return(normalizePath(file.path(dirname(script_path), "../.."), winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
