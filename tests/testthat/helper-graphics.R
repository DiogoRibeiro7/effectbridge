# Route graphics output to a temp directory to avoid creating tracked artifacts
# such as tests/testthat/Rplots.pdf during automated test runs.
with_plot_sandbox <- function(code) {
  old_wd <- getwd()
  sandbox_dir <- tempfile("plot-sandbox-")
  dir.create(sandbox_dir)
  setwd(sandbox_dir)

  on.exit({
    setwd(old_wd)
    unlink(sandbox_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  # Ensure a graphics device is explicitly open in the sandbox.
  plot_path <- tempfile("plot-device-", tmpdir = sandbox_dir, fileext = ".pdf")
  grDevices::pdf(plot_path)
  dev_id <- grDevices::dev.cur()

  on.exit({
    open_devs <- grDevices::dev.list()
    if (!is.null(open_devs) && dev_id %in% open_devs) {
      grDevices::dev.off(which = dev_id)
    }
  }, add = TRUE)

  force(code)
}
