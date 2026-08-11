.spamtp_vignette_root <- function() {
  candidates <- c(".", "..")
  hit <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))]
  if (!length(hit)) {
    stop("Cannot locate the SpaMTP package root for vignette downloads.")
  }
  normalizePath(hit[[1]], mustWork = TRUE)
}

resolve_spamtp_vignette_file <- function(article, filename, url,
                                          expected_md5, attempts = 5L) {
  expected_md5 <- sub("^md5:", "", tolower(expected_md5))
  cache_root <- Sys.getenv(
    "SPAMTP_VIGNETTE_DATA_CACHE",
    unset = file.path(
      .spamtp_vignette_root(), "vignettes", ".knitr-cache", "downloads"
    )
  )
  cached_file <- file.path(cache_root, article, filename)
  partial_file <- paste0(cached_file, ".part")

  is_valid <- function(path) {
    isTRUE(file.exists(path)) && isTRUE(file.info(path)$size > 0) &&
      identical(unname(tools::md5sum(path)), expected_md5)
  }
  if (is_valid(cached_file)) return(cached_file)

  dir.create(dirname(cached_file), recursive = TRUE, showWarnings = FALSE)
  old_timeout <- getOption("timeout", 60)
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(1800, as.numeric(old_timeout)))
  last_error <- NULL

  for (attempt in seq_len(as.integer(attempts))) {
    unlink(partial_file)
    last_error <- tryCatch(
      {
        status <- suppressWarnings(utils::download.file(
          url, partial_file, mode = "wb", method = "libcurl", quiet = TRUE
        ))
        if (!identical(status, 0L)) {
          stop("download.file returned status ", status)
        }
        if (!is_valid(partial_file)) {
          stop("downloaded file failed its MD5 or size check")
        }
        NULL
      },
      error = identity
    )

    if (is.null(last_error)) {
      if (file.exists(cached_file)) unlink(cached_file)
      if (!file.rename(partial_file, cached_file)) {
        stop("Could not finalize downloaded vignette file: ", cached_file)
      }
      return(cached_file)
    }
    if (attempt < attempts) Sys.sleep(min(60, 15 * attempt))
  }

  unlink(partial_file)
  stop(
    "Failed to download and validate ", filename, " after ", attempts,
    " attempts: ", conditionMessage(last_error)
  )
}
