# Restore precompiled R packages and their system dependencies.
local({
  base_url <- paste0(
    "https://raw.githubusercontent.com/fdaPDE/course-materials/main/",
    "courses/%5B2026-Lancaster%5D-",
    "Physics-Informed-Statistical-Learning-for-Spatial-and-Functional-Data/",
    "R/colab_binaries/"
  )

  options(timeout = max(600, getOption("timeout")))

  setup_dir <- tempfile("dependency-setup-", tmpdir = "/content")
  dir.create(setup_dir)

  log_file <- "/content/dependency-setup.log"
  writeLines(character(), log_file)

  run_command <- function(command, args, error_message) {
    output <- suppressWarnings(system2(
      command, args,
      stdout = TRUE,
      stderr = TRUE
    ))

    cat(output, sep = "\n", file = log_file, append = TRUE)
    status <- attr(output, "status")

    if (!is.null(status) && status != 0L) {
      stop(
        error_message, "\n",
        paste(tail(output, 30L), collapse = "\n"),
        "\nFull log: ", log_file,
        call. = FALSE
      )
    }

    invisible(NULL)
  }

  # Download the package bundle and installation metadata.
  filenames <- c(
    "colab_binaries.tar.gz",
    "dependencies.csv",
    "required-packages.txt",
    "system-deps.sh"
  )

  for (filename in filenames) {
    download.file(
      paste0(base_url, filename),
      file.path(setup_dir, filename),
      mode = "wb",
      method = "libcurl",
      quiet = TRUE
    )
  }

  manifest <- read.csv(
    file.path(setup_dir, "dependencies.csv"),
    stringsAsFactors = FALSE
  )

  if (!all(c("Package", "File", "MD5") %in% names(manifest))) {
    stop(
      "The manifest must contain Package, File and MD5 columns.",
      call. = FALSE
    )
  }

  # Extract the bundle and locate the individual package archives.
  extraction_dir <- file.path(setup_dir, "extracted")
  dir.create(extraction_dir)

  run_command(
    "tar",
    c(
      "-xzf", shQuote(file.path(setup_dir, "colab_binaries.tar.gz")),
      "-C", shQuote(extraction_dir)
    ),
    "Failed to extract the dependency bundle."
  )

  archives <- list.files(
    extraction_dir,
    pattern = "\\.(tar\\.gz|tgz)$",
    recursive = TRUE,
    full.names = TRUE
  )

  resolved_archives <- vapply(manifest$File, function(filename) {
    matches <- archives[basename(archives) == filename]

    if (length(matches) != 1L) {
      stop(
        "Expected exactly one archive matching: ", filename,
        call. = FALSE
      )
    }

    matches
  }, character(1))

  # Verify archive integrity before installation.
  checksums <- unname(tools::md5sum(resolved_archives))
  invalid <- is.na(checksums) | checksums != manifest$MD5

  if (any(invalid)) {
    stop(
      "Checksum verification failed for: ",
      paste(manifest$File[invalid], collapse = ", "),
      call. = FALSE
    )
  }

  # Install external system libraries.
  run_command(
    "bash",
    shQuote(file.path(setup_dir, "system-deps.sh")),
    "Failed to install the required system libraries."
  )

  # Restore the precompiled packages into a dedicated library.
  package_library <- "/content/r-library"
  dir.create(package_library, showWarnings = FALSE)

  for (i in seq_along(resolved_archives)) {
    run_command(
      "tar",
      c(
        "-xzf", shQuote(resolved_archives[i]),
        "-C", shQuote(package_library)
      ),
      paste0("Failed to extract package: ", manifest$Package[i])
    )
  }

  .libPaths(c(package_library, .libPaths()))
  Sys.setenv(R_LIBS_USER = package_library)

  # Verify that all required packages can be loaded.
  required_packages <- trimws(readLines(
    file.path(setup_dir, "required-packages.txt"),
    warn = FALSE
  ))
  required_packages <- required_packages[nzchar(required_packages)]

  load_errors <- character()

  for (package in required_packages) {
    tryCatch(
      suppressPackageStartupMessages(loadNamespace(package)),
      error = function(e) {
        load_errors <<- c(
          load_errors,
          paste0(package, ": ", conditionMessage(e))
        )
      }
    )
  }

  if (length(load_errors)) {
    stop(
      "Some dependencies could not be loaded:\n",
      paste(load_errors, collapse = "\n"),
      call. = FALSE
    )
  }

  invisible(NULL)
})
