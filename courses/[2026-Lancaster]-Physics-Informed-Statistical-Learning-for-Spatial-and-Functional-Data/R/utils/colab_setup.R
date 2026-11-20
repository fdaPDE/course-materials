# Configure the environment only when running on Google Colab.
local({
  if (any(nzchar(Sys.getenv(c(
    "COLAB_RELEASE_TAG", "COLAB_BACKEND_VERSION"
  ))))) {
    repository <- "/content/course-materials"

    # Download the data and utilities.
    if (!dir.exists(repository)) {
      stopifnot(system2(
        "git",
        c(
          "clone", "--depth", "1",
          "https://github.com/fdaPDE/course-materials.git",
          shQuote(repository)
        ),
        stdout = FALSE
      ) == 0L)
    }

    setwd(file.path(
      repository,
      "courses",
      "[2026-Lancaster]-Physics-Informed-Statistical-Learning-for-Spatial-and-Functional-Data",
      "R", "notebook"
    ))

    # Install system libraries and restore the precompiled R packages.
    dependencies <- "../colab_binaries"
    package_library <- .libPaths()[1]
    extraction_dir <- tempfile()
    dir.create(extraction_dir)

    stopifnot(system2(
      "bash",
      shQuote(file.path(dependencies, "system-deps.sh")),
      stdout = FALSE
    ) == 0L)

    stopifnot(untar(
      file.path(dependencies, "colab_binaries.tar.gz"),
      exdir = extraction_dir,
      tar = "/usr/bin/tar"
    ) == 0L)

    archives <- list.files(
      extraction_dir,
      pattern = "\\.tar\\.gz$",
      recursive = TRUE,
      full.names = TRUE
    )

    for (archive in archives) {
      stopifnot(untar(
        archive,
        exdir = package_library,
        tar = "/usr/bin/tar"
      ) == 0L)
    }

    unlink(extraction_dir, recursive = TRUE)
  }
  invisible(NULL)
})
