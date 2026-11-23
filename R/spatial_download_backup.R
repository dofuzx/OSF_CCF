# spatial_download_backup.R — Visium spatial data two-channel download (SOP section 4, O3)
# Channel 1 (primary): programmatic download via GEOquery
# Channel 2 (backup):  writes a manual-download checklist template for browser download
# Both channels converge on the same verification step before analysis.
suppressMessages(library(GEOquery))

source("config.R")  # must define:
#   spatial_gses <- c("GSExxxxxx", ...)   # Visium accessions selected at execution time
#   spatial_dir  <- "02_spatial"
#   log_csv      <- "00_meta/download_log.csv"

## --- helpers ---------------------------------------------------------------
geo_urls <- function(gse) {
  stub <- paste0(substr(gse, 1, 5), "nnn")
  list(
    supp = sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/", stub, gse),
    soft = sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/soft/%s_family.soft.gz",
                   stub, gse, gse))
}

## Files a Space Ranger output set must contain for Load10X_Spatial()
visium_required_files <- function(gse) {
  c("filtered_feature_bc_matrix.h5",            # or filtered_feature_bc_matrix/ folder
    "spatial/tissue_positions_list.csv",        # Space Ranger v1; v2 uses tissue_positions.csv
    "spatial/tissue_positions.csv",             # included in checklist; at least one of the two
    "spatial/scalefactors_json.json",
    "spatial/tissue_lowres_image.png",
    paste0(gse, "_family.soft.gz"))             # metadata (sample mapping, patient IDs)
}

write_manual_template <- function(gse, urls, dest) {
  req <- visium_required_files(gse)
  tpl <- data.frame(
    gse = gse,
    filename = req,
    url = rep("", length(req)),           # fill per-sample URL after inspecting the FTP index
    expected_size_mb = "", md5 = "",
    barcode_count = "", geo_statement_count = "", dimension_match = "",
    downloaded_at = "", downloaded_by = "", method = "manual",
    stringsAsFactors = FALSE)
  tpl$url[tpl$filename == paste0(gse, "_family.soft.gz")] <- urls$soft
  tpl$url[tpl$url == ""] <- urls$supp                       # default: supplement folder index
  write.csv(tpl, file.path(dest, paste0(gse, "_manual_checklist.csv")), row.names = FALSE)
  message(sprintf(
    paste0("[", gse, "] Programmatic download failed.\n",
           "MANUAL BACKUP: open %s in a browser, download the files listed in\n%s,\n",
           "then fill md5/barcode_count/dimension_match and re-run this script.\n"),
    urls$supp, file.path(dest, paste0(gse, "_manual_checklist.csv"))))
}

verify_visium_files <- function(gse, dest) {
  req <- visium_required_files(gse)
  pos_file <- file.path(dest, "spatial",
                        if (file.exists(file.path(dest, "spatial/tissue_positions_list.csv")))
                          "tissue_positions_list.csv" else "tissue_positions.csv")
  h5_file <- file.path(dest, "filtered_feature_bc_matrix.h5")
  n_bc <- NA_integer_
  if (file.exists(h5_file)) {
    n_bc <- tryCatch({
      suppressMessages(requireNamespace("rhdf5", quietly = TRUE))
      length(names(rhdf5::h5read(h5_file, "matrix/barcodes")))
    }, error = function(e) NA_integer_)
  }
  checks <- list(
    all_required_present = all(file.exists(file.path(dest, req)) |
                                 file.exists(pos_file) | req == basename(pos_file)),
    nonempty = all(file.info(file.path(dest, req[file.exists(file.path(dest, req))]))$size > 0),
    barcode_count = n_bc)
  checks
}

append_log <- function(gse, dest, method, checks) {
  row <- data.frame(time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    accession = gse, target = dest, method = method,
                    all_required_present = checks$all_required_present,
                    nonempty = checks$nonempty,
                    barcode_count = checks$barcode_count,
                    dimension_match = NA,   # compare barcode_count vs GEO statement; fill before analysis
                    stringsAsFactors = FALSE)
  write.table(row, log_csv, append = file.exists(log_csv),
              sep = ",", row.names = FALSE, col.names = !file.exists(log_csv))
}

## --- main loop -------------------------------------------------------------
for (gse in spatial_gses) {
  dest <- file.path(spatial_dir, gse)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  urls <- geo_urls(gse)

  method <- "programmatic"
  ok <- tryCatch({
    suppressMessages(getGEOSuppFiles(gse, baseDir = dest, makeDirectory = FALSE,
                                     fetch_files = TRUE, timeout = 300))
    all(file.exists(file.path(dest, visium_required_files(gse))) |
          file.exists(file.path(dest, "spatial/tissue_positions.csv")))
  }, error = function(e) FALSE)
  if (!ok) { method <- "manual"; write_manual_template(gse, urls, dest) }

  checks <- verify_visium_files(gse, dest)
  append_log(gse, dest, method, checks)
}

## After any manual download: re-run this script to verify + log, then load:
# seu <- Load10X_Spatial(data.dir = file.path(spatial_dir, gse), filename = "filtered_feature_bc_matrix.h5")
