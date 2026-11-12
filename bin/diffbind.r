#!/usr/bin/env Rscript
#
# DiffBind for a single contrast: two groups of one antibody.
#
# Usage:
#   diffbind.R --samplesheet=FILE --prefix=NAME --numerator=GROUP --denominator=GROUP
#              [--peak_type=narrow|broad] [--cores=N] [--fdr=0.05] [--min_overlap=2]
#
# The sample sheet (written by the DIFFBIND module) has columns
#   SampleID, Tissue, Condition, Treatment, Replicate, [ControlID, bamControl,] bamReads, Peaks
# GROUP is Tissue_Condition_Treatment, e.g. PFC_WT_DMSO. log2FC = numerator / denominator.
#

suppressPackageStartupMessages(library(DiffBind))

## ------------------------------------------------------------------ arguments
parse_args <- function(x) {
  x  <- sub("^--", "", x)
  kv <- regmatches(x, regexpr("=", x), invert = TRUE)      # split on the first "=" only
  setNames(vapply(kv, function(p) if (length(p) > 1) p[2] else "TRUE", ""),
           vapply(kv, `[`, "", 1))
}
opt     <- parse_args(commandArgs(trailingOnly = TRUE))
get_opt <- function(name, default = NULL) if (name %in% names(opt)) opt[[name]] else default

sheet_file  <- get_opt("samplesheet")
prefix      <- get_opt("prefix", "diffbind")
numerator   <- get_opt("numerator")
denominator <- get_opt("denominator")
peak_type   <- get_opt("peak_type", "narrow")
cores       <- as.integer(get_opt("cores", "1"))
fdr         <- as.numeric(get_opt("fdr", "0.05"))
min_overlap <- as.numeric(get_opt("min_overlap", "2"))

for (a in c("sheet_file", "numerator", "denominator")) {
  if (is.null(get(a))) stop("Missing required argument for ", a)
}
if (!peak_type %in% c("narrow", "broad")) stop("--peak_type must be 'narrow' or 'broad'")

## ------------------------------------------------------------------ sample sheet
samples <- read.csv(sheet_file, colClasses = "character", na.strings = "")

has_ctrl <- "bamControl" %in% colnames(samples)
message(if (has_ctrl) "Controls present: subtracting control reads" else "No controls: using raw IP counts")

samples$PeakCaller <- if (peak_type == "broad") "bed" else "narrow"

factor_cols    <- intersect(c("Tissue", "Condition", "Treatment"), colnames(samples)) # Group name from the factor columns present (blank factors are left out of the sheet),
samples$Factor <- do.call(paste, c(unname(as.list(samples[factor_cols])), sep = "_")) # built the same way as in Nextflow, e.g. PFC_WT_DMSO or PFC_WT

# Replicate numbers: keep if complete, otherwise number samples within each group
rep_num <- suppressWarnings(as.integer(samples$Replicate))
if (any(is.na(rep_num))) rep_num <- ave(seq_len(nrow(samples)), samples$Factor, FUN = seq_along)
samples$Replicate <- as.integer(rep_num)

# The sheet must contain exactly the two groups being compared, each replicated
n_per_group <- table(samples$Factor)
print(n_per_group)
if (!setequal(names(n_per_group), c(numerator, denominator))) {
  stop("Sample sheet groups (", paste(names(n_per_group), collapse = ", "),
       ") do not match numerator/denominator (", numerator, ", ", denominator, ")")
}
if (any(n_per_group < 2)) stop("Each group needs >= 2 replicates")

staged  <- unique(c(samples$bamReads, samples$bamControl, samples$Peaks))   # bamControl is NULL without controls
missing <- staged[!file.exists(staged)]
if (length(missing)) stop("Files in the sample sheet were not staged: ", paste(missing, collapse = ", "))

# Save the completed sheet, so the published copy is exactly what DiffBind used
write.csv(samples, paste0(prefix, ".samplesheet.csv"), row.names = FALSE)

## ------------------------------------------------------------------ DiffBind
dba_obj <- dba(sampleSheet = samples, minOverlap = min_overlap)
dba_obj$config$cores       <- cores
dba_obj$config$RunParallel <- cores > 1
dba_obj$config$th          <- fdr
# Peaks are blacklist-filtered upstream, and DiffBind can't recognise the genome from
# RefSeq-style chromosome names, so skip its automatic blacklist/greylist.
dba_obj$config$doBlacklist <- FALSE
dba_obj$config$doGreylist  <- FALSE

count_args <- list(dba_obj, minOverlap = min_overlap, bParallel = cores > 1, bSubControl = has_ctrl)
if (peak_type == "broad") count_args$summits <- FALSE        # don't recenter broad peaks on summits
dba_obj <- do.call(dba.count, count_args)
dba_obj <- dba.normalize(dba_obj)

message("Contrast: ", numerator, " vs ", denominator)
dba_obj <- dba.contrast(dba_obj, design = "~ Factor", contrast = c("Factor", numerator, denominator))
dba_obj <- dba.analyze(dba_obj, bBlacklist = FALSE, bGreylist = FALSE)
print(dba.show(dba_obj, bContrasts = TRUE))
saveRDS(dba_obj, paste0(prefix, ".dba.rds"))

## ------------------------------------------------------------------ outputs
res <- as.data.frame(dba.report(dba_obj, contrast = 1, th = 1))     # th = 1: every consensus site
write.csv(res, paste0(prefix, ".results.csv"), row.names = FALSE)
message(sprintf("%d of %d sites at FDR < %g", sum(res$FDR < fdr, na.rm = TRUE), nrow(res), fdr))

scores <- dba.peakset(dba_obj, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
write.csv(scores, paste0(prefix, ".consensus_scores.csv"), row.names = FALSE)

safe_plot <- function(expr) tryCatch(expr, error = function(e) message("Plot skipped: ", conditionMessage(e)))
pdf(paste0(prefix, ".plots.pdf"), width = 8, height = 8)
safe_plot(dba.plotHeatmap(dba_obj))                 # sample correlation
safe_plot(dba.plotPCA(dba_obj, label = DBA_ID))
safe_plot(dba.plotMA(dba_obj, contrast = 1))
safe_plot(dba.plotVolcano(dba_obj, contrast = 1))
invisible(dev.off())
