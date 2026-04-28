# ============================================================
# collect_opt_results.R
# ------------------------------------------------------------
# Durchsucht rekursiv alle .rds-Dateien, deren Name mit
# "opt_results" beginnt, zieht aus Dateiname UND Dateiinhalt
# die wichtigsten Infos (Algorithmus, Dataset, Gewichte,
# beste Parameter, bester Zielwert, Iterationen, Konvergenz)
# und schreibt eine einzige zusammengefasste Tabelle.
#
# Nutzung:
#   cd ~/code/mt_dsnow
#   Rscript calibration/collect_opt_results.R
# Oder die Konstante ROOT unten anpassen.
# ============================================================

suppressPackageStartupMessages({
  # Keine Zusatzpakete noetig - alles Base R.
})

# -----------------------------------------------------------
# Pfade
# -----------------------------------------------------------
ROOT     <- "~/code/mt_dsnow/calibration"                     # Startordner fuer die Suche
OUT_CSV  <- "~/code/mt_dsnow/calibration/opt_results_summary.csv"
OUT_RDS  <- "~/code/mt_dsnow/calibration/opt_results_summary.rds"

ROOT    <- normalizePath(ROOT, mustWork = TRUE)
message("Suche unter: ", ROOT)

# -----------------------------------------------------------
# 1) Alle relevanten Dateien finden
# -----------------------------------------------------------
files <- list.files(
  ROOT,
  pattern    = "^opt_results.*\\.rds$",
  recursive  = TRUE,
  full.names = TRUE
)

# Archiv-Verzeichnis ausschliessen (enthaelt veraltete Runs)
files <- files[!grepl("/archieve/", files, ignore.case = TRUE)]
files <- files[!grepl("/Archive/",  files, ignore.case = TRUE)]

message("Gefundene Dateien: ", length(files))
if (length(files) == 0) stop("Keine passenden .rds-Dateien gefunden.")

# -----------------------------------------------------------
# 2) Hilfsfunktionen
# -----------------------------------------------------------

# "0p3" -> 0.3, "0" -> 0, "1" -> 1
parse_weight <- function(s) {
  if (is.na(s) || s == "") return(NA_real_)
  as.numeric(sub("p", ".", s, fixed = TRUE))
}

# Dataset aus Pfad ableiten
infer_dataset <- function(path) {
  if (grepl("calibration_Win21",    path)) return("Win21")
  if (grepl("calibration_SNOWPACK", path)) return("SNOWPACK")
  NA_character_
}

# Algorithmus + Gewichte aus Dateinamen holen
parse_filename <- function(path) {
  fname <- basename(path)
  stem  <- sub("\\.rds$", "", fname)

  # Algorithmus: "_DE" irgendwo im Stamm -> DE, sonst Nelder-Mead
  algo <- if (grepl("(^|_)DE(_|$)", stem)) "DE" else "Nelder-Mead"

  # Dataset-Tag im Namen (z.B. "Win21") - wird weiter unten mit Pfad abgeglichen
  tag <- NA_character_
  m <- regmatches(stem, regexpr("Win21", stem))
  if (length(m) && nchar(m)) tag <- "Win21"

  # Gewichte aus "SWE_NRMSE_<w>__RHO_NRMSE_<w>__SWE_NBIAS_<w>"
  get_w <- function(key) {
    rx <- paste0(key, "_([0-9p]+)(?:__|$)")
    m  <- regmatches(stem, regexec(rx, stem))[[1]]
    if (length(m) >= 2) parse_weight(m[2]) else NA_real_
  }

  list(
    file         = fname,
    stem         = stem,
    algorithm    = algo,
    dataset_tag  = tag,
    w_SWE_NRMSE  = get_w("SWE_NRMSE"),
    w_RHO_NRMSE  = get_w("RHO_NRMSE"),
    w_SWE_NBIAS  = get_w("SWE_NBIAS")
  )
}

# Versucht, aus einem beliebigen Optimizer-Objekt die Essenz zu ziehen.
# Unterstuetzt: optim(), optimx(), DEoptim(), eigene Listen mit
# "par"/"bestmem"/"value"/"bestval"/"fn" etc.
extract_summary <- function(obj) {
  out <- list(
    best_value   = NA_real_,
    best_par     = NA_character_,
    n_par        = NA_integer_,
    iterations   = NA_integer_,
    convergence  = NA_character_,
    obj_class    = paste(class(obj), collapse = "/"),
    obj_names    = NA_character_
  )
  if (!is.null(names(obj))) out$obj_names <- paste(names(obj), collapse = "|")

  # Best-Vektor der Parameter zusammenbauen (als "name=value" String)
  format_par <- function(p) {
    if (is.null(p) || length(p) == 0) return(NA_character_)
    if (is.null(names(p))) names(p) <- paste0("p", seq_along(p))
    paste(sprintf("%s=%.6g", names(p), as.numeric(p)), collapse = "; ")
  }

  # ---- DEoptim ----
  if (inherits(obj, "DEoptim") || "optim" %in% names(obj) && is.list(obj$optim)) {
    # DEoptim: $optim$bestmem, $optim$bestval, $optim$iter
    opt <- obj$optim %||% list()
    out$best_value  <- as.numeric(opt$bestval %||% NA)
    out$best_par    <- format_par(opt$bestmem)
    out$n_par       <- length(opt$bestmem %||% NULL)
    out$iterations  <- as.integer(opt$iter %||% NA)
    out$convergence <- as.character(opt$nfeval %||% NA)  # nur zusaetzlich
    return(out)
  }

  # ---- optim() ----
  if (is.list(obj) && all(c("par", "value") %in% names(obj))) {
    out$best_value  <- as.numeric(obj$value)
    out$best_par    <- format_par(obj$par)
    out$n_par       <- length(obj$par)
    out$iterations  <- as.integer(obj$counts["function"] %||% obj$counts[1] %||% NA)
    out$convergence <- as.character(obj$convergence %||% NA)
    if (!is.null(obj$message)) out$convergence <- paste0(out$convergence, ":", obj$message)
    return(out)
  }

  # ---- Liste mit iterativem Log (dataframe mit "value"/"loss"/"obj") ----
  if (is.data.frame(obj)) {
    # Spalte fuer Zielwert finden
    cand <- intersect(c("value", "loss", "objective", "obj", "rmse", "NRMSE", "total_loss"),
                      names(obj))
    if (length(cand)) {
      col <- cand[1]
      ix  <- which.min(obj[[col]])
      out$best_value <- as.numeric(obj[[col]][ix])
      par_cols <- setdiff(names(obj), cand)
      if (length(par_cols)) {
        best_row <- obj[ix, par_cols, drop = FALSE]
        out$best_par <- paste(sprintf("%s=%.6g", names(best_row),
                                      as.numeric(unlist(best_row))),
                              collapse = "; ")
        out$n_par <- length(par_cols)
      }
      out$iterations <- nrow(obj)
      return(out)
    }
  }

  # ---- Liste, die DEoptim/optim wrappt ----
  if (is.list(obj)) {
    # Irgendein Eintrag ist selbst ein DEoptim/optim-Ergebnis?
    for (nm in names(obj)) {
      child <- obj[[nm]]
      if (is.list(child) && all(c("par", "value") %in% names(child))) {
        res <- extract_summary(child)
        res$obj_names <- paste0(out$obj_names, "/via:", nm)
        return(res)
      }
      if (is.list(child) && !is.null(child$optim) &&
          all(c("bestmem", "bestval") %in% names(child$optim))) {
        res <- extract_summary(child)
        res$obj_names <- paste0(out$obj_names, "/via:", nm)
        return(res)
      }
    }
  }

  out
}

# `%||%`-Operator (Basis-R kennt ihn erst ab 4.4)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# -----------------------------------------------------------
# 3) Schleife ueber alle Dateien
# -----------------------------------------------------------
rows <- vector("list", length(files))

for (i in seq_along(files)) {
  path <- files[i]
  meta <- parse_filename(path)
  meta$path    <- path
  meta$dataset <- meta$dataset_tag %||% infer_dataset(path)
  meta$size_kb <- round(file.info(path)$size / 1024, 1)
  meta$mtime   <- format(file.info(path)$mtime, "%Y-%m-%d %H:%M")

  obj <- tryCatch(readRDS(path),
                  error = function(e) { message("ERR ", path, ": ", e$message); NULL })

  if (is.null(obj)) {
    summary_i <- list(best_value=NA_real_, best_par=NA_character_, n_par=NA_integer_,
                      iterations=NA_integer_, convergence="READ_ERROR",
                      obj_class=NA_character_, obj_names=NA_character_)
  } else {
    summary_i <- extract_summary(obj)
  }

  rows[[i]] <- c(meta, summary_i)
  message(sprintf("[%2d/%2d] %-55s best=%-10.4g algo=%-12s",
                  i, length(files), basename(path),
                  summary_i$best_value, meta$algorithm))
}

# -----------------------------------------------------------
# 4) In data.frame umwandeln
# -----------------------------------------------------------
df <- do.call(rbind, lapply(rows, function(r) {
  as.data.frame(r, stringsAsFactors = FALSE)
}))

# Spalten ordnen
col_order <- c(
  "dataset", "algorithm",
  "w_SWE_NRMSE", "w_RHO_NRMSE", "w_SWE_NBIAS",
  "best_value", "best_par", "n_par",
  "iterations", "convergence",
  "obj_class", "obj_names",
  "file", "size_kb", "mtime", "path"
)
col_order <- intersect(col_order, names(df))
df <- df[, c(col_order, setdiff(names(df), col_order))]

# Sortieren: pro Dataset+Algo nach Gewicht
df <- df[order(df$dataset, df$algorithm,
               df$w_SWE_NRMSE, df$w_RHO_NRMSE, df$w_SWE_NBIAS), ]

# -----------------------------------------------------------
# 5) Speichern + Zusammenfassung ausgeben
# -----------------------------------------------------------
write.csv(df, OUT_CSV, row.names = FALSE)
saveRDS(df, OUT_RDS)

message("\n=== Zusammenfassung ===")
print(df[, intersect(c("dataset","algorithm","w_SWE_NRMSE","w_RHO_NRMSE",
                       "w_SWE_NBIAS","best_value","iterations","file"),
                     names(df))], row.names = FALSE)

message("\nCSV:  ", normalizePath(OUT_CSV, mustWork = FALSE))
message("RDS:  ", normalizePath(OUT_RDS, mustWork = FALSE))
