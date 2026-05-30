`%||%` <- function(a, b) {
	if (is.null(a) || length(a) == 0) b else a
}

files <- c(
	"opt_results__SWE_NRMSE_0p3__RHO_NRMSE_0p5__SWE_NBIAS_0p2.rds",
	"opt_results__SWE_NRMSE_0p3__RHO_NRMSE_0p7__SWE_NBIAS_0p0.rds",
	"opt_results__SWE_NRMSE_0p5__RHO_NRMSE_0p5__SWE_NBIAS_0p0.rds",
	"opt_results__SWE_NRMSE_0p6__RHO_NRMSE_0p2__SWE_NBIAS_0p2.rds",
	"opt_results__SWE_NRMSE_0p7__RHO_NRMSE_0p0__SWE_NBIAS_0p3.rds",
    
	"opt_results__SWE_NRMSE_0p7__RHO_NRMSE_0p3__SWE_NBIAS_0p0.rds",
	"opt_results__SWE_NRMSE_1__RHO_NRMSE_0__SWE_NBIAS_0.rds"
)

script_file <- sys.frame(1)$ofile %||% ""
base_dir <- if (nzchar(script_file)) {
	dirname(normalizePath(script_file, mustWork = FALSE))
} else {
	getwd()
}

show_df <- function(df, label) {
	cat("\n", strrep("=", 90), "\n", sep = "")
	cat("Object:", label, "\n")
	cat("Rows:", nrow(df), " Columns:", ncol(df), "\n")
	print(utils::head(df, 10))
	if (interactive()) {
		utils::View(df, title = label)
	}
}

for (f in files) {
	path <- file.path(base_dir, f)

	if (!file.exists(path)) {
		warning("Missing file: ", path)
		next
	}

	obj <- readRDS(path)
	obj_name <- tools::file_path_sans_ext(basename(f))

	if (is.data.frame(obj)) {
		show_df(obj, obj_name)
	} else if (is.list(obj)) {
		cat("\n", strrep("-", 90), "\n", sep = "")
		cat("File:", basename(f), "contains a list with", length(obj), "elements\n")

		for (i in seq_along(obj)) {
			elem <- obj[[i]]
			elem_name <- names(obj)[i] %||% paste0("element_", i)
			label <- paste0(obj_name, "::", elem_name)

			if (is.data.frame(elem)) {
				show_df(elem, label)
			} else {
				cat("\nObject:", label, "(class:", paste(class(elem), collapse = ", "), ")\n")
				print(utils::str(elem))
			}
		}
	} else {
		cat("\n", strrep("-", 90), "\n", sep = "")
		cat("File:", basename(f), "\n")
		cat("Class:", paste(class(obj), collapse = ", "), "\n")
		print(utils::str(obj))
	}
}
