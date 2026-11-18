# Independent numerical anchors and an extracted-download report render.
library(readr)
library(dplyr)
assert_close <- function(actual, expected, tolerance = 1e-8) {
  stopifnot(isTRUE(all.equal(as.numeric(actual), as.numeric(expected),
                           tolerance = tolerance)))
}
amounts <- read_csv("data/amounts.csv", show_col_types = FALSE)
assert_close(amounts$mass_g / amounts$molar_mass_g_mol, c(0.010, 0.015))
extent <- min(amounts$mass_g / amounts$molar_mass_g_mol / amounts$stoich)
assert_close(extent * 180.16, 1.8016)
assert_close(100 * 1.20 / (extent * 180.16), 66.60746003552397)
assert_close(1 / (1 + 10^(4.8 - c(3.8,4.8,5.8))), c(1/11,1/2,10/11))
tlc <- read_csv("data/tlc.csv", show_col_types = FALSE)
assert_close(tlc$spot_distance_cm / tlc$front_distance_cm, c(0.24,0.64,0.24,0.64))
ee <- read_csv("data/enantiomers.csv", show_col_types = FALSE)
assert_close(abs(ee$R_percent - ee$S_percent), c(40,0,90))
conf <- read_csv("data/conformers.csv", show_col_types = FALSE)
w <- exp(-conf$delta_G_J_mol / (8.314462618 * 298.15))
assert_close(sum(w/sum(w)), 1)
stopifnot(w[1] > w[2], all(w > 0))
ir <- read_csv("data/ir.csv", show_col_types = FALSE)
stopifnot(all(ir$transmittance_percent >= 0), all(ir$transmittance_percent <= 100),
          length(unique(ir$sample_id)) == 4)

project_dir <- tempfile("orgchem-student-")
dir.create(project_dir)
unzip("downloads/organic-chemistry-student.zip", exdir = project_dir)
report_path <- file.path(project_dir,"report.Rmd")
# render() knits relative to the Rmd's own folder, proving the extracted folder
# contains all required inputs rather than relying on this checkout's data.
expression <- paste0("rmarkdown::render(", deparse(report_path), ", quiet=TRUE)")
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c("--vanilla", "-e", shQuote(expression)))
stopifnot(status == 0)
result <- file.path(project_dir,"report.html")
stopifnot(file.exists(result), file.info(result)$size > 10000)
unlink(project_dir,recursive=TRUE)
cat("Numerical anchors and a clean-session report from the extracted student ZIP pass.\n")
