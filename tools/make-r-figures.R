# Regenerate plots from the teaching tables. No data are measured observations.
library(readr)
library(dplyr)
library(ggplot2)
dir.create("assets", showWarnings = FALSE)
theme_set(theme_minimal(base_size = 15))
ion <- data.frame(pH = seq(1, 9, by = 0.05)) |>
  mutate(fraction = 1 / (1 + 10^(4.8 - pH)))
p <- ggplot(ion, aes(pH, fraction)) +
  geom_line(colour = "#176b91", linewidth = 1.2) +
  geom_vline(xintercept = 4.8, linetype = "dashed") +
  labs(x = "pH", y = "Fraction present as A−",
       title = "Aqueous ionisation model", subtitle = "Supplied pKa = 4.8; dilute-solution approximation")
ggsave("assets/ionisation.png", p, width = 8, height = 4.5, dpi = 160)
ir <- read_csv("data/ir.csv", show_col_types = FALSE) |>
  filter(sample_id %in% c("ethanol", "ethyl_ethanoate"))
p <- ggplot(ir, aes(wavenumber_cm1, transmittance_percent)) +
  geom_line(colour = "#176b91") + facet_wrap(~sample_id, ncol = 1) +
  scale_x_reverse() + labs(x = expression("Wavenumber ("*cm^{-1}*")"),
       y = "Transmittance (%)", title = "Simulated teaching IR traces")
ggsave("assets/ir-comparison.png", p, width = 8, height = 5.5, dpi = 160)
cat("Built two R plots in assets/.\n")
