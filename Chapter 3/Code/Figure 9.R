library(ggplot2)

# Read the actual Figure 9 data
d <- read.csv("quantile_error_summary.csv")

# Set the ordering used in the paper
d$dose <- factor(
  d$dose,
  levels = c(1, 5, 10, 50, 100),
  labels = paste0("dose: ", c(1, 5, 10, 50, 100))
)

d$model <- factor(
  d$model,
  levels = c(
    "ODE-only",
    "SDE-only",
    "MTBP-to-ODE hybrid",
    "MTBP-to-SDE-to-ODE hybrid"
  ),
  labels = c(
    "ODE only",
    "SDE only",
    "MTBP--ODE",
    "MTBP--SDE--ODE"
  )
)

d$quantile <- factor(
  d$quantile,
  levels = c(0.025, 0.05, 0.25, 0.5,
             0.75, 0.95, 0.975)
)

# Reproduce Figure 9
fig9 <- ggplot(
  d,
  aes(x = quantile, y = model, fill = error_days)
) +
  geom_tile(
    colour = "white",
    linewidth = 0.25
  ) +
  facet_grid(. ~ dose) +
  scale_fill_gradient2(
    low = "#98332F",
    mid = "white",
    high = "#B6ADD2",
    midpoint = 0,
    limits = c(-4, 2),
    breaks = c(-3, -2, -1, 0, 1),
    name = "Error\n(days)"
  ) +
  labs(
    x = "Reference distribution quantile",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(
      angle = 45, hjust = 1, size = 9
    ),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_blank(),
    panel.spacing.x = grid::unit(0.15, "lines"),
    legend.key.height = grid::unit(1.3, "cm")
  )

print(fig9)

ggsave(
  "Figure_9.png",
  fig9,
  width = 12,
  height = 4.2,
  dpi = 300
)

ggsave(
  "Figure_9.pdf",
  fig9,
  width = 12,
  height = 4.2
)
