# -------------------------------------------------------------------------
# Estimated Recovery Durations
# -------------------------------------------------------------------------
# August 2026
# Graphing evidence of different recovery durations and our best guess estimate
# Trying to demonstrate the consideration that went into what duration was
# assigned to the min and max recovery index
# -------------------------------------------------------------------------

# Load packages -----------------------------------------------------------

library(ggplot2)
library(patchwork)

# Min and Max Recovery Index Scores ---------------------------------------

x_min <- 0.15895
x_max <- 3.25544

x_values <- seq(
  x_min,
  x_max,
  length.out = 200
)

# Our best guess estimate of recovery times -------------------------------

# Range from 1 month to 36 months (3 years)

estimate_data <- data.frame(
  Recovery_Index = x_values,
  Recovery_Duration = seq(
    1,
    36,
    length.out = length(x_values)
  )
)

# Evidence of different recovery durations -------------------------------

# Evidence has been coloured by the dimensions of recovery they relate to.

evidence_dots <- data.frame(
  Recovery_Index = c(
    x_min,
    x_min,
    x_min,
    x_min,
    x_min,
    x_min,
    x_max,
    x_max,
    x_max,
    x_min,
    x_max,
    x_max,
    x_max,
    x_max,
    x_max,
    x_max,
    x_max,
    x_max
  ),
  Recovery_Duration = c(
    1,
    3,
    4,
    7,
    10,
    10,
    12,
    12,
    12,
    16,
    18,
    24,
    24,
    24,
    36,
    36,
    48,
    84
  ),
  Dimension_of_Recovery = c(
    "Roading",
    "Roading",
    "Roading",
    "Health",
    "Insurance",
    "Roading",
    "Health",
    "Agriculture",
    "Roading",
    "Insurance",
    "Recovery Projects",
    "Recovery Projects",
    "Agriculture",
    "Health",
    "Recovery Projects",
    "Recovery Projects",
    "Roading",
    "Roading"
  ),
  No_ID = c(
    1,
    3,
    2,
    4,
    5,
    6,
    8,
    9,
    10,
    7,
    20,
    11,
    12,
    13,
    14,
    15,
    16,
    17
  )
)

# Create fixed jittered positions -----------------------------------------
# The same jittered coordinates are used for both the dots and their IDs.
# so each number remains centred on its corresponding dot.

set.seed(123)
evidence_dots$Plot_X <- evidence_dots$Recovery_Index +
  runif(
    nrow(evidence_dots),
    -0.1,
    0.1
  )

evidence_dots$Plot_Y <- evidence_dots$Recovery_Duration +
  runif(
    nrow(evidence_dots),
    -1.5,
    1.5
  )
# Colours -----------------------------------------------------------------

recovery_colours <- c(
  "Roading" = "#E69F00",
  "Health" = "darkgreen",
  "Recovery Projects" = "#CC79A7",
  "Agriculture" = "#D55E00",
  "Insurance" = "#56B4E9"
)

# Main Recovery Duration Plot ---------------------------------------------

recovery_duration_plot <- ggplot(
  estimate_data,
  aes(
    x = Recovery_Index,
    y = Recovery_Duration
  )
) +
  
  # Estimated recovery duration ------------------------------------------

geom_line(
  linewidth = 1.2,
  colour = "black"
) +
  
  # Evidence points -------------------------------------------------------

geom_point(
  data = evidence_dots,
  aes(
    x = Plot_X,
    y = Plot_Y,
    colour = Dimension_of_Recovery
  ),
  inherit.aes = FALSE,
  size = 5
) +
  
  # Evidence IDs ----------------------------------------------------------

geom_text(
  data = evidence_dots,
  aes(
    x = Plot_X,
    y = Plot_Y,
    label = No_ID
  ),
  inherit.aes = FALSE,
  colour = "black",
  size = 3.2,
  fontface = "bold"
) +
  
  
  # X-axis ----------------------------------------------------------------

scale_x_continuous(
  limits = c(0, 4),
  breaks = seq(0, 4, 1),
  expand = c(0, 0)
) +
  
  # Y-axis ----------------------------------------------------------------

scale_y_continuous(
  limits = c(0, 96),
  breaks = seq(0, 96, 12),
  expand = c(0, 0)
) +
  
  # Evidence colours ------------------------------------------------------

scale_colour_manual(
  values = recovery_colours
) +
  
  # Axis labels ------------------------------------------------------------

labs(
  x = "Recovery Index",
  y = "Recovery Duration (Months)"
) +
  
  # Theme -----------------------------------------------------------------

theme_classic() +
  
  theme(
    # Axis text
    axis.text = element_text(
      size = 16,
      colour = "black"
    ),
    
    # Axis titles
    axis.title = element_text(
      size = 17,
      colour = "black"
    ),
    
    # Axis lines
    axis.line = element_line(
      colour = "black",
      linewidth = 0.8
    ),
    
    # Legend in upper-right of main plot
    legend.position = c(
      0.2,
      0.8
    ),
    
    # Remove legend title
    legend.title = element_blank(),
    
    # Legend text
    legend.text = element_text(
      size = 13,
      colour = "black"
    ),
    
    # White legend background with black outline
    legend.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.7
    ),
    
    # Spacing inside legend
    legend.margin = margin(
      6,
      8,
      6,
      8
    ),
    
    # Spacing between legend items
    legend.spacing.x = unit(
      0.15,
      "cm"
    )
    
  )

# Absolute Value Secondary Axis -------------------------------------------

absolute_axis <- ggplot() +
  
  geom_segment(
    aes(
      x = x_min,
      xend = x_max,
      y = 0,
      yend = 0
    ),
    linewidth = 0.8,
    colour = "black"
  ) +
  
  geom_segment(
    aes(
      x = c(x_min, x_max),
      xend = c(x_min, x_max),
      y = 0,
      yend = -0.08
    ),
    linewidth = 0.8,
    colour = "black"
  ) +
  
  geom_text(
    aes(
      x = c(x_min, x_max),
      y = -0.17,
      label = sprintf(
        "%.2f",
        c(x_min, x_max)
      )
    ),
    size = 4.5,
    colour = "black"
  ) +
  
  annotate(
    "text",
    x = x_max + 0.10,
    y = 0,
    label = "Absolute Value",
    hjust = 0,
    vjust = 0.5,
    size = 4.5,
    colour = "black"
  ) +
  
  scale_x_continuous(
    limits = c(0, 4),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(-0.27, 0.05),
    expand = c(0, 0)
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  theme_void()

# Normalised Secondary Axis -----------------------------------------------

normalised_axis <- ggplot() +
  
  geom_segment(
    aes(
      x = x_min,
      xend = x_max,
      y = 0,
      yend = 0
    ),
    linewidth = 0.8,
    colour = "black"
  ) +
  
  geom_segment(
    aes(
      x = c(x_min, x_max),
      xend = c(x_min, x_max),
      y = 0,
      yend = -0.08
    ),
    linewidth = 0.8,
    colour = "black"
  ) +
  
  geom_text(
    aes(
      x = c(x_min, x_max),
      y = -0.17,
      label = c(
        "0",
        "1"
      )
    ),
    size = 4.5,
    colour = "black"
  ) +
  
  annotate(
    "text",
    x = x_max + 0.10,
    y = 0,
    label = "Normalised",
    hjust = 0,
    vjust = 0.5,
    size = 4.5,
    colour = "black"
  ) +
  
  scale_x_continuous(
    limits = c(0, 4),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(-0.27, 0.05),
    expand = c(0, 0)
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  theme_void()

# Combine plots -----------------------------------------------------------

recovery_duration_plot <-
  recovery_duration_plot /
  absolute_axis /
  normalised_axis +
  plot_layout(
    heights = c(
      8,
      1,
      1
    )
  )

# Display plot ------------------------------------------------------------

recovery_duration_plot

# Save plot ----------------------------------------------------------------

ggsave(
  filename = "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Recovery Index/Estimated Recovery Duration.png",
  plot = recovery_duration_plot,
  width = 14,
  height = 9,
  units = "in",
  dpi = 600,
  bg = "white"
)

