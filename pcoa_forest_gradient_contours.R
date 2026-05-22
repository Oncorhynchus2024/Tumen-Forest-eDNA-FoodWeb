# =========================================================
# PCoA ordination with forest cover gradient contour lines
# =========================================================
#
# Purpose:
#   This script performs a Principal Coordinates Analysis (PCoA)
#   based on Bray-Curtis dissimilarities and overlays smoothed
#   contour lines of a forest cover gradient variable.
#
# Expected input files in the working directory:
#   1. spe.txt
#      - Tab-delimited species/ASV abundance table.
#      - First column should contain ASV/taxon IDs.
#      - Remaining columns should be samples/sites.
#
#   2. group.txt
#      - Tab-delimited sample grouping file.
#      - Required columns:
#          site
#          group
#
#   3. FI.txt
#      - Tab-delimited forest cover gradient file.
#      - Required columns:
#          site
#          FI
#
# Output:
#   protist_pcoa_forest_gradient_contours.pdf
#
# Notes:
#   - The exact sampling coordinates are not required for this script.
#   - The variable "FI" can be replaced by another forest gradient variable
#     if the input file and variable name are changed accordingly.
# =========================================================


# =========================================================
# 0. Install and load required packages
# =========================================================

# List of required R packages
required_packages <- c(
  "ggplot2",
  "vegan",
  "dplyr",
  "RColorBrewer",
  "mgcv",
  "metR",
  "ggnewscale"
)

# Install missing packages from CRAN
packages_to_install <- required_packages[
  !required_packages %in% installed.packages()[, "Package"]
]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install, repos = "https://cloud.r-project.org")
}

# Load packages
library(ggplot2)
library(vegan)
library(dplyr)
library(RColorBrewer)
library(mgcv)
library(metR)
library(ggnewscale)

# Set the default plotting theme
theme_set(theme_bw(base_family = "Times New Roman"))


# =========================================================
# 1. Define input and output files
# =========================================================

# Input files
species_file <- "spe.txt"
group_file <- "group.txt"
forest_gradient_file <- "FI.txt"

# Output figure
output_pdf <- "protist_pcoa_forest_gradient_contours.pdf"

# Name of the forest gradient variable in FI.txt
forest_gradient_variable <- "FI"

# Contour levels to display
forest_gradient_breaks <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

# Color scale for the forest gradient contours
# Low values are shown in yellow; high values are shown in blue-purple.
forest_gradient_colors <- c("#FDE725", "#5DC863", "#21908C", "#3B528B", "#440154")


# =========================================================
# 2. Read input data
# =========================================================

# Read the species/ASV abundance table.
# The original table is assumed to have taxa as rows and sites as columns.
species_table <- read.table(
  species_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Transpose the table so that rows represent sites and columns represent taxa.
species_table <- t(species_table)
species_table <- as.data.frame(species_table, check.names = FALSE)

# Convert all abundance values to numeric.
species_table[] <- lapply(species_table, as.numeric)

# Read site grouping information.
group_data <- read.table(
  group_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Read the forest gradient data.
forest_gradient_data <- read.table(
  forest_gradient_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# =========================================================
# 3. Check required columns
# =========================================================

# Check whether the group file contains required columns.
if (!all(c("site", "group") %in% colnames(group_data))) {
  stop("The group file must contain columns named 'site' and 'group'.")
}

# Check whether the forest gradient file contains required columns.
if (!all(c("site", forest_gradient_variable) %in% colnames(forest_gradient_data))) {
  stop(
    paste0(
      "The forest gradient file must contain columns named 'site' and '",
      forest_gradient_variable,
      "'."
    )
  )
}


# =========================================================
# 4. Align species data with metadata
# =========================================================

# Create a site list from the row names of the species table.
site_data <- data.frame(
  site = rownames(species_table),
  stringsAsFactors = FALSE
)

# Merge site IDs with group and forest gradient information.
metadata <- site_data %>%
  left_join(group_data, by = "site") %>%
  left_join(forest_gradient_data, by = "site")

# Remove samples with missing group or forest gradient information.
if (
  any(is.na(metadata$group)) ||
  any(is.na(metadata[[forest_gradient_variable]]))
) {
  warning(
    "Some sites are missing group or forest gradient values and were removed."
  )

  metadata <- metadata %>%
    filter(
      !is.na(group),
      !is.na(.data[[forest_gradient_variable]])
    )
}

# Reorder the species table to match the metadata.
species_table <- species_table[metadata$site, , drop = FALSE]

# Confirm that sample order is identical between the abundance table and metadata.
stopifnot(all(rownames(species_table) == metadata$site))


# =========================================================
# 5. Calculate Bray-Curtis dissimilarity
# =========================================================

# Bray-Curtis dissimilarity is commonly used for community composition data.
bray_curtis_distance <- vegdist(species_table, method = "bray")


# =========================================================
# 6. Run Principal Coordinates Analysis
# =========================================================

# Perform PCoA using classical multidimensional scaling.
pcoa_result <- cmdscale(
  bray_curtis_distance,
  eig = TRUE,
  k = 2
)

# Extract PCoA site scores.
pcoa_scores <- as.data.frame(pcoa_result$points)
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
pcoa_scores$site <- rownames(pcoa_scores)

# Add group and forest gradient information to the PCoA scores.
pcoa_scores <- pcoa_scores %>%
  left_join(metadata, by = "site")


# =========================================================
# 7. Calculate axis explanatory percentages
# =========================================================

# Only positive eigenvalues are used as the denominator.
eigenvalues <- pcoa_result$eig
positive_eigenvalues <- eigenvalues[eigenvalues > 0]

pcoa1_explained <- eigenvalues[1] / sum(positive_eigenvalues) * 100
pcoa2_explained <- eigenvalues[2] / sum(positive_eigenvalues) * 100

x_axis_label <- sprintf("PCoA1 (%.1f%%)", pcoa1_explained)
y_axis_label <- sprintf("PCoA2 (%.1f%%)", pcoa2_explained)


# =========================================================
# 8. Calculate convex hulls for each group
# =========================================================

# Convex hulls are used to outline the distribution range of each group.
hull_points <- pcoa_scores %>%
  group_by(group) %>%
  slice(chull(PCoA1, PCoA2)) %>%
  ungroup()


# =========================================================
# 9. Fit a smooth surface for the forest gradient
# =========================================================

# Choose the basis dimension for the GAM smooth term.
# It should be smaller than the number of samples.
smooth_k <- min(10, nrow(pcoa_scores) - 1)
smooth_k <- max(4, smooth_k)

# Fit a two-dimensional GAM surface using PCoA coordinates.
forest_gradient_gam <- gam(
  as.formula(
    paste0(
      forest_gradient_variable,
      " ~ s(PCoA1, PCoA2, k = ",
      smooth_k,
      ")"
    )
  ),
  data = pcoa_scores,
  method = "REML"
)

# Build a prediction grid across the PCoA ordination space.
x_sequence <- seq(
  min(pcoa_scores$PCoA1) - 0.05,
  max(pcoa_scores$PCoA1) + 0.05,
  length.out = 300
)

y_sequence <- seq(
  min(pcoa_scores$PCoA2) - 0.05,
  max(pcoa_scores$PCoA2) + 0.05,
  length.out = 300
)

prediction_grid <- expand.grid(
  PCoA1 = x_sequence,
  PCoA2 = y_sequence
)

# Predict forest gradient values on the grid.
prediction_grid$forest_gradient_predicted <- predict(
  forest_gradient_gam,
  newdata = prediction_grid
)

# Restrict predicted values to the observed range to reduce unrealistic extrapolation.
observed_forest_gradient_range <- range(
  pcoa_scores[[forest_gradient_variable]],
  na.rm = TRUE
)

prediction_grid$forest_gradient_predicted <- pmin(
  pmax(
    prediction_grid$forest_gradient_predicted,
    observed_forest_gradient_range[1]
  ),
  observed_forest_gradient_range[2]
)


# =========================================================
# 10. Create the PCoA plot with forest gradient contours
# =========================================================

pcoa_plot <- ggplot() +

  # Forest gradient contour lines.
  metR::geom_contour2(
    data = prediction_grid,
    aes(
      x = PCoA1,
      y = PCoA2,
      z = forest_gradient_predicted,
      color = after_stat(level)
    ),
    breaks = forest_gradient_breaks,
    linewidth = 0.9
  ) +

  # Text labels for contour lines.
  metR::geom_text_contour(
    data = prediction_grid,
    aes(
      x = PCoA1,
      y = PCoA2,
      z = forest_gradient_predicted,
      label = after_stat(sprintf("%.1f", level)),
      color = after_stat(level)
    ),
    breaks = forest_gradient_breaks,
    size = 4,
    check_overlap = TRUE
  ) +

  # Continuous color scale for the forest gradient.
  scale_color_gradientn(
    colours = forest_gradient_colors,
    breaks = forest_gradient_breaks,
    limits = c(
      min(forest_gradient_breaks),
      max(forest_gradient_breaks)
    ),
    name = "Forest cover gradient",
    guide = guide_colorbar(order = 2)
  ) +

  # Start a new independent color scale for the sample groups.
  ggnewscale::new_scale_color() +

  # Convex hulls for each group.
  geom_polygon(
    data = hull_points,
    aes(
      x = PCoA1,
      y = PCoA2,
      group = group,
      color = group
    ),
    fill = NA,
    linetype = "dashed",
    linewidth = 0.8
  ) +

  # Site points.
  geom_point(
    data = pcoa_scores,
    aes(
      x = PCoA1,
      y = PCoA2,
      color = group
    ),
    size = 3,
    shape = 17
  ) +

  # Discrete group color scale.
  scale_color_brewer(
    palette = "Dark2",
    name = "Group",
    guide = guide_legend(order = 1)
  ) +

  # Reference lines crossing the origin.
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "black"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "black"
  ) +

  # Axis labels and final layout settings.
  labs(
    x = x_axis_label,
    y = y_axis_label
  ) +
  coord_equal() +
  theme(
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

# Display the figure in the R plotting window.
print(pcoa_plot)


# =========================================================
# 11. Export the figure
# =========================================================

ggsave(
  filename = output_pdf,
  plot = pcoa_plot,
  width = 8,
  height = 6
)
