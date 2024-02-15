# Function to print help and usage instructions.
print_help <- function() {
  cat("
Usage: Rscript density_plotter.R <filename> <chrom_string> [merge] <threshold> [no_size_cutoff]

<filename>      : Path to the TE superfamily density table generated using density_table.py.
<chrom_string>  : Optional. Specifies a string pattern to include in the analysis. Chromosome IDs without this string will be excluded. Useful for excluding scaffolds, contigs, or off-target chromosomes.
[merge]         : Optional. Boolean flag. If 'merge' is specified, all chromosome plots are merged onto a single page.
<threshold>     : Optional. Exclusionary threshold for plotting repeat types. Types that do not exceed this density percentage at any position will be excluded [00-100]. ie, '02' would exclude types that do not exceed 2% density in any window.
[no_size_cutoff]: Optional. Boolean flag. If 'no_size_cutoff' is specified, all chromosomes and scaffolds will be plotted. 

This script plots the distribution of EDTA annotations across chromosomes.

Dependencies: ggplot2, dplyr for plotting. R and its packages can be installed with conda.
For installation:
conda config --add channels conda-forge && conda config --add channels CRAN && conda create -n r_env r-base r-ggplot2 r-dplyr

")
  quit(status = 0)
}

# Get command line arguments.
args <- commandArgs(trailingOnly = TRUE)

# Initialize variables.
chrom_string <- NULL
merge_plots <- FALSE
threshold <- NULL
no_size_cutoff <- FALSE 
filename <- NULL

# Check arguments using a switch statement.
for (arg in args) {
  switch(tolower(arg),
    "merge" = { merge_plots <- TRUE },
    "no_size_cutoff" = { no_size_cutoff <- TRUE },
    {
      # Check if the argument is a number (for threshold).
      numeric_arg <- as.numeric(arg)
      if (!is.na(numeric_arg)) {
        threshold <- numeric_arg
      } else if (is.null(filename)) {
        # Assume the first non-flag, non-numeric argument is the filename.
        filename <- arg
      } else {
        # If not a number and filename is already set, it's assumed to be the chrom_string.
        chrom_string <- arg
      }
    }
  )
}

# Check if filename is not provided or help is requested.
if (is.null(filename) || filename == '-h' || filename == '--help') {
  print_help()
}

filename <- args[1]

# Check if the file has a header.
has_header <- readLines(filename, n = 1)
header <- if (grepl("chrom", has_header, fixed = TRUE)) TRUE else FALSE

repeats <- read.table(filename, sep='\t', header = header, 
                      col.names = if (!header) c('chrom', 'coord', 'density', 'type') else NULL)
                      
# Calculate the 10% size cutoff and apply it if no_size_cutoff is not set.
# This will exclude plotting sequences that are less than 10% the length of the longest chromosome. 
# Aimed at excluding contigs and scaffolds. 
if (!no_size_cutoff) {
  chrom_lengths <- aggregate(coord ~ chrom, data = repeats, max)
  max_length <- max(chrom_lengths$coord)
  length_threshold <- max_length * 0.1
  valid_chroms <- chrom_lengths$chrom[chrom_lengths$coord >= length_threshold]
  repeats <- subset(repeats, chrom %in% valid_chroms)
}

# Exclude specific types.
excluded_types <- c('target_site_duplication', 'repeat_region', 'long_terminal_repeat')
repeats <- subset(repeats, !type %in% excluded_types)

# Replace 'Copia' with 'Ty1' and 'Gypsy' with 'Ty3'.
library(dplyr)
repeats <- repeats %>%
    mutate(type = gsub("Copia_LTR_retrotransposon", "Ty1_LTR_retrotransposon", type)) %>%
    mutate(type = gsub("Gypsy_LTR_retrotransposon", "Ty3_LTR_retrotransposon", type))

# Reorder the type factor to move 'Ty1_LTR_retrotransposon' higher up.
# This is so that Ty1 and Ty3 dont end up plotted as similar colors.
repeats$type <- factor(repeats$type, levels = c("Ty1_LTR_retrotransposon", setdiff(unique(repeats$type), "Ty1_LTR_retrotransposon")))

# Diagnostic print statement.
cat("Status of no_size_cutoff flag: ", no_size_cutoff, "\n")

# Optional chromosome string filtering.
if (!is.null(chrom_string)) {
  repeats <- subset(repeats, grepl(chrom_string, chrom))
}

# Apply the exclusionary threshold.
if (!is.null(threshold)) {
  max_density <- aggregate(density ~ type, data = repeats, max)
  included_types <- max_density[max_density$density * 100 > threshold, "type"]
  repeats <- subset(repeats, type %in% included_types)
}

# Convert chromosome names to a factor with the correct order.
chrom_levels <- unique(repeats$chrom)
chrom_levels_sorted <- chrom_levels[order(as.numeric(gsub("[^0-9]", "", chrom_levels)))]
repeats$chrom <- factor(repeats$chrom, levels = chrom_levels_sorted)

# Plotting.
library(ggplot2)

# Define a shuffled color palette using hcl.colors with the 'Dark 3' palette.
num_categories <- length(unique(repeats$type))
set.seed(46) # Change the seed for a new shuffle of colors. 
color_palette <- sample(hcl.colors(num_categories, "Dark 3"))

plot_border_theme <- theme(
  panel.border = element_rect(colour = "grey", fill=NA, linewidth=1),
  axis.text.x = element_text(),
  axis.ticks.x = element_line(),
  plot.margin = margin(1, 1, 1, 1, "lines") # Adjust plot margin.
)

if (merge_plots) {

  p <- ggplot(repeats, aes(x = coord / 1e6, y = density * 100, color = type)) + 
       geom_line() +
       facet_wrap(~ chrom, scales = "free_x") + 
       labs(title = "Density of Repeat Types on Chromosomes (%)", 
            x = "Position on chromosome (Mb)", 
            y = "Percent repeat type in window (%)", 
            color = "Type") + 
       theme_minimal() +
       plot_border_theme +
       scale_color_manual(values = color_palette) +
       scale_x_continuous(expand = c(0, 0)) + # Tighten grey plot border for x-axis.
       scale_y_continuous(expand = expansion(mult = c(0.005, 0.01))) # Tighten grey plot border for y-axis.

  # Save the plot to a PDF file.
  pdf("chromosome_density_plots_merged.pdf", width = 16, height = 6)
  print(p)
  dev.off()
  
} else {
  
  # Unmerged plot logic.
  pdf("chromosome_density_plots.pdf", width = 12, height = 4)
  for (chr in unique(repeats$chrom)) {
    subset_df <- repeats[repeats$chrom == chr,]
    
    if (nrow(subset_df) > 1) {
      p <- ggplot(subset_df, aes(x = coord / 1e6, y = density * 100, group = interaction(chrom, type), color = type)) + 
           geom_line() + 
           labs(title = paste("Density of Repeat Types on", chr, "(%)"), 
                x = "Position on chromosome (Mb)", 
                y = "Percent repeat type in window (%)", 
                color = "Type") +
           theme_minimal() +
           plot_border_theme +
           scale_color_manual(values = color_palette) +
           scale_x_continuous(expand = c(0, 0)) + # Tighten grey plot border for x-axis.
           scale_y_continuous(expand = expansion(mult = c(0.005, 0.01))) 
      print(p)
    } else {
      cat(paste("Not enough data points to plot for chromosome:", chr, "\n"))
    }
  }
  dev.off()
}
