#!/usr/bin/env Rscript

# Common parameters for all plots

# Figure width
width_cm = 20
width_in = 20

# Figure height
height_cm = 10
height_in = 13

# Figure resolution
dpi_no = 300

# Justification
plot_title_justification <- 0.5

# Title sizes
axis_title_rel_size   <- 1.25
plot_title_rel_size   <- 1.5

# theme_own <- theme(
#     #text = element_text(family = "serif", size = 14),
#     rect = element_blank(),
#     # panel.grid = element_blank(),
#     #title = element_text(color = "#8b0000"),
#     axis.line = element_line(color = "black")
#     )

# ggtitle(title) +
#   theme(plot.title = element_text(hjust = plot_title_justification, 
#                                   size  = plot_title_size), 
#         axis.title=element_text(  size  = axis_title_size)) + 
#   theme(legend.text=element_text( size  = legend_text_size)) +
#   theme(legend.title=element_text(size  = legend_title_size)) +
#   theme(axis.text=element_text(   size  = axis_text_size))