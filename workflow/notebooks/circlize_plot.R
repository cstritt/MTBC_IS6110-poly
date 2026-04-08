library(circlize)
library(ComplexHeatmap)
library(grid)

circlize_plot <- function(genome, bed_list, colors,
                          fontfamily = "Times", fontface = "italic") {
  
  circos.clear()
  circos.par(start.degree = 90)
  
  # Initialize without default axes/labels
  circos.genomicInitialize(genome, plotType = NULL)
  
  # Add sector labels and custom axes
  for (chr in get.all.sector.index()) {
    # sector label
    circos.text(
      mean(CELL_META$xlim), CELL_META$ylim[2] + convert_y(2, "mm"),
      labels = chr,
      sector.index = chr,
      facing = "bending.inside",
      niceFacing = TRUE,
      fontfamily = fontfamily,
      fontface = fontface,
      cex = 0.8
    )
    
    # axis ticks (no labels)
    circos.axis(
      h = "top",
      sector.index = chr,
      labels = FALSE,
      major.tick.length = convert_y(1.5, "mm")
    )
    
    # tick labels
    major.by = CELL_META$major.by
    breaks = seq(from = ceiling(CELL_META$xlim[1]/major.by) * major.by,
                 to   = floor(CELL_META$xlim[2]/major.by) * major.by,
                 by   = major.by)
    for (br in breaks) {
      circos.text(
        br, CELL_META$ylim[1],
        labels = as.character(br),
        sector.index = chr,
        facing = "clockwise",
        adj = c(0, 0.5),
        niceFacing = TRUE,
        fontfamily = fontfamily,
        fontface = fontface,
        cex = 0.5
      )
    }
  }
  
  # Data track
  circos.genomicTrack(
    bed_list,
    stack = TRUE,
    track.height = 0.7,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      i = getI(...)
      circos.genomicRect(region, value, col = colors[i], border = colors[i], ...)
    }
  )
  
  # Legend
  lgd = Legend(
    labels = names(colors),
    legend_gp = gpar(fill = colors),
    labels_gp = gpar(
      fontfamily = fontfamily,
      fontface = fontface,
      fontsize = 10
    ),
    title = "Legend",
    title_gp = gpar(
      fontfamily = fontfamily,
      fontface = fontface,
      fontsize = 12
    )
  )
  
  draw(lgd, x = unit(1, "npc") - unit(5, "mm"),
       y = unit(1, "npc") - unit(5, "mm"),
       just = c("right", "top"))
}
