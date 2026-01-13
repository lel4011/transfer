#' Create a multipanel theme compatible with ggplot2 4.0.0 S7 system
#'
#' @description
#' Creates a multipanel theme that is compatible with the new ggplot2 4.0.0 S7 system.
#' This function creates a fresh theme object using the new S7 system instead of relying
#' on the saved S3 theme object.
#'
#' @return A `ggplot2` theme object compatible with ggplot2 4.0.0
#'
#' @export
sccomp_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.background = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), size = 7, colour = "black"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), size = 7, colour = "black"),
      panel.spacing.x = ggplot2::unit(0.1, "lines"),
      axis.text.x = ggplot2::element_text(size = 6, colour = "grey30"),
      axis.text.y = ggplot2::element_text(size = 6, colour = "grey30"),
      strip.text.x = ggplot2::element_text(size = 7, colour = "grey10"),
      strip.text.y = ggplot2::element_text(size = 7, colour = "grey10"),
      legend.key.size = ggplot2::unit(5, 'mm'),
      legend.title = ggplot2::element_text(size = 7, colour = "black"),
      legend.text = ggplot2::element_text(size = 6, colour = "black"),
      strip.clip = "off",
      plot.title = ggplot2::element_text(size = 7, colour = "black"),
      axis.line.x = ggplot2::element_line(linewidth = 0.2, colour = "black"),
      axis.line.y = ggplot2::element_line(linewidth = 0.2, colour = "black"),
      axis.ticks.x = ggplot2::element_line(linewidth = 0.2, colour = "grey20"),
      axis.ticks.y = ggplot2::element_line(linewidth = 0.2, colour = "grey20"),
      axis.ticks.length = ggplot2::unit(2.75, "pt")
    )
}
