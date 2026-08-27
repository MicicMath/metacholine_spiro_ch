packageManage = function(...) {
  libs = unlist(list(...))
  req = unlist(lapply(libs, require, character.only = TRUE))
  has = libs[req == TRUE]
  lapply(has, library, character.only = TRUE)
  need = libs[req == FALSE]
  if (length(need) > 0) {
    install.packages(need)
    lapply(need, library, character.only = TRUE)
  }
}

packageManage(
  "ggplot2",
  "gridExtra",
  "ggrepel",
  "scales",
  "grid",
  "rspiro",
  "glmnet",
  "ranger",
  "e1071",
  "parallel",
  "boot",
  "dplyr",
  "MASS",
  "car",
  "caret",
  "Boruta",
  "DescTools",
  "performance",
  "jtools",
  "circlize",
  "qvalue",
  "dendextend",
  "viridis",
  "ggsci",
  "cowplot",
  "ggpubr",
  "rstatix",
  "reshape2",
  "ROCR"
)

#######
# setup
#######

p_rnd = function(x) {
  p_val = ifelse(x < 0.001,
    "<0.001",
    sprintf("%.3f", round(x, 3))
  )
}

n_rnd = function(x, n = 1) {
  form = paste0("%.", n, "f")
  sprintf(form, x)
}

sig_symbol = function(p) {
  if (p > 0.05) {
    return("ns")
  }
  if (p <= 0.05 & p > 0.01) {
    return("*")
  }
  if (p <= 0.01 & p > 0.001) {
    return("**")
  }
  if (p <= 0.001 & p > 0.0001) {
    return("***")
  }
  if (p <= 0.0001) {
    return("****")
  }
}

eff_interpretation = function(x) {
  if (x < 0.1) {
    return("negligible")
  }
  if (x >= 0.1 & x < 0.3) {
    return("small")
  }
  if (x >= 0.3 & x < 0.5) {
    return("medium")
  }
  if (x >= 0.5) {
    return("large")
  }
}

ggtheme_no_legend = ggplot2::theme(
  legend.position = "none",
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(),
  plot.title = ggplot2::element_text(size = 14, hjust = 0),
  plot.subtitle = ggplot2::element_text(size = 12),
  axis.title = ggplot2::element_text(size = 12),
  axis.text.x = ggplot2::element_text(size = 11),
  axis.text.y = ggplot2::element_text(size = 11)
)

ggtheme_no_legend_45 = ggplot2::theme(
  legend.position = "none",
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(),
  plot.title = ggplot2::element_text(size = 14, hjust = 0),
  plot.subtitle = ggplot2::element_text(size = 12),
  axis.title = ggplot2::element_text(size = 12),
  axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 11),
  axis.text.y = ggplot2::element_text(size = 11)
)

ggtheme_bottom_right = theme(
  legend.position = "inside",
  legend.position.inside = c(.99, .01),
  legend.justification = c("right", "bottom"),
  # legend.background = element_blank(),
  legend.background = element_rect(fill = "white"),
  legend.spacing.y = unit(0.2, "cm"),
  legend.key.size = unit(.2, "cm"),
  legend.box.background = element_rect(colour = "black", linewidth = 0.0),
  legend.key.height = unit(0.1, "cm"),
  legend.title = element_blank(),
  legend.text = element_text(size = 9, hjust = 0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_text(size = 14, hjust = 0),
  plot.subtitle = element_text(size = 12, hjust = 0.5),
  axis.title = element_text(size = 12),
  axis.text.x = element_text(size = 11),
  axis.text.y = element_text(size = 11)
)

ggtheme_bottom_left = theme(
  legend.position = "inside",
  legend.position.inside = c(.01, .01),
  legend.justification = c("left", "bottom"),
  # legend.background = element_blank(),
  legend.background = element_rect(fill = "white"),
  legend.spacing.y = unit(0.2, "cm"),
  legend.key.size = unit(.2, "cm"),
  legend.box.background = element_rect(colour = "black", linewidth = 0.0),
  legend.key.height = unit(0.1, "cm"),
  legend.title = element_blank(),
  legend.text = element_text(size = 9, hjust = 0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_text(size = 14, hjust = 0),
  plot.subtitle = element_text(size = 12, hjust = 0.5),
  axis.title = element_text(size = 12),
  axis.text.x = element_text(size = 11),
  axis.text.y = element_text(size = 11)
)

ggtheme_top_right = theme(
  legend.position = "inside",
  legend.position.inside = c(.99, .99),
  legend.justification = c("right", "top"),
  # legend.background = element_blank(),
  legend.background = element_rect(fill = "white"),
  legend.spacing.y = unit(.5, "cm"),
  legend.key.size = unit(.2, "cm"),
  legend.box.background = element_rect(colour = "black", linewidth = 0.2),
  legend.key.height = unit(0.1, "cm"),
  legend.title = element_blank(),
  legend.text = element_text(size = 9, hjust = 0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_text(size = 14, hjust = 0),
  plot.subtitle = element_text(size = 12, hjust = 0.5),
  axis.title = element_text(size = 12),
  axis.text.x = element_text(size = 11),
  axis.text.y = element_text(size = 11)
)

ggtheme_top_left = theme(
  legend.position = "inside",
  legend.position.inside = c(.01, .99),
  legend.justification = c("left", "top"),
  # legend.background = element_blank(),
  legend.background = element_rect(fill = "white"),
  legend.spacing.y = unit(0.2, "cm"),
  legend.key.size = unit(.2, "cm"),
  legend.box.background = element_rect(colour = "black", linewidth = 0.0),
  legend.key.height = unit(0.1, "cm"),
  legend.title = element_blank(),
  legend.text = element_text(size = 9, hjust = 0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_text(size = 14, hjust = 0),
  plot.subtitle = element_text(size = 12, hjust = 0.5),
  axis.title = element_text(size = 12),
  axis.text.x = element_text(size = 11),
  axis.text.y = element_text(size = 11)
)

