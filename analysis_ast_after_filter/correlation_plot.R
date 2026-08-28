source("./data_prep.R")

dat = final_ast_hc_cf_dat

#################
# subset to group
#################

# dat = dat[dat$diagnosis_simple == "HC", , drop = FALSE]

#####
# cor
#####

M = cor(dat[, sensors], use = "pairwise.complete.obs")

hc = hclust(dist(1 - M), method = "average")
ord = hc$order

# M = M[ord, ord]

M_long = reshape2::melt(M)

M_long = M_long[
  as.numeric(M_long$Var1) < as.numeric(M_long$Var2),
]

plt_cor = ggplot(M_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white", linewidth = 1) +
  
  geom_text(
    aes(label = sprintf("%.2f", value)),
    size = 3,
    color = "black"
  ) +
  
  scale_fill_gradient2(
    low = "#BB4444",
    mid = "#FFFFFF",
    high = "#4477AA",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation"
  ) +
  
  coord_fixed() +
  scale_x_discrete(position = "top") +
  
  theme_minimal(base_size = 12) +
  theme(
    # legend.position = "bottom",
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = .5, vjust = 1),
    axis.text.y = element_text()
  )
plt_cor

########
# hclust
########

M = cor(dat[, sensors], use = "pairwise.complete.obs")
d = as.dist(1 - abs(M))
hc = hclust(d, method = "complete")

plot(hc)
