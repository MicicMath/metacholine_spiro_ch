source("./data_prep.R")

dat = final_ast_hc_cf_dat

#################
# subset to group
#################

# dat = dat[dat$diagnosis_simple == "AST", , drop = FALSE]

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

#############################################
# permutation test between asthma and healthy
#############################################

M_ast = cor(
  dat[dat$diagnosis_simple == "AST", sensors],
  use = "pairwise.complete.obs"
)

M_hc = cor(
  dat[dat$diagnosis_simple == "HC", sensors],
  use = "pairwise.complete.obs"
)

# observed global difference between correlation matrices
idx = upper.tri(M_ast)

T_obs = sum((M_ast[idx] - M_hc[idx])^2)

# permutation test
set.seed(101)

B = 10000

diagnosis = dat$diagnosis_simple

T_perm = numeric(B)

for (b in 1:B) {

  diagnosis_perm = sample(diagnosis)

  M_ast_perm = cor(
    dat[diagnosis_perm == "AST", sensors],
    use = "pairwise.complete.obs"
  )

  M_hc_perm = cor(
    dat[diagnosis_perm == "HC", sensors],
    use = "pairwise.complete.obs"
  )

  T_perm[b] = sum(
    (M_ast_perm[idx] - M_hc_perm[idx])^2
  )
}

# permutation p-value
p_value = (sum(T_perm >= T_obs) + 1) / (B + 1)

T_obs
p_value
