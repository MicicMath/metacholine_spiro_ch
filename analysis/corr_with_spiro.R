source("./data_prep.R")

dat = final_ast_hc_cf_dat
dat = dat[dat$diagnosis != "CF", , drop = FALSE]

###############################################################
# restrict to mid and high metacholine response only for asthma
###############################################################

mid_high_metach_response_mask = dat$metacholine_response %in% c("Low", "Mid", "High")

dat = dat[mid_high_metach_response_mask, , drop = FALSE]

table(dat$diagnosis_simple)

###############
# scatter plots
###############

plts = lapply(
  sensors,
  function(x) {
    plt = ggplot(dat, aes(x = !!sym(x), y = fev1_z)) +
      geom_point(shape = 21, size = 3) +
      geom_smooth(method = "lm", formula= y ~ x) +
      theme_classic()
    return(plt)
  }
)

plot_grid(
  plotlist = plts,
  nrow = 4,
  align = "hv"
)

##########
# cor test
##########

cors = lapply(
  sensors,
  function(x) {
    cor_test = cor.test(dat[[x]], dat$fev1_z, method = "pearson", na.rm = TRUE)
    res = data.frame(
      sensor = x,
      est = round(cor_test$estimate, 2),
      p_val = round(cor_test$p.value, 3) 
    )
    rownames(res) = NULL
    res
  }
)
cors = do.call(rbind, cors)
cors
