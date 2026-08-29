source("./data_prep.R")

dat = final_ast_hc_cf_dat

dat = dat[order(dat$date_time_exhal), ]
dat$record = 1:nrow(dat)
dat$record_str = paste0("p", dat$record)

traces = lapply(
  sensors,
  function(x) {
    trace = ggplot(
      dat,
      aes(
        # x = reorder(record_str, record),
        x = date_time_exhal,
        y = !!sym(x),
        color = diagnosis_simple,
        fill = diagnosis_simple
      )
    ) +
      geom_point(size = 1, stroke = 1) +
      # geom_text_repel(aes(label = sample_base_id)) +
      # scale_color_viridis_d(option = "D", end = 0.9) +
      # scale_fill_viridis_d(option = "D", end = 0.9) +
      scale_shape_manual(values = c(1, 3)) +
      scale_color_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2") +
      labs(
        title = x,
        x = "",
        y = "Intensity [a.u.]"
      ) +
      theme_bw() +
      # facet_wrap(~ meta_label, nrow = 2, scales = "fixed") +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank()
      )
  }
)
names(traces) = sensors

plot_grid(
  plotlist = traces,
  ncol = 4,
  align = "hv"
)

################################################
# regression over idagnosis and measurement days
################################################

# measurement time in days since first measurement
dat$time_days = as.numeric(
  difftime(
    dat$date_time_exhal,
    min(dat$date_time_exhal, na.rm = TRUE),
    units = "days"
  )
)

dat$diagnosis_simple = relevel(
  factor(dat$diagnosis_simple),
  ref = "HC"
)

res_date = lapply(
  sensors,
  function(s) {

    fit = lm(
      dat[[s]] ~ diagnosis_simple + time_days + factor(operator_exhal),
      data = dat
    )

    sm = summary(fit)$coefficients

    data.frame(
      sensor = s,
      beta_diagnosis = sm["diagnosis_simpleAST", "Estimate"],
      p_diagnosis = sm["diagnosis_simpleAST", "Pr(>|t|)"],
      beta_time = sm["time_days", "Estimate"],
      p_time = sm["time_days", "Pr(>|t|)"]
    )
  }
)

res_date = do.call(rbind, res_date)

res_date$p_diagnosis_FDR = p.adjust(
  res_date$p_diagnosis,
  method = "BH"
)

res_date$p_time_FDR = p.adjust(
  res_date$p_time,
  method = "BH"
)

res_date


dat = dat[
  dat$operator_exhal %in% c("BHI", "OTH", "SSC"),
]

res_date = lapply(
  sensors,
  function(s) {

    fit = lm(
      dat[[s]] ~ diagnosis_simple + time_days + factor(operator_exhal),
      data = dat
    )

    sm = summary(fit)$coefficients
    dr = drop1(fit, test = "F")

    data.frame(
      sensor = s,
      beta_diagnosis = sm["diagnosis_simpleAST", "Estimate"],
      p_diagnosis = sm["diagnosis_simpleAST", "Pr(>|t|)"],
      beta_time = sm["time_days", "Estimate"],
      p_time = sm["time_days", "Pr(>|t|)"],
      p_operator = dr["factor(operator_exhal)", "Pr(>F)"]
    )
  }
)

res_date = do.call(rbind, res_date)

res_date$p_diagnosis_FDR = p.adjust(
  res_date$p_diagnosis,
  method = "BH"
)

res_date$p_time_FDR = p.adjust(
  res_date$p_time,
  method = "BH"
)

res_date$p_operator_FDR = p.adjust(
  res_date$p_operator,
  method = "BH"
)

res_date
