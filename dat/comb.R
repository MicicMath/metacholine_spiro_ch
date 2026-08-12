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
  "dplyr",
  "MASS"
)

########
# import
########

sensors = c(
  "S1",
  "S3",
  "S4",
  "S5",
  "S6",
  "S7",
  "S1BH",
  "S2BH",
  "S3BH",
  "S4BH",
  "S5BH",
  "S6BH",
  "S7BH"
)

#############
# import rdpc
#############

rdcp_dat = read.csv("./rdcp_exhal_analysis.csv")
colnames(rdcp_dat)[colnames(rdcp_dat) == "measurement_id_exhal"] = "m_id"

############
# old vs old
############

e_dat_prev_0 = read.csv("./cad_and_ch_combined_cleaned.csv")
e_dat_prev_0 = e_dat_prev_0[e_dat_prev_0$Country == "CH_", ]
colnames(e_dat_prev_0)[colnames(e_dat_prev_0) == "Measurement_ID"] = "m_id"
e_dat_prev_0 = e_dat_prev_0[, c("m_id", sensors)]

e_dat_prev_1 = read.csv("./dat_ch_enose_20260105.csv")

e_dat_prev_1 = e_dat_prev_1[e_dat_prev_1$m_id %in% e_dat_prev_0$m_id, ]

e_dat_prev_1 = e_dat_prev_1[match(e_dat_prev_0$m_id, e_dat_prev_1$m_id), ]

table(e_dat_prev_1$m_id %in% rdcp_dat$m_id)

#####
# new
#####

e_dat_new = read.csv("./dat_ch_enose_20260808.csv")
colnames(e_dat_new)[colnames(e_dat_new) == "id"] = "m_id"

table(e_dat_new$m_id %in% rdcp_dat$m_id)

e_dat_new = e_dat_new[e_dat_new$m_id %in% rdcp_dat$m_id, ]

###################
# combine with rdcp
###################

cf_hc_mask_rdcp = (rdcp_dat$diagnosis %in% c("HC", "PUL_CF")) & (rdcp_dat$visit_exhal == 1) & ((rdcp_dat$diagnosis_status == 1) | is.na(rdcp_dat$diagnosis_status))

rdcp_dat_cf_hc = rdcp_dat[cf_hc_mask_rdcp, ]

e_dat = rbind(e_dat_prev_1, e_dat_new)

final_cf_hc_dat = merge(e_dat, rdcp_dat_cf_hc, by = "m_id", sort = FALSE)

# write.csv(final_cf_hc_dat, "../analysis_combined_new_data/tables/rdcp_cf_hc_swiss_meta_data.csv", row.names = FALSE)

######################################################
# transform back to old structure and combine with cad
######################################################

cad_and_ch = read.csv("./cad_and_ch_combined_cleaned.csv")
dat_cad = cad_and_ch[cad_and_ch$Country == "CAD_", ]
colnames(dat_cad)[colnames(dat_cad) == "Measurement_ID"] = "m_id"

cols_generic = c(
  "m_id",
  "label",
  "Country",
  "diagnosis",
  sensors
)

dat_cad = dat_cad[, cols_generic]

cols_ch = c(
  "m_id",
  "diagnosis",
  sensors
)

dat_ch = final_cf_hc_dat[, cols_ch]

dat_ch$Country = "CH_"
dat_ch$diagnosis = ifelse(dat_ch$diagnosis == "PUL_CF", "CF", "HC")
dat_ch$label = paste0(dat_ch$Country, dat_ch$diagnosis)

dat_ch = dat_ch[, cols_generic]

dat_tot = rbind(dat_cad, dat_ch)

# write.csv(dat_tot, "../dat/cad_and_ch_combined_cleaned_reduced_and_extended_20260808.csv", row.names = FALSE)
