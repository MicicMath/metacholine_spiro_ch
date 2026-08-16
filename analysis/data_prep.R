source("./setup.R")

#############
# import rdpc
#############

rdcp_dat = read.csv("../dat/rdcp_exhal_analysis.csv")
colnames(rdcp_dat)[colnames(rdcp_dat) == "measurement_id_exhal"] = "m_id"

###########################
# subset to relevant labels
###########################

rdcp_exhal_dat = rdcp_dat[
  !is.na(rdcp_dat$m_id) &
    rdcp_dat$diagnosis %in% c(
      "HC",
      "PUL_ASTHMA_ALLERGIC",
      "PUL_ASTHMA_NONALLERGIC",
      "PUL_CF"
    ), ,
  drop = FALSE
]

###############################################
# import eNose dat and filter by measurement id
###############################################

e_dat_0 = read.csv("../dat/dat_ch_enose_20260105.csv")
e_dat_1 = read.csv("../dat/dat_ch_enose_20260808.csv")

colnames(e_dat_0)[colnames(e_dat_0) == "id"] = "m_id"
colnames(e_dat_1)[colnames(e_dat_1) == "id"] = "m_id"

e_dat_0 = e_dat_0[e_dat_0$m_id %in% rdcp_exhal_dat$m_id, ]
e_dat_1 = e_dat_1[e_dat_1$m_id %in% rdcp_exhal_dat$m_id, ]

########################
# merge with REDCap data
########################

e_dat = rbind(e_dat_0, e_dat_1)

final_ast_hc_cf_dat = merge(
  e_dat,
  rdcp_exhal_dat,
  by = "m_id",
  sort = FALSE
)

#########################
# subset to relevant cols
#########################

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

cols = c(
  "patient_id",
  "m_id",
  "diagnosis",
  "diagnosis_status",
  sensors,
  "visit_exhal",
  "metacholine_test_result"
)

final_ast_hc_cf_dat = final_ast_hc_cf_dat[, cols]

###################
# relabel diagnosis
###################

final_ast_hc_cf_dat$diagnosis = ifelse(
  final_ast_hc_cf_dat$diagnosis == "PUL_ASTHMA_ALLERGIC",
  "AST_AL",
  ifelse(final_ast_hc_cf_dat$diagnosis == "PUL_ASTHMA_NONALLERGIC",
    "AST_NOAL",
    ifelse(final_ast_hc_cf_dat$diagnosis == "PUL_CF",
      "CF", "HC"
    )
  )
)

##########################################
# relabel metacholine provocation response
##########################################

final_ast_hc_cf_dat$metacholine_response = ifelse(
  (grepl("AST", final_ast_hc_cf_dat$diagnosis)) & (final_ast_hc_cf_dat$metacholine_test_result == 0),
  "Neg", ifelse(
    (grepl("AST", final_ast_hc_cf_dat$diagnosis)) & (final_ast_hc_cf_dat$metacholine_test_result == 1),
    "Low", ifelse(
      (grepl("AST", final_ast_hc_cf_dat$diagnosis)) & (final_ast_hc_cf_dat$metacholine_test_result == 2),
      "Mid", ifelse(
        (grepl("AST", final_ast_hc_cf_dat$diagnosis)) & (final_ast_hc_cf_dat$metacholine_test_result == 3),
        "High", "No_test"
      )
    )
  )
)

final_ast_hc_cf_dat$metacholine_response = ifelse(
  is.na(final_ast_hc_cf_dat$metacholine_response),
  "No_test",
  final_ast_hc_cf_dat$metacholine_response
)

final_ast_hc_cf_dat$metacholine_response = factor(
  final_ast_hc_cf_dat$metacholine_response,
  levels = c("No_test", "Neg", "Low", "Mid", "High")
)

###############################################
# diagnosis confirmed filter and only 1st visit
###############################################

# diag_confirmed_mask = !is.na(final_ast_hc_cf_dat$diagnosis_status) & (final_ast_hc_cf_dat$diagnosis_status %in% c(0, 1))
diag_confirmed_mask = !is.na(final_ast_hc_cf_dat$diagnosis_status) & final_ast_hc_cf_dat$diagnosis_status == 1

hc_mask = final_ast_hc_cf_dat$diagnosis == "HC"

first_visit_mask = !is.na(final_ast_hc_cf_dat$visit_exhal) &
  final_ast_hc_cf_dat$visit_exhal == 1

final_ast_hc_cf_dat = final_ast_hc_cf_dat[
  (diag_confirmed_mask | hc_mask) & first_visit_mask,
  ,
  drop = FALSE
]

##########################################################
# remove non allergic asthma and relabel to simpler labels
##########################################################

final_ast_hc_cf_dat = final_ast_hc_cf_dat[!(final_ast_hc_cf_dat$diagnosis == "AST_NOAL"), , drop = FALSE]

final_ast_hc_cf_dat$diagnosis_simple = ifelse(
  final_ast_hc_cf_dat$diagnosis == "AST_AL",
  "AST",
  final_ast_hc_cf_dat$diagnosis
)

######
# test
######

temp = final_ast_hc_cf_dat[final_ast_hc_cf_dat$diagnosis_status == 0, ]
temp = temp[temp$metacholine_response %in% c("Low", "Mid", "High"), ]
