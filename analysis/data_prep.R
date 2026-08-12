source("./setup.R")

########
# setup
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

cols = c(
  "patient_id",
  "m_id",
  "diagnosis",
  "diagnosis_status",
  sensors,
  "visit_exhal",
  "metacholine_test_result"
)

#############
# import rdpc
#############

rdcp_dat = read.csv("../dat/rdcp_exhal_analysis.csv")
colnames(rdcp_dat)[colnames(rdcp_dat) == "measurement_id_exhal"] = "m_id"

##############
# import e_dat
##############

e_dat_0 = read.csv("../dat/dat_ch_enose_20260105.csv")
colnames(e_dat_0)[colnames(e_dat_0) == "id"] = "m_id"

table(e_dat_0$m_id %in% rdcp_dat$m_id)

e_dat_0 = e_dat_0[e_dat_0$m_id %in% rdcp_dat$m_id, ]

e_dat_1 = read.csv("../dat/dat_ch_enose_20260808.csv")
colnames(e_dat_1)[colnames(e_dat_1) == "id"] = "m_id"

table(e_dat_1$m_id %in% rdcp_dat$m_id)

e_dat_1 = e_dat_1[e_dat_1$m_id %in% rdcp_dat$m_id, ]

e_dat = rbind(e_dat_0, e_dat_1)

##############################
# subset to asthma and healthy
##############################

ast_hc_mask = (rdcp_dat$diagnosis %in% c("HC", "PUL_ASTHMA_ALLERGIC", "PUL_ASTHMA_NONALLERGIC"))

rdcp_dat_ast_hc = rdcp_dat[ast_hc_mask, , drop = FALSE]

###################
# combine with rdcp
###################

final_ast_hc_dat = merge(e_dat, rdcp_dat_ast_hc, by = "m_id", sort = FALSE)

# write.csv(final_ast_hc_dat, "../dat/merged_e_dat_rdcp.csv", row.names = FALSE)

