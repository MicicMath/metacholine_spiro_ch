# Data checks: measurement dates and possible duplicate participant

Two issues were identified while reviewing the asthma-versus-healthy-control breath analysis. The findings below need checking against the original records.

## 1. Five healthy-control measurements have suspicious dates

These measurements belong to earlier eNose batches, but their recorded dates fall after the dates in the batch filenames. The suspicious dates are already present in both REDCap CSV exports; they were not introduced by sorting the trace plot or joining the sensor data.

Neighbouring measurement records suggest the following corrections:

| Measurement ID | Recorded date | Suggested date — to confirm |
|---|---|---|
| 130456 | 2026-10-29 | 2025-10-29 |
| 130494 | 2026-11-05 | 2025-11-05 |
| 130501 | 2026-11-06 | 2025-11-06 |
| 130558 | 2026-11-13 | 2025-11-13 |
| 131505 | 2026-12-06 | 2026-06-12 |

The first four belong to batch 20260105 and appear to have an incorrect year. Measurement 131505 belongs to batch 20260808 and appears to have its day and month reversed: the surrounding measurements 131504, 131506 and 131507 are dated 12 June 2026.

**Please verify these five dates against the acquisition records before correcting them.** The suggested dates are inferred, not confirmed. The current training/test split is based on eNose batch membership, not these dates.

## 2. Two test measurements may belong to the same child

| Measurement ID | Recorded date and time | Recorded visit |
|---|---|---|
| 131689 | 2026-07-08 14:50 | 1 |
| 131824 | 2026-08-05 15:00 | 1 |

These records share the same patient ID, date of birth and sex, but have different REDCap record IDs. Both are asthma measurements in test batch 20261027.

The analysis keeps records labelled “visit 1.” Because both records carry that label, both pass the filter and enter the test set. The filter does not independently identify the earliest measurement for each child.

**Please confirm whether these records represent the same child.** If they do, the child currently contributes two measurements to an analysis intended to use one first visit per child. The duplicate registration and/or visit labels need correcting. If they represent different children, the shared patient ID needs checking.

## Source references

- Current clinical data: `dat/rdcp_measurements_merged.csv`
- Older clinical export containing the same five suspicious dates: `dat/rdcp_exhal_analysis.csv`
- Data preparation and first-visit filtering: `analysis_ast_after_filter/data_prep.R`
- Chronological plots: `analysis_ast_after_filter/traces.R`
