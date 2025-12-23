#* 10: Generate Manuscript Tables
#+ 6.0a: Examine Normality
#! Did manual examination of this via below TernTables runs
#! T1 has 12/16 normality fail, T2 has 2/3, and T3 has 7/10
#! But, just going to use consider normality = TRUE to do test dynamically
#+ 10.0b: Set Conflicts
conflicts_prefer(purrr::compose)
#+ 10.1: Table 1 (Recipient Preoperative Characteristics) 
ST1 <- ternG(
  data = T1_data,
  exclude_vars = "Patient",
  group_var = "severe_PGD",
  descriptive = TRUE,
  output_docx = "Outputs/Tables/T1.docx",
  consider_normality = TRUE,
  show_test = FALSE,
  round_intg = FALSE,
  insert_subheads = TRUE
)
#+ 10.2: Table 2 (Donor Characteristics) 
ST2 <- ternG(
  data = T2_data,
  vars = NULL,
  exclude_vars = "Patient",
  group_var = "severe_PGD",
  descriptive = TRUE,
  output_docx = "Outputs/Tables/T2.docx",
  consider_normality = TRUE,
  show_test = FALSE,
  round_intg = TRUE,
  insert_subheads = TRUE
)
#+ 10.3: Table 3 (Table 3. Procurement/Surgical Factors and Perioperative/Post-Transplant Outcomes) 
ST3 <- ternG(
  data = T3_data,
  vars = NULL,
  exclude_vars = "Patient",
  group_var = "severe_PGD",
  descriptive = TRUE,
  output_docx = "Outputs/Tables/T3.docx",
  consider_normality = TRUE,
  show_test = FALSE,
  round_intg = TRUE,
  insert_subheads = TRUE
)