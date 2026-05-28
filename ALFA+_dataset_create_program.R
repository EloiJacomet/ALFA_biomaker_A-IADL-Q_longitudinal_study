#
# ============================================================================= #
# ---- ALFA+ LONGITUDINAL DATASET BUILDING SCRIPT ----
# ============================================================================= #
# Purpose:
#   - Build a clean longitudinal ALFA+ dataset from the ALFA+ cohort by integrating:
#     - Functional data
#     - Cognitive data
#     - Genetic data
#     - Fluid Biomarkers data
# 
# Description:
# This script performs data ingestion, cleaning, harmonization and merging of
# multiple modules (IADL, EDUCATION, ATN, BIOMARKERS, APOE, SUCOG, PACC,
# FRAILTY, etc.) into a unified longitudinal dataset at the participant-visit 
# level.
#
# The Pipeline includes:
#   - Conversion of the general dataset to long format
#   - Module-specific preprocessing and variable selection
#   - Quality control checks (missing values, duplicate keys, join validations...)
#   - Merging of modules
#   - Construction of baseline and longitudinal variables
#
# Inputs:
#   - ALFA+ CSV data release files (organized in module-specific folders)
#   - Biomarker data file (CSF and plasma)
#   - Genetic data (APOE)
#   - Frailty index data
#
# Outputs:
#   - Longitudinal dataset (CSV and RData file format)
#   - Optional metadata and QC outputs (if enabled)
#
# Reproducibility:
#   - Developed in R version 4.5.2
#   - Key Packages:
#     - dplyr (1.1.4)
#     - tidyr (1.3.1)
#     - readr (2.1.6)
#     - lubridate (1.9.4)
#     - janitor (2.2.1)
#
# Code: Eloi Jacomet & Federica Anastasi
# Date: April 6th 2026 
# NOTE: Dataset creation for ALFA+ IADL Study
# ============================================================================= #

rm(list = ls(all.names = TRUE))
invisible(gc())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(lubridate)
  library(janitor)
  library(tibble)
  library(readxl)
  library(here)
})

# ============================================================================= #
# 1. CONFIGURATION (EDIT ONLY THIS) ----
# ============================================================================= #

## ---- 1.1 Base directory ----

# Set working directory (try not to hardcode if possible)
setwd(here())
wd <- getwd()

## ---- 1.2 Release folder/date ----

release_folder    <- file.path(wd)
release_docs_date <- "20260406"

## ---- 1.3 Folders ----

dir_general   <- file.path(release_folder, "csv")
dir_docs      <- file.path(release_folder, "docs")

## ---- 1.4 Adjust external modules folder structure ----

# They may contain subfolders
dir_atn   <- file.path(release_folder, "ATN_status")
dir_biomk <- file.path(release_folder, "Biomarkers_Alfa_Cohort_ALL")
dir_apoe  <- file.path(release_folder, "APOE_Genetic")
dir_sucog   <- file.path(release_folder, "SUCOG")
dir_pacc   <- file.path(release_folder, "PACC")
dir_frailty   <- file.path(release_folder, "Frailty_Index")

## ---- 1.5 Modules to use ----

use_iadl      <- TRUE
use_education <- TRUE
use_atn       <- dir.exists(dir_atn)
use_biomk     <- dir.exists(dir_biomk)
use_apoe      <- dir.exists(dir_apoe)
use_sucog      <- dir.exists(dir_sucog)
use_pacc      <- dir.exists(dir_pacc)
use_frailty      <- dir.exists(dir_frailty)


## ---- 1.6 Detect prefix automatically from CSV names ----

files_csv <- list.files(dir_general, pattern="\\.csv(\\.csv)?$", full.names = FALSE)
stopifnot(length(files_csv) > 0)

probe_general <- files_csv[grepl("AlfaPlus_general_", files_csv)]
probe <- if (length(probe_general)) probe_general[1] else files_csv[1]
prefix <- sub("(AlfaPlus_general_|\\d{4}_).*", "", probe)  # everything before form id or general tag

prefix_rx <- paste0("^", gsub("\\+", "\\\\+", prefix)) # escape + for regex

## ---- 1.6 File name patterns (allow .csv or .csv.csv) ----

pattern_general <- paste0(prefix_rx, "AlfaPlus_general_.*\\.csv(\\.csv)?$")
pattern_iadl    <- paste0(prefix_rx, "1051_ADL_.*\\.csv(\\.csv)?$")
pattern_edu     <- paste0(prefix_rx, "1020_DEMOG_\\d{14}\\.csv(\\.csv)?$")
pattern_atn     <- "Aplus_ATN_.*\\.csv(\\.csv)?$"
pattern_biomk   <- "Biomarkers_Alfa_Cohort_.*\\.csv(\\.csv)?$"
pattern_apoe    <- "APOE_.*\\.csv(\\.csv)?$"
pattern_sucog     <- paste0(prefix_rx, "1023_SUCOG_\\d{14}\\.csv(\\.csv)?$")
pattern_pacc     <- "alfaplus_pacc_ID_Visit.csv(\\.csv)?$"
pattern_frailty     <- "ALFA_clean.csv(\\.csv)?$"


## ---- 1.7 Other options ----

baseline_visit   <- "V1"
decimal_mark_csv <- "."
join_keys        <- c("IdParticipante", "Visita")

## ---- 1.8. Biomarkers selection ----

biomk_assays_keep <- c("Abeta1-40","Abeta1-42","NFL","GFAP","pTau181","pTau217", "pTau231","tTau")
biomk_ntk_all <- FALSE     # keep NTK CSF/PLASMA combos automatically
biomk_combo_keep <- c(
  "Abeta1.40_CSF_Elecsys_UGOT_Roche",
  "Abeta1.40_PLASMA_Elecsys_UGOT_Roche",
  "Abeta1.42_CSF_Elecsys_UGOT_Roche",
  "Abeta1.42_PLASMA_Elecsys_UGOT_Roche",
  "GFAP_CSF_Elecsys_UGOT_Roche",
  "GFAP_PLASMA_Elecsys_UGOT_Roche",
  "NFL_CSF_Elecsys_UGOT_Roche",
  "NFL_PLASMA_Elecsys_UGOT_Roche",
  "pTau181_CSF_Elecsys_UGOT_Roche",
  "pTau181_PLASMA_Elecsys_UGOT_Roche",
  "tTau_CSF_Elecsys_UGOT_Roche",
  "pTau231_CSF_ELISA_UGOT_ADx.NeuroSciences",
  "pTau217_CSF_MSD_LILLY_in.house"
)  # only used if biomk_ntk_all == FALSE

# Check final configuration summary
cat("\nCONFIG SUMMARY\n",
    "wd: ", wd, "\n",
    "release_folder: ", release_folder, "\n",
    "release_docs_date: ", release_docs_date, "\n",
    "prefix: ", prefix, "\n",
    "dir_general: ", dir_general, "\n",
    "use_atn: ", use_atn, "  use_biomk: ", use_biomk, "  use_apoe: ", use_apoe,"  use_sucog: ", use_sucog,
    "  use_pacc: ", use_pacc,"  use_frailty: ", use_frailty, "\n", 
    sep = "")

stopifnot(dir.exists(dir_general))
stopifnot(dir.exists(dir_docs))



# ============================================================================= #
#  ---- 2. GENERIC HELPERS ----
# ============================================================================= #

## ---- 2.1. Read .csv files ----

read_alfaplus_csv <- function(path) {
  read_delim(
    path,
    delim = ";",
    show_col_types = FALSE,
    trim_ws = TRUE,
    na = c("", "NA"),
    locale = locale(decimal_mark = decimal_mark_csv)
  )
}

## ---- 2.2. Get files from directory ----

get_latest_file <- function(dir, pattern, recursive = FALSE) {
  stopifnot(dir.exists(dir))
  files <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = recursive)
  if (!length(files)) stop("No files in '", dir, "' with pattern '", pattern, "'.", call.=FALSE)
  files[which.max(file.mtime(files))]
}

## ---- 2.3. Simple and efficient Quality Control ----

qc_checkpoint <- function(df, label, keys = c("IdParticipante","Visita"),
                          must_have = NULL, show = 5) {
  cat("\n================ QC:", label, "================\n")
  cat("Rows:", nrow(df), "Cols:", ncol(df), "\n")
  cat("Columns:\n"); print(names(df))
  cat("Head:\n"); print(utils::head(df, show))
  
  if (!is.null(must_have)) {
    miss <- setdiff(must_have, names(df))
    if (length(miss)) stop(label, " missing columns: ", paste(miss, collapse=", "), call.=FALSE)
  }
  
  if (!is.null(keys) && all(keys %in% names(df))) {
    na_keys  <- sum(!complete.cases(df[, keys]))
    dup_keys <- sum(duplicated(df[, keys]))
    cat("Key NA rows:", na_keys, "\n")
    cat("Duplicate keys:", dup_keys, "\n")
    if (na_keys > 0) stop(label, " has NA in join keys.", call.=FALSE)
    if (dup_keys > 0) warning(label, " has duplicated keys (join may inflate rows).")
  }
  
  invisible(df)
}

## ---- 2.4. Join summary before merging ----

qc_join_plan <- function(master_keys_df, module_df, label, keys = c("IdParticipante","Visita")) {
  cat("\n------------- JOIN PLAN:", label, "-------------\n")
  stopifnot(all(keys %in% names(master_keys_df)))
  stopifnot(all(keys %in% names(module_df)))
  
  master_keys <- distinct(master_keys_df[, keys])
  module_keys <- distinct(module_df[, keys])
  matched <- inner_join(master_keys, module_keys, by = keys)
  
  cat("Master unique keys:", nrow(master_keys), "\n")
  cat("Module unique keys:", nrow(module_keys), "\n")
  cat("Matched keys:", nrow(matched), " (",
      round(100*nrow(matched)/nrow(master_keys), 1), "% of master)\n", sep="")
  
  dup_keys <- sum(duplicated(module_df[, keys]))
  if (dup_keys > 0) warning(label, " duplicated keys -> join will inflate rows unless you deduplicate.")
  invisible(TRUE)
}

## ---- 2.5. Safe merge of modules ----

safe_left_join <- function(x, y, by, label) {
  n_before <- nrow(x)
  out <- left_join(x, y, by = by)
  if (nrow(out) != n_before) {
    stop(label, ": join inflated rows (", n_before, " -> ", nrow(out),
         "). Module likely has duplicate keys.", call.=FALSE)
  }
  out
}

## ---- 2.6. Change delimiter for certain .csv files ----
convert_csv_comma_to_semicolon <- function(path_in, path_out = path_in) {
  first_line <- readLines(path_in, n = 1, encoding = "UTF-8")
  if (!grepl(",", first_line)) {
    message("\nThe file appears to already use ';'. Will not convert.\n")
    return(invisible(NULL))
  }
  df <- readr::read_csv(
    path_in,
    show_col_types = FALSE,
    na = c("", "NA")
  )
  readr::write_delim(
    df,
    path_out,
    delim = ";",
    na = ""
  )
  message("\nFile converted with delim ';':\n", path_out)
}

# ============================================================================= #
#  ---- 3. Create General dataset and convert to Long format ----
# ============================================================================= #

## ---- 3.1. Extract general file and quick quality control ----

file_general <- get_latest_file(dir_general, pattern_general)
cat("\nUsing GENERAL file:\n  ", file_general, "\n", sep="")
alfaplus_raw <- read_alfaplus_csv(file_general)
qc_checkpoint(alfaplus_raw, "GENERAL raw", keys = NULL)

## ---- 3.2. Convert general dataset to long format ----

make_alfaplus_long <- function(alfaplus_raw) {
  alfaplus_raw %>%
    pivot_longer(
      cols = starts_with("Date_"),
      names_to = "Visita_raw",
      values_to = "Fecha"
    ) %>%
    filter(!is.na(Fecha)) %>%
    mutate(
      Visita = str_extract(Visita_raw, "V\\d(?:_\\d)?") %>%
        str_replace("_", "_"),
      Fecha = suppressWarnings(dmy(Fecha))
    ) %>%
    filter(!is.na(Visita)) %>%
    select(IdParticipante, Visita, Fecha, everything())
}

alfaplus_long <- make_alfaplus_long(alfaplus_raw) %>%
  mutate(
    IdParticipante = suppressWarnings(as.numeric(IdParticipante)),
    Visita = as.factor(Visita)
  )

## ---- 3.3. General variable selection ----

vars_keep <- c(
  "IdParticipante", "Visita", "Fecha",
  "DOB", "Sex", "Type", "Study_state",
  "Age_V1", "Age_V2", "Age_V2_5", "Age_V3", "Age_V3_5", "Dx_dt", "Dx", "SCD_def"
)

alfaplus_long <- alfaplus_long %>%
  select(any_of(vars_keep)) %>%
  arrange(IdParticipante, Visita)

## ---- 3.4. Quick QC of general long dataset ----

qc_checkpoint(alfaplus_long, "GENERAL long", must_have = c("IdParticipante","Visita","Fecha"))
master_keys_df <- distinct(select(alfaplus_long, all_of(join_keys)))


# ============================================================================= #
#  ---- 4. MODULE HELPERS ----
# ============================================================================= #

## ---- 4.1 IADL ----

clean_iadl <- function(path_iadl) {
  adl_raw <- read_alfaplus_csv(path_iadl)
  adl_named <- adl_raw %>% janitor::row_to_names(row_number = 1)
  
  cols_keep <- c(
    "IdParticipante", "Visita", "FechaResultado",
    "glb_Form_nREASND", "AMSFormOL", "F02_6_007C", "AIADL_TraitScore"
  )
  
  adl_named %>%
    select(any_of(cols_keep)) %>%
    mutate(
      FechaResultado   = suppressWarnings(dmy(FechaResultado)),
      F02_6_007C       = suppressWarnings(parse_number(F02_6_007C)),
      AIADL_TraitScore = suppressWarnings(parse_number(AIADL_TraitScore)),
      glb_Form_nREASND = as.factor(glb_Form_nREASND),
      AMSFormOL        = suppressWarnings(parse_number(AMSFormOL)),
      IdParticipante   = suppressWarnings(as.numeric(IdParticipante)),
      Visita           = as.factor(Visita)
    )
}

## ---- 4.2 EDUCATION ----

clean_education <- function(path_edu) {
  edu_raw <- read_alfaplus_csv(path_edu)
  edu_named <- edu_raw %>% janitor::row_to_names(row_number = 1)
  
  edu_named %>%
    select(IdParticipante, Visita, YearsEducation = F01_010B
) %>%
    mutate(
      IdParticipante = suppressWarnings(as.numeric(IdParticipante)),
      Visita         = as.factor(Visita),
      YearsEducation = suppressWarnings(as.numeric(YearsEducation))
    )
}

## ---- 4.3 ATN ----

clean_atn <- function(atn_input, baseline_visit = "V1") {
  
  atn_df <- if (is.character(atn_input)) read_alfaplus_csv(atn_input) else atn_input
  
  # normalize keys
  atn_df <- atn_df %>%
    rename(
      IdParticipante = idparticipante,
      Visita         = Visit
    )
  
  # build visit-level + baseline variables
  out <- atn_df %>%
    transmute(
      IdParticipante = suppressWarnings(as.numeric(IdParticipante)),
      Visita         = as.factor(Visita),
      
      A_status = suppressWarnings(as.integer(ab42_ab40_rslt)),
      T_status = suppressWarnings(as.integer(pTau181_rslt)),
      N_status = suppressWarnings(as.integer(tTau_rslt)),  # because you actually HAVE it
      
      AT_status = dplyr::case_when(
        A_status == 0 & T_status == 0 ~ "A-T-",
        A_status == 0 & T_status == 1 ~ "A-T+",
        A_status == 1 & T_status == 0 ~ "A+T-",
        A_status == 1 & T_status == 1 ~ "A+T+",
        TRUE ~ NA_character_
      ),
      
      # keep numeric ATN code too (useful)
      ATN_code = suppressWarnings(as.integer(ATN_status)),
      
      # and make the label you want (this replaces your missing ATN_label)
      ATN_status = dplyr::recode(
        suppressWarnings(as.integer(ATN_status)),
        `1`="A-T-N-", `2`="A-T+N-", `3`="A-T-N+", `4`="A-T+N+",
        `5`="A+T-N-", `6`="A+T+N-", `7`="A+T-N+", `8`="A+T+N+",
        .default = NA_character_
      )
    ) %>%
    mutate(
      AT_status = factor(AT_status, levels = c("A-T-","A-T+","A+T-","A+T+")),
      ATN_status = factor(
        ATN_status,
        levels = c("A-T-N-","A-T+N-","A-T-N+","A-T+N+",
                   "A+T-N-","A+T+N-","A+T-N+","A+T+N+"),
        ordered = TRUE
      )
    ) %>%
    arrange(IdParticipante, Visita) %>%
    group_by(IdParticipante) %>%
    mutate(
      A_bl = {v <- A_status[Visita == baseline_visit & !is.na(A_status)];
      if (length(v)) v[1] else dplyr::first(A_status[!is.na(A_status)], default = NA_integer_)},
      T_bl = {v <- T_status[Visita == baseline_visit & !is.na(T_status)];
      if (length(v)) v[1] else dplyr::first(T_status[!is.na(T_status)], default = NA_integer_)},
      N_bl = {v <- N_status[Visita == baseline_visit & !is.na(N_status)];
      if (length(v)) v[1] else dplyr::first(N_status[!is.na(N_status)], default = NA_integer_)},
      
      AT_bl = dplyr::case_when(
        A_bl == 0 & T_bl == 0 ~ "A-T-",
        A_bl == 0 & T_bl == 1 ~ "A-T+",
        A_bl == 1 & T_bl == 0 ~ "A+T-",
        A_bl == 1 & T_bl == 1 ~ "A+T+",
        TRUE ~ NA_character_
      ),
      ATN_bl = {v <- ATN_status[Visita == baseline_visit & !is.na(ATN_status)];
      if (length(v)) v[1] else dplyr::first(ATN_status[!is.na(ATN_status)], default = NA)}
    ) %>%
    ungroup() %>%
    mutate(AT_bl = factor(AT_bl, levels = c("A-T-","A-T+","A+T-","A+T+")))
  
  out
}


## ---- 4.4 BIOMARKERS ----

clean_biomarkers <- function(
    path_biomk,
    assays_keep,
    combos_keep = NULL,
    assay_col = "Assay",
    id_col = "idparticipante",
    visit_col = "Visit"
) {
  biomk_raw <- read_alfaplus_csv(path_biomk)
  
  value_col <- "Result"
  if (!value_col %in% names(biomk_raw)) stop("BIOMK: 'Result' column not found.", call.=FALSE)
  
  if (!"IdParticipante" %in% names(biomk_raw)) biomk_raw <- biomk_raw %>% rename(IdParticipante = !!id_col)
  if (!visit_col %in% names(biomk_raw)) stop("BIOMK: Visit column missing.", call.=FALSE)
  biomk_raw <- biomk_raw %>% rename(Visita = !!visit_col)
  
  biomk_raw <- biomk_raw %>%
    mutate(
      IdParticipante = suppressWarnings(as.numeric(IdParticipante)),
      Visita = as.factor(Visita)
    )
  
  all_assays <- sort(unique(biomk_raw[[assay_col]]))
  cat("\n[BIOMK] Assays available:\n"); print(all_assays)
  
  if (!length(assays_keep)) stop("BIOMK: biomk_assays_keep is empty.", call.=FALSE)
  assays_keep2 <- intersect(assays_keep, all_assays)
  if (!length(assays_keep2)) stop("BIOMK: none of assays_keep found in file.", call.=FALSE)
  
  biomk_filt <- biomk_raw %>% filter(.data[[assay_col]] %in% assays_keep2)
  
  if (biomk_ntk_all) {
    biomk_filt <- biomk_filt %>%
      filter(SampleType %in% c("PLASMA","CSF"), grepl("NTK", Kit))
  }
  
  biomk_filt <- biomk_filt %>%
    mutate(
      combo_raw = paste(.data[[assay_col]], SampleType, Platform, LbNam, Company, sep = "_"),
      combo_clean = combo_raw %>% make.names() %>% str_replace_all("_+", "_")
    )
  
  combo_summary <- biomk_filt %>%
    group_by(combo_clean) %>%
    summarise(
      Assay = first(.data[[assay_col]]),
      SampleType = first(SampleType),
      Platform = first(Platform),
      LabName = first(LbNam),
      Company = first(Company),
      n_measures = n(),
      n_participants = n_distinct(IdParticipante),
      n_visits = n_distinct(Visita),
      .groups = "drop"
    )
  
  cat("\n[BIOMK] combos summary:\n")
  print(combo_summary, n = Inf, width = Inf)
  
  combos_use <- if (biomk_ntk_all) combo_summary$combo_clean else intersect(combos_keep, combo_summary$combo_clean)
  if (!length(combos_use)) stop("BIOMK: no combos selected (check filters).", call.=FALSE)
  
  biomk_sel <- biomk_filt %>% filter(combo_clean %in% combos_use)
  
  biomk_meta <- biomk_sel %>%
    select(IdParticipante, Visita, !!assay_col, SampleType, Platform, LbNam, Company, Kit,
           collection_dt, AnaDt, Result, ResultU, comment, combo_clean)
  
  biomk_wide <- biomk_sel %>%
    group_by(IdParticipante, Visita, combo_clean) %>%
    summarise(biomk_value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      id_cols = c(IdParticipante, Visita),
      names_from = combo_clean,
      values_from = biomk_value,
      names_prefix = "BIO_"
    )
  
  list(wide = biomk_wide, meta = biomk_meta, combo_summary = combo_summary, combos_used = combos_use)
}

## ---- 4.5 APOE ----

clean_apoe <- function(apoe_input) {
  
  apoe_df <- if (is.character(apoe_input) && length(apoe_input) == 1) {
    read_alfaplus_csv(apoe_input)
  } else apoe_input
  
  apoe_df %>%
    transmute(
      IdParticipante = suppressWarnings(as.numeric(IdParticipante)),
      APOE_genetic   = as.character(APOE_genetic),
      APOE_final     = suppressWarnings(as.integer(APOE_final)),
      APOE_class     = suppressWarnings(as.integer(APOE_class)),
      APOE_binary    = suppressWarnings(as.integer(APOE_binary))
    ) %>%
    arrange(IdParticipante) %>%
    group_by(IdParticipante) %>%
    summarise(
      APOE_genetic = dplyr::first(na.omit(APOE_genetic), default = NA_character_),
      APOE_final   = dplyr::first(na.omit(APOE_final),   default = NA_integer_),
      APOE_class   = dplyr::first(na.omit(APOE_class),   default = NA_integer_),
      APOE_binary  = dplyr::first(na.omit(APOE_binary),  default = NA_integer_),
      .groups = "drop"
    )
}

## ---- 4.6 SUCOG ----

clean_sucog <- function(path_sucog) {
  
  sucog_raw <- read_alfaplus_csv(path_sucog)
  sucog_named <- sucog_raw %>% janitor::row_to_names(row_number = 1)
  
  cols_keep <- c(
    "IdParticipante",
    "Visita",
    "FechaResultado",
    "sucog_decr2yrs",
    "ME_Total_SUCOG",
    "LE_Total_SUCOG",
    "FE_Total_SUCOG",
    "F02_3_036C"
  )
  
  sucog_named %>%
    select(IdParticipante, 
           Visita, 
           FechaResultado_SUCOG = FechaResultado, 
           sucog_decr2yrs,
           ME_Total_SUCOG, 
           LE_Total_SUCOG, 
           FE_Total_SUCOG, 
           Total_SUCOG = F02_3_036C) %>%
    mutate(
      FechaResultado_SUCOG = suppressWarnings(lubridate::dmy(FechaResultado_SUCOG)),
      sucog_decr2yrs = suppressWarnings(as.integer(sucog_decr2yrs)),
      ME_Total_SUCOG        = suppressWarnings(as.integer(ME_Total_SUCOG)),
      LE_Total_SUCOG        = suppressWarnings(as.integer(LE_Total_SUCOG)),
      FE_Total_SUCOG        = suppressWarnings(as.integer(FE_Total_SUCOG)),
      Total_SUCOG        = suppressWarnings(as.integer(Total_SUCOG)),
      IdParticipante = suppressWarnings(as.numeric(IdParticipante)),
      Visita         = as.factor(Visita)
    )
}

## ---- 4.7 PACC ----

# Change csv delimiter to ;
convert_csv_comma_to_semicolon(file.path(dir_pacc, "alfaplus_pacc_ID_Visit.csv"))
                               
clean_pacc <- function(path_pacc) {
  
  pacc_raw <- read_alfaplus_csv(path_pacc)
  pacc_named <- pacc_raw
  
  cols_keep <- c(
    "IdParticipante",
    "Visita",
    "pacc"
  )
  
  pacc_named %>%
    select(IdParticipante, Visita, PACC = pacc) %>%
    mutate(
      PACC        = suppressWarnings(as.numeric(PACC)),
      IdParticipante = suppressWarnings(as.numeric(IdParticipante)),
      Visita         = as.factor(Visita)
    )
}


## ---- 4.8 Frailty Index ----

# Change csv delimiter to ;
convert_csv_comma_to_semicolon(file.path(dir_frailty, "ALFA_clean.csv"))

clean_frailty <- function(path_frailty) {

  frailty_raw <- read_alfaplus_csv(path_frailty)
  frailty_named <- frailty_raw

  cols_keep <- c(
    "ParticipantID",
    "FI"
  )

  frailty_named %>%
    select(IdParticipante = ParticipantID, Frailty_Index = FI) %>%
    mutate(
      Frailty_Index        = suppressWarnings(as.numeric(Frailty_Index)),
      IdParticipante = suppressWarnings(as.numeric(IdParticipante)),
    )
}

# ============================================================================= #
#  ---- 5. LOAD + QC BEFORE MERGE ----
# ============================================================================= #

source_meta <- list()

## ---- 5.1 IADL ----

iadl_clean <- NULL
if (use_iadl) {
  file_iadl <- get_latest_file(dir_general, pattern_iadl)
  cat("\nUsing IADL file:\n  ", file_iadl, "\n", sep="")
  iadl_clean <- clean_iadl(file_iadl)
  
  qc_checkpoint(iadl_clean, "IADL", must_have = c("IdParticipante","Visita","AIADL_TraitScore"))
  qc_join_plan(master_keys_df, iadl_clean, "IADL", join_keys)
  
  source_meta[["IADL"]] <- tibble(component="IADL", file_path=file_iadl,
                                  file_mtime=file.info(file_iadl)$mtime, file_size=file.info(file_iadl)$size,
                                  n_rows_used=nrow(iadl_clean))
}

View(iadl_clean)

dup_keys <- iadl_clean %>%
  count(IdParticipante, Visita, name = "n") %>%
  filter(n > 1)
dup_keys

iadl_dup_rows <- iadl_clean %>%
  semi_join(dup_keys, by = c("IdParticipante", "Visita")) %>%
  arrange(IdParticipante, Visita, FechaResultado)
iadl_dup_rows

iadl_clean <- iadl_clean %>%
  arrange(IdParticipante, Visita, desc(FechaResultado)) %>%
  distinct(IdParticipante, Visita, .keep_all = TRUE)
sum(duplicated(iadl_clean[, c("IdParticipante","Visita")]))


## ---- 5.2 EDUCATION ----

education_clean <- NULL
if (use_education) {
  file_edu <- get_latest_file(dir_general, pattern_edu)
  cat("\nUsing EDUCATION file:\n  ", file_edu, "\n", sep="")
  education_clean <- clean_education(file_edu)
  
  qc_checkpoint(education_clean, "EDUCATION", must_have = c("IdParticipante","Visita","YearsEducation"))
  qc_join_plan(master_keys_df, education_clean, "EDUCATION", join_keys)
  
  source_meta[["EDUCATION"]] <- tibble(component="EDUCATION", file_path=file_edu,
                                       file_mtime=file.info(file_edu)$mtime, file_size=file.info(file_edu)$size,
                                       n_rows_used=nrow(education_clean))
}

dup_keys_edu <- education_clean %>%
  count(IdParticipante, Visita, name = "n") %>%
  filter(n > 1)

dup_keys_edu

edu_dup_rows <- education_clean %>%
  semi_join(dup_keys_edu, by = c("IdParticipante","Visita")) %>%
  arrange(IdParticipante, Visita)

edu_dup_rows
View(edu_dup_rows)

education_clean <- education_clean %>%
  distinct(IdParticipante, Visita, YearsEducation, .keep_all = TRUE) %>%
  distinct(IdParticipante, Visita, .keep_all = TRUE)

## ---- 5.3 ATN  ----

atn_clean <- NULL
if (use_atn) {
  file_atn <- get_latest_file(dir_atn, pattern_atn, recursive = TRUE)
  cat("\nUsing ATN file:\n  ", file_atn, "\n", sep="")
  
  # clean_atn can take a path OR a df (depending on how you defined it)
  # If your clean_atn expects a df, keep the 2 lines below.
  atn_raw  <- read_alfaplus_csv(file_atn)
  atn_clean <- clean_atn(atn_raw, baseline_visit = baseline_visit)
  
  # QC: must_have should match the OUTPUT columns of clean_atn()
  qc_checkpoint(
    atn_clean, "ATN",
    must_have = c("IdParticipante","Visita","A_status","T_status","N_status","ATN_status","ATN_bl")
  )
  
  cat("\nATN_status distribution:\n")
  print(table(atn_clean$ATN_status, useNA = "ifany"))
  
  cat("\nVisits in ATN:\n")
  print(table(atn_clean$Visita, useNA = "ifany"))
  
  qc_join_plan(master_keys_df, atn_clean, "ATN", join_keys)
  
  source_meta[["ATN"]] <- tibble(
    component  = "ATN",
    file_path  = file_atn,
    file_mtime = file.info(file_atn)$mtime,
    file_size  = file.info(file_atn)$size,
    n_rows_used = nrow(atn_clean)
  )
}

# --- baseline-only table: 1 row per participant (fills across all visits later)
atn_baseline <- atn_clean %>%
  distinct(IdParticipante, A_bl, T_bl, N_bl, AT_bl, ATN_bl)

stopifnot(!anyDuplicated(atn_baseline$IdParticipante))

## ---- 5.4 BIOMARKERS ----

biomk_clean <- NULL
if (use_biomk) {
  file_biomk <- get_latest_file(dir_biomk, pattern_biomk, recursive = TRUE)
  cat("\nUsing BIOMK file:\n  ", file_biomk, "\n", sep="")
  biomk_clean <- clean_biomarkers(
    path_biomk   = file_biomk,
    assays_keep  = biomk_assays_keep,
    combos_keep  = biomk_combo_keep
  )
  
  qc_checkpoint(biomk_clean$wide, "BIOMK wide", must_have = c("IdParticipante","Visita"))
  if (!any(grepl("^BIO_", names(biomk_clean$wide)))) {
    stop("BIOMK wide has no BIO_ columns. Check assays/filters.", call.=FALSE)
  }
  qc_join_plan(master_keys_df, biomk_clean$wide, "BIOMK wide", join_keys)
  
  source_meta[["BIOMK"]] <- tibble(component="BIOMK", file_path=file_biomk,
                                   file_mtime=file.info(file_biomk)$mtime, file_size=file.info(file_biomk)$size,
                                   n_rows_used=nrow(biomk_clean$wide))
}

## ---- 5.5 APOE ----

apoe_clean <- NULL
if (use_apoe) {
  # try recursive, APOE often sits under dated/CSV subfolder
  file_apoe <- tryCatch(
    get_latest_file(dir_apoe, pattern_apoe, recursive = TRUE),
    error = function(e) {
      cat("\n[APOE] No match with pattern_apoe. Listing files (first 50):\n")
      print(head(list.files(dir_apoe, recursive = TRUE), 50))
      stop(e$message, call.=FALSE)
    }
  )
  
  cat("\nUsing APOE file:\n  ", file_apoe, "\n", sep="")
  apoe_clean <- clean_apoe(file_apoe)
  
  # If Visita exists, treat like visit-level; if not, participant-level only
  if ("Visita" %in% names(apoe_clean)) {
    qc_checkpoint(apoe_clean, "APOE", must_have = c("IdParticipante","Visita"))
    qc_join_plan(master_keys_df, apoe_clean, "APOE", join_keys)
  } else {
    qc_checkpoint(apoe_clean, "APOE (participant-level)", keys = c("IdParticipante"),
                  must_have = c("IdParticipante"))
    # join plan: match by IdParticipante only
    cat("\n[JOIN PLAN] APOE participant-level\n")
    cat("Master unique participants:",
        n_distinct(master_keys_df$IdParticipante), "\n")
    cat("APOE participants:",
        n_distinct(apoe_clean$IdParticipante), "\n")
  }
  
  source_meta[["APOE"]] <- tibble(component="APOE", file_path=file_apoe,
                                  file_mtime=file.info(file_apoe)$mtime, file_size=file.info(file_apoe)$size,
                                  n_rows_used=nrow(apoe_clean))
}

## ---- 5.6 SUCOG ----

sucog_clean <- NULL
if (use_sucog) {
  file_sucog <- get_latest_file(dir_sucog, pattern_sucog)
  cat("\nUsing SUCOG file:\n  ", file_sucog, "\n", sep = "")
  
  sucog_clean <- clean_sucog(file_sucog)
  
  qc_checkpoint(
    sucog_clean,
    "SUCOG",
    must_have = c("IdParticipante", "Visita", "sucog_decr2yrs","ME_Total_SUCOG", "LE_Total_SUCOG", "FE_Total_SUCOG", "Total_SUCOG")
  )
  
  qc_join_plan(master_keys_df, sucog_clean, "SUCOG", join_keys)
  
  source_meta[["SUCOG"]] <- tibble(
    component   = "SUCOG",
    file_path   = file_sucog,
    file_mtime  = file.info(file_sucog)$mtime,
    file_size   = file.info(file_sucog)$size,
    n_rows_used = nrow(sucog_clean)
  )
}

dup_keys_sucog <- sucog_clean %>%
  count(IdParticipante, Visita, name = "n") %>%
  filter(n > 1)

dup_keys_sucog

sucog_dup_rows <- sucog_clean %>%
  semi_join(dup_keys_sucog, by = c("IdParticipante", "Visita")) %>%
  arrange(IdParticipante, Visita)

sucog_dup_rows
View(sucog_dup_rows)

sucog_clean <- sucog_clean %>%
  arrange(IdParticipante, Visita) %>%
  distinct(IdParticipante, Visita, .keep_all = TRUE)

## ---- 5.7 PACC ----

pacc_clean <- NULL
if (use_pacc) {
  file_pacc <- get_latest_file(dir_pacc, pattern_pacc)
  cat("\nUsing PACC file:\n  ", file_pacc, "\n", sep = "")
  
  pacc_clean <- clean_pacc(file_pacc)
  
  qc_checkpoint(
    pacc_clean,
    "PACC",
    must_have = c("IdParticipante", "Visita", "PACC")
  )
  
  qc_join_plan(master_keys_df, pacc_clean, "PACC", join_keys)
  
  source_meta[["PACC"]] <- tibble(
    component   = "PACC",
    file_path   = file_pacc,
    file_mtime  = file.info(file_pacc)$mtime,
    file_size   = file.info(file_pacc)$size,
    n_rows_used = nrow(pacc_clean)
  )
}

dup_keys_pacc <- pacc_clean %>%
  count(IdParticipante, Visita, name = "n") %>%
  filter(n > 1)

dup_keys_pacc

pacc_dup_rows <- pacc_clean %>%
  semi_join(dup_keys_pacc, by = c("IdParticipante", "Visita")) %>%
  arrange(IdParticipante, Visita)

pacc_dup_rows
View(pacc_dup_rows)

pacc_clean <- pacc_clean %>%
  arrange(IdParticipante, Visita) %>%
  distinct(IdParticipante, Visita, .keep_all = TRUE)

## ---- 5.8 Frailty Index ----

frailty_clean <- NULL
if (use_frailty) {
  file_frailty <- get_latest_file(dir_frailty, pattern_frailty)
  cat("\nUsing PACC file:\n  ", file_frailty, "\n", sep = "")

  frailty_clean <- clean_frailty(file_frailty)

  qc_checkpoint(
    frailty_clean,
    "FRAILTY",
    must_have = c("IdParticipante", "Frailty_Index")
  )

  source_meta[["FRAILTY"]] <- tibble(
    component   = "FRAILTY",
    file_path   = file_frailty,
    file_mtime  = file.info(file_frailty)$mtime,
    file_size   = file.info(file_frailty)$size,
    n_rows_used = nrow(frailty_clean)
  )
}

dup_keys_frailty <- frailty_clean %>%
  count(IdParticipante, name = "n") %>%
  filter(n > 1)

dup_keys_frailty

frailty_dup_rows <- frailty_clean %>%
  semi_join(dup_keys_frailty, by = "IdParticipante") %>%
  arrange(IdParticipante)

frailty_dup_rows
View(frailty_dup_rows)

frailty_clean <- frailty_clean %>%
  arrange(IdParticipante) %>%
  distinct(IdParticipante, .keep_all = TRUE)

# ============================================================================= #
#  ---- 6. BUILD DATASET (MERGE ONLY AFTER QC) ----
# ============================================================================= #

alfa_long <- alfaplus_long

## --- 6.1 IADL ----
if (use_iadl && !is.null(iadl_clean)) {
  alfa_long <- safe_left_join(alfa_long, iadl_clean, by = join_keys, label = "IADL")
}

## ---- 6.2 EDUCATION ----
if (use_education && !is.null(education_clean)) {
  alfa_long <- safe_left_join(alfa_long, education_clean, by = join_keys, label = "EDUCATION")
}

## ---- 6.3 ATN visit-level (by Id + Visita): statuses that vary with visit (NA if missing that visit) ----
if (use_atn && !is.null(atn_clean)) {
  atn_visit <- atn_clean %>%
    select(IdParticipante, Visita, A_status, T_status, N_status, AT_status, ATN_code, ATN_status)
  
  alfa_long <- safe_left_join(alfa_long, atn_visit, by = join_keys, label = "ATN visit")
}

## ---- 6.4 ATN baseline (by Id only): fills across ALL visits ----
if (use_atn && exists("atn_baseline") && !is.null(atn_baseline)) {
  n_before <- nrow(alfa_long)
  alfa_long <- left_join(alfa_long, atn_baseline, by = "IdParticipante")
  stopifnot(nrow(alfa_long) == n_before)
}

## ---- 6.5 BIOMARKERS ----
if (use_biomk && !is.null(biomk_clean)) {
  alfa_long <- safe_left_join(alfa_long, biomk_clean$wide, by = join_keys, label = "BIOMK wide")
}

## ---- 6.6 APOE ----
if (use_apoe && !is.null(apoe_clean)) {
  if ("Visita" %in% names(apoe_clean)) {
    alfa_long <- safe_left_join(alfa_long, apoe_clean, by = join_keys, label = "APOE")
  } else {
    n_before <- nrow(alfa_long)
    alfa_long <- left_join(alfa_long, apoe_clean, by = "IdParticipante")
    stopifnot(nrow(alfa_long) == n_before)
  }
}
## ---- 6.7 SUCOG ----
if (use_sucog && !is.null(sucog_clean)) {
  alfa_long <- safe_left_join(alfa_long, sucog_clean, by = join_keys, label = "SUCOG")
}

## ---- 6.8 PACC ----
if (use_pacc && !is.null(pacc_clean)) {
  alfa_long <- safe_left_join(alfa_long, pacc_clean, by = join_keys, label = "PACC")
}

## ---- 6.9 Frailty Index ----
if (use_frailty && !is.null(frailty_clean)) {
  alfa_long <- safe_left_join(alfa_long, frailty_clean, by = "IdParticipante", label = "FRAILTY_INDEX")
}


cat("\nFinal alfa_long:\n")
glimpse(alfa_long)
cat("\nN participants:", n_distinct(alfa_long$IdParticipante), "\n")
cat("Visits:\n"); print(table(alfa_long$Visita))


# ============================================================================= #
#  ---- 7. FINAL RENAMES / TYPES ----
# ============================================================================= #

alfa_long <- alfa_long %>%
  rename(
    alfa_id = IdParticipante,
    visit   = Visita,
    sex     = Sex,
    dob     = DOB,
    date    = Fecha,
    alfa_type = Type,
    iadl_date = glb_Form_nREASND,                       
    iadl_score = F02_6_007C,
    iadl_trait_score = AIADL_TraitScore,
    iadl_online = AMSFormOL,
    iadl_na_motiv = glb_Form_nREASND,
    diagnostic_date = Dx_dt,
    diagnostic = Dx
  ) %>%
  mutate(
    alfa_id = suppressWarnings(as.numeric(alfa_id)),
    visit   = as.factor(visit),
    sex     = as.factor(sex)
  )

# convert common date vars if character
date_like_vars <- c("date", "dob", "iadl_date")
for (v in date_like_vars) {
  if (v %in% names(alfa_long) && is.character(alfa_long[[v]])) {
    alfa_long[[v]] <- suppressWarnings(lubridate::dmy(alfa_long[[v]]))
  }
}

cat("\nSummary of classes:\n")
print(sapply(alfa_long, class))



# ============================================================================= #
#  ---- 8. SAVE ----
# ============================================================================= #

out_dir <- file.path(release_folder, "derived")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

components_tag <- paste(
  c(
    if (use_iadl)      "IADL",
    if (use_education) "EDU",
    if (use_atn)       "ATN",
    if (use_biomk)     "BIOMK",
    if (use_apoe)      "APOE",
    if (use_sucog)       "SUCOG",
    if (use_pacc)       "PACC",
    if (use_frailty)       "FRAILTY"
  ),
  collapse = "_"
)
if (!nzchar(components_tag)) components_tag <- "GENERAL_ONLY"

build_date <- format(Sys.Date(), "%Y%m%d")

base_name <- paste0(
  "ALFAplus_long_",
  components_tag, "_",
  release_docs_date,
  "_build_", build_date
)

out_file_csv <- file.path(out_dir, paste0(base_name, ".csv"))
write_csv(alfa_long, out_file_csv)
cat("\nSaved CSV to:\n  ", out_file_csv, "\n", sep="")

out_file_rdata <- file.path(out_dir, paste0(base_name, ".RData"))
save(alfa_long, file = out_file_rdata)
cat("Saved RData to:\n  ", out_file_rdata, "\n", sep="")

sources_file <- file.path(out_dir, paste0(base_name, "_sources.csv"))
write_csv(sources_log, sources_file)
cat("Saved sources log to:\n  ", sources_file, "\n", sep="")

if (use_biomk && !is.null(biomk_clean)) {
  biomk_meta_file <- file.path(out_dir, paste0("ALFAplus_BIOMK_metadata_", release_docs_date, "_build_", build_date, ".csv"))
  write_csv(biomk_clean$meta, biomk_meta_file)
  cat("Saved BIOMK metadata to:\n  ", biomk_meta_file, "\n", sep="")
  
  biomk_combo_file <- file.path(out_dir, paste0("ALFAplus_BIOMK_combos_", release_docs_date, "_build_", build_date, ".csv"))
  write_csv(biomk_clean$combo_summary, biomk_combo_file)
  cat("Saved BIOMK combo summary to:\n  ", biomk_combo_file, "\n", sep="")
}

cat("\nDone.\n")