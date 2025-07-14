#' Map temporal ICD-10 codes into outcomes for Cox regression using the CCSR convention
#'
#' @param dataframe_ID_BL A dataframe with a participant ID column. Dataframe gets returned with ICD-10 information is added to it
#' @param icd_table_long A long table with ICD-10 information, including the ID, ICD-10 code and date on which it was diagnosed
#' @param tbl_mortality_contact A dataframe containing participant ID, baseline visit date and last contact/death date columns.
#' @param censor_time Time at which events are censored
#' @param save_table Decide if the dataframe is written into the current directory
#' @param one_to_one_mapping If false, will map all CCSR-ICD-10 pairs instead of the single default CCSR-ICD-10 pair
#' @return Returns a dataframe, with all outcomes (event yes/no + time to event/censoring) and prevalences at baseline accoding to CCSR convention
#' @export


time_to_event_CCSR <- function(dataframe_ID_BL,
                        icd_table_long,
                        tbl_mortality_contact,
                        censor_time = 5,
                        one_to_one_mapping = TRUE,
                        save_table = FALSE) {

  ### create a loop going through all Phecodes, grabbing all outcomes with associated ICD10 codes

  DXCCSR <- read.csv("N:/Transfer/ing1m/PheCode_CCSR_Outcomes_GHS_MyoVasc/DXCCSR_v2025-1.csv")

  # get ICD10 format to match our data
  DXCCSR[] <- lapply(DXCCSR, function(x) gsub("'", "", x))
  DXCCSR$ICD10_reformat <- sub("(.{3})", "\\1.", DXCCSR$X.ICD.10.CM.CODE.)
  DXCCSR$tissue <- substr(DXCCSR$X.Default.CCSR.CATEGORY.IP., 1, 3)


  if (!one_to_one_mapping) {
    DXCCSR_long <- rbind(as.matrix(DXCCSR[c("ICD10_reformat", "X.CCSR.CATEGORY.1.", "CCSR.CATEGORY.1.DESCRIPTION")]),
                                    as.matrix(DXCCSR[c("ICD10_reformat", "X.CCSR.CATEGORY.2.", "CCSR.CATEGORY.2.DESCRIPTION")]),
                                    as.matrix(DXCCSR[c("ICD10_reformat", "X.CCSR.CATEGORY.3.", "CCSR.CATEGORY.3.DESCRIPTION")]),
                                    as.matrix(DXCCSR[c("ICD10_reformat", "CCSR.CATEGORY.4", "CCSR.CATEGORY.4.DESCRIPTION")]),
                                    as.matrix(DXCCSR[c("ICD10_reformat", "CCSR.CATEGORY.5", "CCSR.CATEGORY.5.DESCRIPTION")]),
                                    as.matrix(DXCCSR[c("ICD10_reformat", "CCSR.CATEGORY.6", "CCSR.CATEGORY.6.DESCRIPTION")]))

    DXCCSR_long <- as.data.frame(DXCCSR_long)
    DXCCSR_long$X.CCSR.CATEGORY.1.[DXCCSR_long$X.CCSR.CATEGORY.1.==" "] <- NA
    DXCCSR_long <- na.omit(DXCCSR_long)

    colnames(DXCCSR_long) <- c("ICD10_reformat", "X.Default.CCSR.CATEGORY.IP.", "Default.CCSR.CATEGORY.DESCRIPTION.IP")

    DXCCSR <- DXCCSR_long
  }

  DXCCSR_codes <- unique(DXCCSR$X.Default.CCSR.CATEGORY.IP.)
  DXCCSR_traits <- unique(DXCCSR$Default.CCSR.CATEGORY.DESCRIPTION.IP)

  fails <- list() # store outcomes that failed for whatever reason

  pb <- txtProgressBar(min = 1, max = length(DXCCSR_traits), style = 3)

  for (i in 1:length(DXCCSR_codes)) {
    icd <- DXCCSR[DXCCSR$X.Default.CCSR.CATEGORY.IP. == DXCCSR_codes[i],]$ICD10_reformat
    trait <- DXCCSR_codes[i]

    skip_to_next <- FALSE

    tryCatch(dataframe_ID_BL <- create_tte_icd(dataframe_ID_BL = dataframe_ID_BL,
                                       icd_table_long = icd_table_long,
                                       tbl_mortality_contact = tbl_mortality_contact,
                                       censor_time = censor_time,
                                       icd_codes = icd,
                                       variable_name = trait),
             error = function(e) { skip_to_next <<- TRUE})

    if(skip_to_next) { fails <- c(fails, trait) }

    if(skip_to_next) { next }

    #print(table(dataframe_ID_BL[paste0(trait, "_event")]))
    setTxtProgressBar(pb, i)

  }

  # subset for successful outcomes

  cols <- c("ID",
            sort(c(paste0(DXCCSR_codes, "_time"),
                   paste0(DXCCSR_codes, "_event"),
                   paste0(DXCCSR_codes, "_hist"))))
  cols <- cols[cols %in% colnames(dataframe_ID_BL)]
  dataframe_ID_BL_DXCCSR <- dataframe_ID_BL[cols]
  colnames(dataframe_ID_BL_DXCCSR)[2:length(cols)] <- paste("CCSR_", cols[2:length(cols)], sep = "")

  if (save_table) {
    write.csv(dataframe_ID_BL_DXCCSR, "tte_CCSR.csv")
  }

  print(paste("Failed for: ", paste(fails, collapse = " ")))
  return(dataframe_ID_BL_DXCCSR)

}
