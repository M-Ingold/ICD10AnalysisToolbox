#' Convert temporal ICD-10 codes into usable information for Cox regression
#'
#' @param dataframe_ID_BL A dataframe with a participant ID column. Dataframe gets returned with ICD-10 information is added to it
#' @param icd_table_long A long table with ICD-10 information, including the ID, ICD-10 code and date on which it was diagnosed
#' @param tbl_mortality_contact A dataframe containing participant ID, baseline visit date and last contact/death date columns.
#' @param icd_codes One or multiple ICD-10 codes of interest
#' @param variable_name Desired column names. Defaults to concatenated icd_codes
#' @param censor_time Time at which events are censored
#' @return Returns dataframe_ID_BL, with the outcome (event yes/no + time to event/censoring) and prevalence at baseline added.
#' @export

create_tte_icd <- function(dataframe_ID_BL,
                           icd_table_long,
                           tbl_mortality_contact,
                            icd_codes = c('I21'),
                            variable_name = NULL,
                            censor_time = 5) {

  ### Define new variable name
  if(is.null(variable_name)) variable_name <- tolower(paste0(icd_codes, collapse = '_'))

  ### Prepare matching string for ICD-10 codes
  icd_string <- paste0(icd_codes, collapse = '|')

  ### Select ICD-10 event dates
  icd_table_long <- icd_table_long[grepl(icd_string, icd_table_long$ICD10),]
  icd_table_long <- na.omit(icd_table_long[,c('ID', 'date_event')])
  dat <- icd_table_long[order(icd_table_long$ID, icd_table_long$date_event),]
  icd_table_long <- unique(icd_table_long)
  icd_table_long$date_event <- as.Date(icd_table_long$date_event)
  dat <- merge(icd_table_long[,c('ID', 'date_event')], tbl_mortality_contact[, c('ID', 'date_baseline', 'date_contact_loss')], all.x = T, all.y = T, by = 'ID')
  # dat <- merge(tbl_mortality_contact, dat, all.x = T, all.y = T, by = 'ID')
  dat$event <- ifelse(dat$date_event > dat$date_baseline, 1, 0)
  dat$hist <- ifelse(dat$date_event < dat$date_baseline, 1, 0)

  tbl_mortality_contact$time_contact_loss <- as.numeric(difftime(tbl_mortality_contact$date_contact_loss, tbl_mortality_contact$date_baseline, units = "days"))/365.25

  tte <- data.frame(ID=tbl_mortality_contact$ID, event=NA, hist=NA, time=NA)

  for (i in 1:nrow(tte)) {
    id <- tte$ID[i]
    if(!id %in% icd_table_long$ID) {
      tte$event[i] <- 0
      tte$hist[i] <- 0
      tte$time[i] <- tbl_mortality_contact$time_contact_loss[i]
    } else {
      dat_temp <- dat[dat$ID == id,]

      tte$event[i] <- max(dat_temp$event)
      tte$hist[i] <- max(dat_temp$hist)
      tte$time[i] <- tbl_mortality_contact$time_contact_loss[i]
      if(1 %in% dat_temp$event) {
        tte$time[i] <- min(as.numeric(difftime(dat_temp$date_event[dat_temp$event == 1], dat_temp$date_baseline[dat_temp$event == 1]))/365.25)
      }
    }
  }


  ### Censor after specific time
  if(!is.null(censor_time)) {
    tte$event[tte$time >= censor_time] <- 0
    tte$time[tte$time >= censor_time] <- censor_time
  }


  tte[,paste0(variable_name, '_event')] <- tte$event
  tte[,paste0(variable_name, '_hist')]  <- tte$hist
  tte[,paste0(variable_name, '_time')]  <- tte$time

  #tte$v11_ID01 <- tte$ID

  dataframe_ID_BL <- merge(dataframe_ID_BL, tte[, c('ID', paste0(variable_name, c('_event', '_hist', '_time')))],
                  all.x = T, suffixes = c('', '_icd'))

  return(dataframe_ID_BL)
}




