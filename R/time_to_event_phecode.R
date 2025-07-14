#' Map temporal ICD-10 codes into outcomes for Cox regression using the phecode convention
#'
#' @param dataframe_ID_BL A dataframe with a participant ID column. Dataframe gets returned with ICD-10 information is added to it
#' @param icd_table_long A long table with ICD-10 information, including the ID, ICD-10 code and date on which it was diagnosed
#' @param tbl_mortality_contact A dataframe containing participant ID, baseline visit date and last contact/death date columns.
#' @param censor_time Time at which events are censored
#' @param save_table Decide if the dataframe is written into the current directory
#' @return Returns a dataframe, with all outcomes (event yes/no + time to event/censoring) and prevalences at baseline accoding to phecode convention
#' @export

time_to_event_phecode <- function(dataframe_ID_BL,
         icd_table_long,
         tbl_mortality_contact,
         censor_time = 5,
         save_table = FALSE) {
### create a loop going through all Phecodes, grabbing all outcomes with associated ICD10 codes

phecodes <- read.csv("N:/Transfer/ing1m/PheCode_CCSR_Outcomes_GHS_MyoVasc/Phecode_map_v1_2_icd9_icd10cm_09_30_2024.csv")

phecodes <- phecodes[phecodes$Flag==10,]

phecodes$Phecode_char <- gsub("\\.", "_", as.character(phecodes$Phecode))
phecodes_list <- unique(phecodes$Phecode_char)

fails <- list() # store outcomes that failed for whatever reason

pb <- txtProgressBar(min = 1, max = length(phecodes_list), style = 3)

for (i in 1:length(phecodes_list)) {
  icd <- phecodes[phecodes$Phecode_char == phecodes_list[i],]$ICD
  trait <- phecodes_list[i]

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

  # Print progress
  setTxtProgressBar(pb, i)

}

# subset for phecode outcomes

cols <- c("ID",
          sort(c(paste0(phecodes_list, "_time"),
                 paste0(phecodes_list, "_event"),
                 paste0(phecodes_list, "_hist"))))
cols <- cols[cols %in% colnames(dataframe_ID_BL)]
dataframe_ID_BL_phecodes <- dataframe_ID_BL[cols]

colnames(dataframe_ID_BL_phecodes)[2:length(cols)] <- paste("phecode_", cols[2:length(cols)], sep = "")

if (save_table) {
  write.csv(dataframe_ID_BL_phecodes, "tte_phecodes.csv")
}

print(paste("Failed for: ", paste(fails, collapse = " ")))
return(dataframe_ID_BL_phecodes)

}
