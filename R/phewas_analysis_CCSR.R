#' Determine the association of a variable of interest by running Cox regression on time to event information and logistic regression on prevalence at baseline.
#' Generates PheWAS plots to assess the strenth of associations
#'
#' @param cohort_df A dataframe including participant ID and adjustment columns
#' @param outcome_df Result from CCSR_tte function
#' @param adj_vars List of variables to adjust models for
#' @param target Analysis variable column name
#' @param exclude_history If TRUE, only considers events with no history before baseline (default)
#' @param save_as Prefix for saved plot names
#' @param n_event_cutoff Diseases with fewer events are omitted
#' @param exclude_categories CCSR categories to omit
#' @param skip_colors Omit colors from the palette, useful to sync plots between cohorts
#' @return Returns two dataframes with regression results and generates PheWAS plots
#' @export
#' @importFrom PheWAS "phenotypeManhattan"

phewas_analysis_CCSR <- function(#CCSR_codes = CCSR_codes,
                        cohort_df, outcome_df = CCSR_outcomes, target,
                   adj_vars = c("age", "sex"), save_as="CCSR", n_event_cutoff =10,
                   exclude_categories = c("Congenital malformations", "Injury and poisoning consequences", "Health status factors",
                                          "Pregnancy and childbirth", "Perinatal conditions", "Unacceptable principal diagnosis"),
                   one_to_one_mapping = TRUE,
                   exclude_history = TRUE, return_tables = TRUE, annotation_size = 5, skip_colors = NULL){

  cohort_df$target.scaled <- scale(cohort_df[[target]])[, ]
  cohort_df <- cohort_df[!is.na(cohort_df$target.scaled),]

  outcome_df <- outcome_df[outcome_df$ID %in% cohort_df$ID,]

  if (exclude_history == TRUE) {
    outcome_df <- outcome_df |>
      dplyr::mutate(dplyr::across(dplyr::ends_with("_event"),
                    ~ ifelse(get(stringr::str_replace(dplyr::cur_column(), "_event$", "_hist")) == 1,
                             NA,
                             .)))
  }


  cohort_eps_CCSR_codes <- merge(cohort_df, outcome_df, by = "ID")


  # grab and adjust phecodes file to match variable names
  DXCCSR <- read.csv("N:/Transfer/ing1m/PheCode_CCSR_Outcomes_GHS_MyoVasc/DXCCSR_v2025-1.csv")

  # get ICD10 format to match our data
  DXCCSR[] <- lapply(DXCCSR, function(x) gsub("'", "", x))
  DXCCSR$ICD10_reformat <- sub("(.{3})", "\\1.", DXCCSR$X.ICD.10.CM.CODE.)
  DXCCSR$category <- substr(DXCCSR$X.Default.CCSR.CATEGORY.IP., 1, 3)

  if (!one_to_one_mapping) {
    DXCCSR_long <- rbind(as.matrix(DXCCSR[c("ICD10_reformat", "X.CCSR.CATEGORY.1.", "CCSR.CATEGORY.1.DESCRIPTION", "category")]),
                         as.matrix(DXCCSR[c("ICD10_reformat", "X.CCSR.CATEGORY.2.", "CCSR.CATEGORY.2.DESCRIPTION", "category")]),
                         as.matrix(DXCCSR[c("ICD10_reformat", "X.CCSR.CATEGORY.3.", "CCSR.CATEGORY.3.DESCRIPTION", "category")]),
                         as.matrix(DXCCSR[c("ICD10_reformat", "CCSR.CATEGORY.4", "CCSR.CATEGORY.4.DESCRIPTION", "category")]),
                         as.matrix(DXCCSR[c("ICD10_reformat", "CCSR.CATEGORY.5", "CCSR.CATEGORY.5.DESCRIPTION", "category")]),
                         as.matrix(DXCCSR[c("ICD10_reformat", "CCSR.CATEGORY.6", "CCSR.CATEGORY.6.DESCRIPTION", "category")]))

    DXCCSR_long <- as.data.frame(DXCCSR_long)
    DXCCSR_long$X.CCSR.CATEGORY.1.[DXCCSR_long$X.CCSR.CATEGORY.1.==" "] <- NA
    DXCCSR_long <- na.omit(DXCCSR_long)

    colnames(DXCCSR_long) <- c("ICD10_reformat", "X.Default.CCSR.CATEGORY.IP.", "Default.CCSR.CATEGORY.DESCRIPTION.IP", "category")

    DXCCSR <- DXCCSR_long
  }


  # get better category descriptions
  description_mapping <- c(
    BLD = "Blood and immune disorders",
    CIR = "Circulatory system diseases",
    DEN = "Dental diseases",
    DIG = "Digestive system diseases",
    EAR = "Ear and mastoid diseases",
    END = "Endocrine and metabolic diseases",
    EXT = "External causes morbidity",
    EYE = "Eye and adnexa diseases",
    FAC = "Health status factors",
    GEN = "Genitourinary system diseases",
    INF = "Infectious and parasitic diseases",
    INJ = "Injury and poisoning consequences",
    MAL = "Congenital malformations",
    MBD = "Mental and behavioral disorders",
    MUS = "Musculoskeletal system diseases",
    NEO = "Neoplasms",
    NVS = "Nervous system diseases",
    PNL = "Perinatal conditions",
    PRG = "Pregnancy and childbirth",
    RSP = "Respiratory system diseases",
    SKN = "Skin and subcutaneous diseases",
    SYM = "Symptoms and abnormal findings",
    XXX = "Unacceptable principal diagnosis"
  )

  # Add the new column with descriptions
  DXCCSR <- DXCCSR |>
    dplyr::mutate(category = description_mapping[category])

  DXCCSR_codes <- unique(DXCCSR$X.Default.CCSR.CATEGORY.IP.)
  DXCCSR_traits <- unique(DXCCSR$Default.CCSR.CATEGORY.DESCRIPTION.IP)

  # Exclude categories as desired
  DXCCSR <- DXCCSR[!DXCCSR$category %in% exclude_categories,]


  traits <- unique(DXCCSR_codes)
  traits_in_df <- gsub("_event", "", colnames(cohort_eps_CCSR_codes))
  traits_in_df <- gsub("CCSR_", "", traits_in_df)

  traits <- traits[traits %in% traits_in_df]

  # exclude events < 10
  n_event <- sapply(cohort_eps_CCSR_codes[paste0("CCSR_", traits, "_event")], function(x) sum(x==1, na.rm=T))
  exclude <- names(which(n_event < n_event_cutoff))
  exclude <- gsub("CCSR_", "", exclude)
  traits <- setdiff(traits, gsub("_event", "", exclude))


  if (!length(traits)==0) {


    input <- data.frame(HR=0, negOR=0, LCI=0, UCI=0, p=0, nevent=0)


    for (i in 1:length(traits)) {

      form <- reformulate(c('target.scaled', adj_vars),
                          paste0('survival::Surv(CCSR_', traits[i], '_time,', 'CCSR_',traits[i], '_event)'))

      input[i,] <- c(summary(survival::coxph(form, cohort_eps_CCSR_codes))$conf.int[1,], summary(survival::coxph(form, cohort_eps_CCSR_codes))[["coefficients"]][1,5], summary(survival::coxph(form, cohort_eps_CCSR_codes))$nevent)


    }

  # merge with phecodes dataframe for additional info
  input$DXCCSR_code <- traits
  #input$PhecodeString <- gsub("_", " ", input$PhecodeString)

  DXCCSR_unique <- DXCCSR[!duplicated(DXCCSR$X.Default.CCSR.CATEGORY.IP.),]

  input <- merge(input, DXCCSR_unique[c("X.Default.CCSR.CATEGORY.IP.", "category","Default.CCSR.CATEGORY.DESCRIPTION.IP")], by.x= "DXCCSR_code", by.y="X.Default.CCSR.CATEGORY.IP.")

  input$p_adj <- p.adjust(input$p, method = "fdr")
  outcome_table <- input

  # adjust for plotting
  input$OR <- input$HR
  input$phenotype <- input$Default.CCSR.CATEGORY.DESCRIPTION.IP
  input$group <- input$category
  input$groupnum <- as.numeric(as.factor(input$category))
  input$color <- as.factor(input$groupnum)
  input$description <- input$phenotype



  colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                         "#000000","#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                         "#000000", "#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  if (!is.null(skip_colors)) {colorBlindBlack8 <- colorBlindBlack8[-skip_colors]} # harmonize colors between cohorts, for example

  png(paste0(save_as, "_outcomes.png"), width = 2000, height = 1000, res=120)
  print(
    PheWAS::phenotypeManhattan(input, annotate.phenotype.description = T, use.color = T,
                               annotate.level = sort(input$p)[15], # annotate top 15
                               #annotate.level = max(input$p[input$p_adj<0.05]), # annotate all FDR significant
                               significant.line= 0.05/length(input$p), # needs to be specified, otherwise it takes the suggestive line p-value/number of tests
                               suggestive.line = ifelse(max(input$p[input$p_adj<0.05]) < 0.05/length(input$p),
                                                        0.05/length(input$p)-0.1, max(input$p[input$p_adj<0.05])), # blue line at largest FDR-significant p-value or just below Bonferroni line
                               color.palette = colorBlindBlack8, sort.by.category.value=T,
                               OR.direction = T, size.x.labels=14, size.y.labels=14,
                               max.y = ifelse(-log10(min(input$p, na.rm = T)) > 5, -log10(min(input$p, na.rm = T))+0.1, 5),
                               annotate.size=annotation_size)+
      ggplot2::theme(plot.margin = ggplot2::margin(0.5,1.5,0.1,0.1, "cm")
    )
  )
  dev.off()


  pdf(paste0(save_as, "_outcomes.pdf"), width = 15, height = 7.5)
  print(
    PheWAS::phenotypeManhattan(input, annotate.phenotype.description = T, use.color = T,
                               annotate.level = sort(input$p)[15], # annotate top 15
                               #annotate.level = max(input$p[input$p_adj<0.05]), # annotate all FDR significant
                               significant.line= 0.05/length(input$p), # needs to be specified, otherwise it takes the suggestive line p-value/number of tests
                               suggestive.line = ifelse(max(input$p[input$p_adj<0.05]) < 0.05/length(input$p),
                                                        0.05/length(input$p)-0.1, max(input$p[input$p_adj<0.05])), # blue line at largest FDR-significant p-value or just below Bonferroni line
                               color.palette = colorBlindBlack8, sort.by.category.value=T,
                               OR.direction = T, size.x.labels=14, size.y.labels=14,
                               max.y = ifelse(-log10(min(input$p, na.rm = T)) > 5, -log10(min(input$p, na.rm = T))+0.1, 5),
                               annotate.size=annotation_size)+
      ggplot2::theme(plot.margin = ggplot2::margin(0.5,1.5,0.1,0.1, "cm")
      )
  )
  dev.off()

  }

  # Phecode prevalence ----

  traits <- unique(DXCCSR_codes)
  traits_in_df <- gsub("_hist", "", colnames(cohort_eps_CCSR_codes))
  traits_in_df <- gsub("CCSR_", "", traits_in_df)
  traits <- traits[traits %in% traits_in_df]

  # exclude events < 10
  n_event <- sapply(cohort_eps_CCSR_codes[paste0('CCSR_', traits, "_hist")], function(x) sum(x==1, na.rm=T))
  exclude <- names(which(n_event < n_event_cutoff))
  exclude <- gsub("_hist", "", exclude)
  exclude <- gsub("CCSR_", "", exclude)

  traits <- setdiff(traits, gsub("_hist", "", exclude))

  if (!length(traits)==0) {


  # phecode_prevalence <- matrix(nrow = length(traits), ncol = 5)
  # list_marker <- list()


    model5 <- fast_glm2(response = paste0("CCSR_", traits, "_hist"),
                         response_name = traits,
                         Xvar = 'target.scaled',
                         dataset = cohort_eps_CCSR_codes,
                         adj_vars = adj_vars,
                         family = 'poisson',
                         robust = TRUE,
                         sort = FALSE,
                         exponentiate = TRUE)



    # merge with CCSR_codes dataframe for additional info
    model5$DXCCSR_code <- traits

    DXCCSR_unique <- DXCCSR[!duplicated(DXCCSR$X.Default.CCSR.CATEGORY.IP.),]

    input <- merge(model5, DXCCSR_unique[c("X.Default.CCSR.CATEGORY.IP.", "category","Default.CCSR.CATEGORY.DESCRIPTION.IP")], by.x= "DXCCSR_code", by.y="X.Default.CCSR.CATEGORY.IP.")

    input$p_adj <- p.adjust(input$p, method = "fdr")
    prevalence_table <- input

    input$OR <- input$estimate
    input$phenotype <- input$Default.CCSR.CATEGORY.DESCRIPTION.IP
    input$group <- input$category
    input$groupnum <- as.numeric(as.factor(input$category))
    input$color <- as.factor(input$groupnum)
    input$description <- input$phenotype

    colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                           "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                           "#000000","#E69F00", "#56B4E9", "#009E73",
                           "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                           "#000000", "#E69F00", "#56B4E9", "#009E73",
                           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    png(paste0(save_as, "_comorbidities.png"), width = 2000, height = 1000, res=120)
    print(
      PheWAS::phenotypeManhattan(input, annotate.phenotype.description = T, use.color = T,
                                 annotate.level = sort(input$p)[15], # annotate top 15
                                 #annotate.level = max(input$p[input$p_adj<0.05]), # annotate all FDR significant
                                 significant.line= 0.05/length(input$p), # needs to be specified, otherwise it takes the suggestive line p-value/number of tests
                                 suggestive.line = ifelse(max(input$p[input$p_adj<0.05]) < 0.05/length(input$p),
                                                          0.05/length(input$p)-0.01, max(input$p[input$p_adj<0.05])), # blue line at largest FDR-significant p-value or just below Bonferroni line
                                 color.palette = colorBlindBlack8, sort.by.category.value=T,
                                 OR.direction = T, size.x.labels=14, size.y.labels=14,
                                 max.y = ifelse(-log10(min(input$p, na.rm = T)) > 5, -log10(min(input$p, na.rm = T))+0.1, 5)
      )

    )
    dev.off()

  }

  if (return_tables) {
    return(list(outcome_table, prevalence_table))
  }

}


