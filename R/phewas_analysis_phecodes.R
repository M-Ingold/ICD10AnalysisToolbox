#' Determine the association of a variable of interest by running Cox regression on time to event information and logistic regression on prevalence at baseline.
#' Generates PheWAS plots to assess the strenth of associations
#'
#' @param cohort_df A dataframe including participant ID and adjustment columns
#' @param outcome_df Result from phecode_tte function
#' @param adj_vars List of variables to adjust models for
#' @param target Analysis variable column name
#' @param exclude_history If TRUE, only considers events with no history before baseline (default)
#' @param save_as Prefix for saved plot names
#' @param n_event_cutoff Diseases with fewer events are omitted
#' @param exclude_categories Phecode categories to omit
#' @param skip_colors Omit colors from the palette, useful to sync plots between cohorts
#' @param add_events Add a user-generated event to the analysis. Needs to exist in the outcome_df or cohort_df dataframe
#' @return Returns two dataframes with regression results and generates PheWAS plots
#' @export
#' @importFrom PheWAS "phenotypeManhattan"

phewas_analysis_phecodes <- function(cohort_df = ghs_eps, outcome_df = phecode_outcomes, adj_vars = c("age", "sex"),
                   exclude_categories = c("congenital anomalies", "injuries & poisonings",
                                          "pregnancy complications", "symptoms", "other"),
                   save_as="phecode", n_event_cutoff =10, target = target.scaled, exclude_history = TRUE,
                   add_events = NULL, add_events_name = NULL, add_events_category = NULL,
                   return_tables = FALSE, annotation_size = 5, skip_colors = NULL, plot_pdf=FALSE){

  source("N:/Transfer/ing1m/biosignature_pipeline/fast_glm2.R")

  cohort_df$target.scaled <- scale(cohort_df[[target]])[, ]
  cohort_df <- cohort_df[!is.na(cohort_df$target.scaled),]

  outcome_df <- outcome_df[outcome_df$ID %in% cohort_df$ID,]

  if (exclude_history == TRUE) {
    outcome_df <- outcome_df |>
      dplyr::mutate(across(ends_with("_event"),
                    ~ ifelse(get(stringr::str_replace(dplyr::cur_column(), "_event$", "_hist")) == 1,
                             NA,
                             .)))
  }



  cohort_eps_phecodes <- merge(cohort_df, outcome_df, by = "ID")


  phecodes <- read.csv("N:/Transfer/ing1m/pipeline_data/Phecode_map_v1_2_icd9_icd10cm_09_30_2024.csv")

  phecodes <- phecodes[phecodes$Flag==10,]
  phecodes$Phecode_char <- paste0("phecode_", gsub("\\.", "_", as.character(phecodes$Phecode)))

  # Add uncategorized traits to other. This includes e.g. "Other ill-defined and unknown causes of morbidity and mortality",
  # "Complications of surgical and medical procedures" or "Other tests"
  phecodes$PhecodeCategory[is.na(phecodes$PhecodeCategory)] <- "other"

  # Exclude categories as desired
  phecodes <- phecodes[!phecodes$PhecodeCategory %in% exclude_categories,]

  traits <- unique(phecodes$Phecode_char)
  traits_in_df <- gsub("_event", "", colnames(cohort_eps_phecodes))

  traits <- traits[traits %in% traits_in_df]



  # exclude events < 10
  n_event <- sapply(cohort_eps_phecodes[paste0(traits, "_event")], function(x) sum(x==1, na.rm=T))
  exclude <- names(which(n_event < n_event_cutoff))
  traits <- setdiff(traits, gsub("_event", "", exclude))
  traits <- c(traits, add_events)

  if (!length(traits)==0) {


  input <- data.frame(HR=0, negOR=0, LCI=0, UCI=0, p=0, nevent=0)


  for (i in 1:length(traits)) {

    form <- reformulate(c('target.scaled', adj_vars),
                        paste0('survival::Surv(', traits[i], '_time,', traits[i], '_event)'))

    input[i,] <- c(summary(survival::coxph(form, cohort_eps_phecodes))$conf.int[1,], summary(survival::coxph(form, cohort_eps_phecodes))[["coefficients"]][1,5], summary(survival::coxph(form, cohort_eps_phecodes))$nevent)


  }

  # merge with phecodes dataframe for additional info
  input$Phecode_char <- traits

  phecodes_unique <- phecodes[!duplicated(phecodes$Phecode_char),]

  input <- merge(input, phecodes_unique[c("Phecode_char", "Phecode","PhecodeString","PhecodeCategory")], by="Phecode_char", all.x = T)

  # handle manual addition of outcomes
  if (!is.null(add_events)) {
    add_phecode <- data.frame(Phecode_char=add_events,
                              PhecodeString=add_events_name,
                              PhecodeCategory=add_events_category)
    input[1:length(add_events),8:10] <- add_phecode

  }

  input$p_adj <- p.adjust(input$p, method = "fdr")
  outcome_table <- input

  # adjust for plotting
  input$OR <- input$HR
  input$phenotype <- input$PhecodeString
  input$group <- input$PhecodeCategory
  input$groupnum <- as.numeric(as.factor(input$PhecodeCategory))
  input$color <- as.factor(input$groupnum)
  input$description <- input$phenotype

  colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                         "#000000","#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                         "#000000", "#E69F00", "#56B4E9", "#009E73",
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  if (!is.null(skip_colors)) {colorBlindBlack8 <- colorBlindBlack8[-skip_colors]}

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

  if (plot_pdf) {
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


  }


  # Phecode prevalence ----

  traits <- unique(phecodes$Phecode_char)
  traits_in_df <- gsub("_hist", "", colnames(cohort_eps_phecodes))
  traits <- traits[traits %in% traits_in_df]

  # exclude events < 10
  n_event <- sapply(cohort_eps_phecodes[paste0(traits, "_hist")], function(x) sum(x==1, na.rm=T))
  exclude <- names(which(n_event < n_event_cutoff))
  traits <- setdiff(traits, gsub("_hist", "", exclude))

  if (!length(traits)==0) {


  # phecode_prevalence <- matrix(nrow = length(traits), ncol = 5)
  # list_marker <- list()


    model5 <- fast_glm2(response = paste0(traits, "_hist"),
                         response_name = traits,
                         Xvar = 'target.scaled',
                         dataset = cohort_eps_phecodes,
                         adj_vars = adj_vars,
                         family = 'poisson',
                         robust = TRUE,
                         sort = FALSE,
                         exponentiate = TRUE)



    # merge with phecodes dataframe for additional info
    model5$Phecode_char <- traits

    phecodes_unique <- phecodes[!duplicated(phecodes$Phecode_char),]

    input <- merge(model5, phecodes_unique[c("Phecode_char", "Phecode","PhecodeString","PhecodeCategory")], by="Phecode_char")

    prevalence_table <- input

    input$OR <- input$estimate
    input$phenotype <- input$PhecodeString
    input$group <- input$PhecodeCategory
    input$groupnum <- as.numeric(as.factor(input$PhecodeCategory))
    input$color <- as.factor(input$groupnum)
    input$description <- input$phenotype

    input$p_adj <- p.adjust(input$p, method = "fdr")


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


