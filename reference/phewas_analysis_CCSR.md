# Determine the association of a variable of interest by running Cox regression on time to event information and logistic regression on prevalence at baseline. Generates PheWAS plots to assess the strenth of associations

Determine the association of a variable of interest by running Cox
regression on time to event information and logistic regression on
prevalence at baseline. Generates PheWAS plots to assess the strenth of
associations

## Usage

``` r
phewas_analysis_CCSR(
  cohort_df,
  outcome_df = CCSR_outcomes,
  target,
  adj_vars = c("age", "sex"),
  save_as = "CCSR",
  n_event_cutoff = 10,
  exclude_categories = c("Congenital malformations", "Injury and poisoning consequences",
    "Health status factors", "Pregnancy and childbirth", "Perinatal conditions",
    "Unacceptable principal diagnosis"),
  one_to_one_mapping = TRUE,
  exclude_history = TRUE,
  return_tables = TRUE,
  annotation_size = 5,
  skip_colors = NULL
)
```

## Arguments

- cohort_df:

  A dataframe including participant ID and adjustment columns

- outcome_df:

  Result from CCSR_tte function

- target:

  Analysis variable column name

- adj_vars:

  List of variables to adjust models for

- save_as:

  Prefix for saved plot names

- n_event_cutoff:

  Diseases with fewer events are omitted

- exclude_categories:

  CCSR categories to omit

- exclude_history:

  If TRUE, only considers events with no history before baseline
  (default)

- skip_colors:

  Omit colors from the palette, useful to sync plots between cohorts

## Value

Returns two dataframes with regression results and generates PheWAS
plots
