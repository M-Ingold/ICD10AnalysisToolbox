# Determine the association of a variable of interest by running Cox regression on time to event information and logistic regression on prevalence at baseline. Generates PheWAS plots to assess the strenth of associations

Determine the association of a variable of interest by running Cox
regression on time to event information and logistic regression on
prevalence at baseline. Generates PheWAS plots to assess the strenth of
associations

## Usage

``` r
phewas_analysis_phecodes(
  cohort_df = ghs_eps,
  outcome_df = phecode_outcomes,
  adj_vars = c("age", "sex"),
  exclude_categories = c("congenital anomalies", "injuries & poisonings",
    "pregnancy complications", "symptoms", "other"),
  save_as = "phecode",
  n_event_cutoff = 10,
  target = target.scaled,
  exclude_history = TRUE,
  add_events = NULL,
  add_events_name = NULL,
  add_events_category = NULL,
  return_tables = FALSE,
  annotation_size = 5,
  skip_colors = NULL,
  plot_pdf = FALSE
)
```

## Arguments

- cohort_df:

  A dataframe including participant ID and adjustment columns

- outcome_df:

  Result from phecode_tte function

- adj_vars:

  List of variables to adjust models for

- exclude_categories:

  Phecode categories to omit

- save_as:

  Prefix for saved plot names

- n_event_cutoff:

  Diseases with fewer events are omitted

- target:

  Analysis variable column name

- exclude_history:

  If TRUE, only considers events with no history before baseline
  (default)

- add_events:

  Add a user-generated event to the analysis. Needs to exist in the
  outcome_df or cohort_df dataframe

- skip_colors:

  Omit colors from the palette, useful to sync plots between cohorts

## Value

Returns two dataframes with regression results and generates PheWAS
plots
