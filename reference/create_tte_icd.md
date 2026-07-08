# Convert temporal ICD-10 codes into usable information for Cox regression

Convert temporal ICD-10 codes into usable information for Cox regression

## Usage

``` r
create_tte_icd(
  dataframe_ID_BL,
  icd_table_long,
  tbl_mortality_contact,
  icd_codes = c("I21"),
  variable_name = NULL,
  censor_time = 5,
  no_FU_as_NA = FALSE
)
```

## Arguments

- dataframe_ID_BL:

  A dataframe with a participant ID column. Dataframe gets returned with
  ICD-10 information is added to it

- icd_table_long:

  A long table with ICD-10 information, including the ID, ICD-10 code
  and date on which it was diagnosed

- tbl_mortality_contact:

  A dataframe containing participant ID, baseline visit date and last
  contact/death date columns.

- icd_codes:

  One or multiple ICD-10 codes of interest

- variable_name:

  Desired column names. Defaults to concatenated icd_codes

- censor_time:

  Time at which events are censored

## Value

Returns dataframe_ID_BL, with the outcome (event yes/no + time to
event/censoring) and prevalence at baseline added.
