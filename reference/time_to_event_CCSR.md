# Map temporal ICD-10 codes into outcomes for Cox regression using the CCSR convention

Map temporal ICD-10 codes into outcomes for Cox regression using the
CCSR convention

## Usage

``` r
time_to_event_CCSR(
  dataframe_ID_BL,
  icd_table_long,
  tbl_mortality_contact,
  censor_time = 5,
  one_to_one_mapping = TRUE,
  save_table = FALSE
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

- censor_time:

  Time at which events are censored

- one_to_one_mapping:

  If false, will map all CCSR-ICD-10 pairs instead of the single default
  CCSR-ICD-10 pair

- save_table:

  Decide if the dataframe is written into the current directory

## Value

Returns a dataframe, with all outcomes (event yes/no + time to
event/censoring) and prevalences at baseline accoding to CCSR convention
