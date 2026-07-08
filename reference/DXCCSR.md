# ICD-10-CM to Clinical Classifications Software Refined (CCSR) mapping

Mapping of ICD-10-CM codes to Clinical Classifications Software Refined
(CCSR) codes, as well as descriptions and categories, adapted from CCSR
v2025.1

## Usage

``` r
DXCCSR
```

## Format

### `DXCCSR`

A data frame with 90,768 rows and 6 columns:

- X.ICD.10.CM.CODE.:

  ICD-10 code

- ICD.10.CM.CODE.DESCRIPTION:

  ICD-10 description

- X.Default.CCSR.CATEGORY.IP.:

  Default CCSR ID

- Default.CCSR.CATEGORY.DESCRIPTION.IP:

  Default CCSR description

- X.CCSR.CATEGORY.1-5:

  Assigned category or categories, descending in order of relevance

- Rationale.for.Default.Assignment:

  Rationale for default CCSR code assignment

## Source

US Agency for Healthcare Research and Quality
<https://hcup-us.ahrq.gov/toolssoftware/ccsr/dxccsr.jsp> (Fiscal Year
2025 version)

## References

HCUP Clinical Classifications Software Refined (CCSR) for ICD-10-CM
diagnoses, v2025.1. Healthcare Cost and Utilization Project (HCUP).
Agency for Healthcare Research and Quality, Rockville, MD.
hcup-us.ahrq.gov/toolssoftware/ccsr/dxccsr.jsp. Accessed January 2025.
