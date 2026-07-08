# ICD-10-CM to PheCode mapping

Mapping of ICD-10-CM codes to Phecodes, as well as descriptions and
categories, adapted from Phecode Map 1.2

## Usage

``` r
phecodes
```

## Format

### `phecodes`

A data frame with 90,768 rows and 6 columns:

- ICD:

  ICD-10 code

- Flag:

  ICD version (always 10)

- ICDString:

  ICD-10 description

- Phecode:

  Phecode as number

- PhecodeString:

  Phecode description

- PhecodeCategory:

  Disease category

## Source

PheWAS Resources
<https://phewascatalog.org/phewas/_w_7b1256e92bb84559bf1c493fd51cf540/data/Phecode_map_v1_2_icd9_icd10cm.csv.zip>
(downloaded June 2026)

## References

Denny, J., Bastarache, L., Ritchie, M. et al. Systematic comparison of
phenome-wide association study of electronic medical record data and
genome-wide association study data. Nat Biotechnol 31, 1102–1111 (2013).
([doi:10.1038/nbt.2749](https://doi.org/10.1038/nbt.2749) )
