# Truncate suitability predictions based on an extrapolation value

Exclusion of suitability predictions in environmental conditions assumed
to be extrapolative.

## Usage

``` r
extra_truncate(suit, extra, threshold = 50, trunc_value = 0)
```

## Arguments

- suit:

  SpatRaster with suitability values

- extra:

  SpatRaster with extrapolation values preferable measured with
  extra_eval function

- threshold:

  numeric. Vector with one or more extrapolation values used for
  truncate suitability Default 50%

- trunc_value:

  numeric. Numeric value to be used to those cells assumed to be
  extrapolative

## Value

A SpatRaster object with truncated suitability values

## Details

Exclusion of suitability predictions in environmental conditions assumed
to be extrapolative. In this function it is possible to use any metric
measuring degree of extrapolation (e.g., MESS-Multivariate Environmental
Similarity Surfaces, EO-Environmental Overlap, MOP-Mobility-Oriented
Parity, EXDET-Extrapolation Detection, or AOA-Area of Applicability).
However, we recommend to use Shape approach (see
[`extra_eval`](https://sjevelazco.github.io/flexsdm/reference/extra_eval.md),
and [Velazco et al., 2023](https://doi.org/10.1111/ecog.06992)).

This function truncates suitability predictions assigning a given value,
generally 0 or NA. Usage trunc_value = NA. Default 0.

To those cells assumed to be extrapolative, i.e., higher than a given
threshold of a given extrapolation metric.

See this [vignette at flexsdm
website](https://sjevelazco.github.io/flexsdm/articles/v06_Extrapolation_example.html)
for further details about Shape metric, model truncation, and tools to
explore model extrapolation.

## See also

[`extra_eval`](https://sjevelazco.github.io/flexsdm/reference/extra_eval.md),
[`p_extra`](https://sjevelazco.github.io/flexsdm/reference/p_extra.md),
[`p_pdp`](https://sjevelazco.github.io/flexsdm/reference/p_pdp.md),
[`p_bpdp`](https://sjevelazco.github.io/flexsdm/reference/p_bpdp.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# see examples in extra_eval function
} # }
```
