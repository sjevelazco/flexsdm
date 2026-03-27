# Constraint of suitability based on extrapolation

Exclusion of suitability values less than a given extrapolation value
(EXPERIMENTAL)

## Usage

``` r
extra_exclude(suit, extra, threshold = 50)
```

## Arguments

- suit:

  SpatRaster with suitability values

- extra:

  SpatRaster with extrapolation values measured in percentage (output
  from extra_eval function)

- threshold:

  numeric. Vector with one or more values used for correct
  extrapolation. Default 50% (DOES THIS FUNCTION SET THE PROJECTED
  SUITABILITY VALUES LESS THAN THE THRESHOLD TO ZERO? UNCLEAR. PLEASE BE
  EXPLICIT)

## Value

A SpatRaster object with corrected suitability values

## See also

[`extra_eval`](https://sjevelazco.github.io/flexsdm/reference/extra_eval.md)

## Examples

``` r
if (FALSE) {
# see examples in extra_eval function
}
```
