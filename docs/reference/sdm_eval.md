# Calculate different model performance metrics

This function calculates threshold dependent and independent model
performance metrics.

## Usage

``` r
sdm_eval(p, a, bg = NULL, thr = NULL)
```

## Arguments

- p:

  numeric. Predicted suitability for presences

- a:

  numeric. Predicted suitability for absences

- bg:

  numeric. Predicted suitability for background points, used for BOYCE
  metric. If not provided (NULL), the Boyce index will be calculated
  using absence data instead. **Note:** Using absence data for the Boyce
  index is not standard practice and may result in inflated performance
  values. It is highly recommended to provide background points for a
  correct calculation. The Boyce index is calculated using the `boyce`
  function, which is an adaptation of the method implemented in the
  `ecospat` package.

- thr:

  character. Threshold criterion used to get binary suitability values
  (i.e. 0,1). Used for threshold-dependent performance metrics. It is
  possible to use more than one threshold type. A vector must be
  provided for this argument. The following threshold criteria are
  available:

  - lpt: The highest threshold at which there is no omission.

  - equal_sens_spec: Threshold at which the Sensitivity and Specificity
    are equal.

  - max_sens_spec: Threshold at which the sum of the Sensitivity and
    Specificity is the highest (aka threshold that maximizes the TSS).

  - max_jaccard: The threshold at which the Jaccard index is the
    highest.

  - max_sorensen: The threshold at which the Sorensen index is the
    highest.

  - max_fpb: The threshold at which FPB (F-measure on
    presence-background data) is the highest.

  - sensitivity: Threshold based on a specified Sensitivity value. Usage
    thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens'
    refers to Sensitivity value. If a sensitivity value is not
    specified, the default value is 0.9

  If more than one threshold type is used, concatenate threshold types,
  e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt',
  'max_sens_spec', 'sensitivity'). Function will use all thresholds if
  no threshold type is specified

## Value

a tibble with next columns

- threshold: threshold names

- thr_value: threshold values

- n_presences: number of presences

- n_absences: number of absences

- from TPR to CRPS: performance metrics

## Details

This function is used for evaluating different models approaches base on
the combination of presence-absences or presence-pseudo-absences and
background point data and suitability predicted by any model or flexsdm
modeling function families (fit\_, esm\_, and tune\_.)

### Performance Metrics Formulas

It calculates the next performance metric:

|  |  |  |
|----|----|----|
| Performance metric | Threshold dependent | Values ranges |
| TPR (True Positive Rate, also called Sensitivity) | yes | 0 - 1 |
| TNR (True Negative Rate, also called Specificity) | yes | 0 - 1 |
| W_TPR_TNR (Weighted TPR-TNR; Li et al. 2020) | yes | 0 - 1 |
| SORENSEN | yes | 0 - 1 |
| JACCARD | yes | 0 - 1 |
| FPB (F-measure on presence-background) | yes | 0 - 2 |
| OR (Omission Rate) | yes | 0 - 1 |
| TSS (True Skill Statistic) | yes | -1 - 1 |
| KAPPA | yes | 0 - 1 |
| MCC (Matthews Correlation Coefficient; Matthews 1975) | yes | -1 - 1 (1 is best) |
| AUC (Area Under Curve) | no | 0 - 1 |
| BOYCE (continuous Boyce index)\* | no | -1 - 1 |
| IMAE (Inverse Mean Absolute Error)\*\* | no | 0 - 1 |
| CRPS (Continuous Ranked Probability Score based on Brier Score, Brier 1950)\*\* | no | 0 - 1 |

\\ The continuous Boyce index is calculated based on presences and
background points. If background points are not provided, it will be
calculated using presences and absences, which is not standard and may
lead to misleading results. The code for calculating this metric is an
adaptation of the `ecospat` package.

\\\* IMAE and CRPS are calculated as 1-(Mean Absolute Error) and
1-(CRPS), respectively, in order to be consistent with the other metrics
where the higher the value of a given performance metric, the greater
the model's.

To define the formulas, the following components of the confusion matrix
are used:

- `tp`: True Positives (presences correctly predicted as presences)

- `tn`: True Negatives (absences correctly predicted as absences)

- `fp`: False Positives (absences incorrectly predicted as presences)

- `fn`: False Negatives (presences incorrectly predicted as absences)

- `np`: Number of presences (`length(p)`)

- `na`: Number of absences (`length(a)`)

The formulas are:

- **TPR (Sensitivity)**: \$\$TPR = \frac{tp}{tp + fn}\$\$

- **TNR (Specificity)**: \$\$TNR = \frac{tn}{tn + fp}\$\$

- **W_TPR_TNR** (Li et al. 2020): \$\$W\\TPR\\TNR = w \cdot TPR + (1-w)
  \cdot TNR\$\$, where \$\$w = \frac{na}{na + np}\$\$

- **SORENSEN**: \$\$Sorensen = \frac{2 \cdot tp}{fn + 2 \cdot tp +
  fp}\$\$

- **JACCARD**: \$\$Jaccard = \frac{tp}{fn + tp + fp}\$\$

- **FPB**: \$\$FPB = 2 \cdot Jaccard\$\$

- **OR**: \$\$OR = 1 - TPR\$\$

- **TSS**: \$\$TSS = TPR + TNR - 1\$\$

- **KAPPA**: \$\$KAPPA = \frac{Pr(a) - Pr(e)}{1 - Pr(e)}\$\$, where
  \$\$Pr(a) = \frac{tp+tn}{tp+tn+fp+fn}\$\$ and \$\$Pr(e) =
  \frac{(tp+fp)(tp+fn) + (fn+tn)(fp+tn)}{(tp+tn+fp+fn)^2}\$\$

- **MCC** (Matthews 1975): \$\$MCC = \frac{(tp \cdot tn) - (fp \cdot
  fn)}{\sqrt{(tp+fp)(tp+fn)(tn+fp)(tn+fn)}}\$\$

- **AUC**: Calculated as the Wilcoxon-Mann-Whitney U statistic, which is
  equivalent to the area under the ROC curve.

- **BOYCE**: The continuous Boyce index, which measures how model
  predictions differ from a random distribution of observed presences
  across the prediction gradient.

- **CRPS** (Brier 1950): For binary outcomes, this is calculated as
  \$\$1 - \frac{\sum(predicted - observed)^2}{N}\$\$, which is 1 minus
  the Brier Score.

- **IMAE**: \$\$IMAE = 1 - \frac{\sum\|predicted - observed\|}{N}\$\$,
  where N is the total number of records.

## References

- Brier GW. (1950) Verification of forecasts expressed in terms of
  probability. Monthly Weather Review 78(1): 1–3.
  https://doi.org/10.1175/1520-0493(1950)078\<0001:VOFEIT\>2.0.CO;2

- Li, J., Liu, H., & Li, L. (2020). A novel performance metric for
  imbalanced learning and its application in credit default prediction.
  Expert Systems with Applications, 152, 113382.
  https://doi.org/10.1016/j.eswa.2020.113382

- Matthews BW. (1975) Comparison of the predicted and observed secondary
  structure of T4 phage lysozyme. Biochim Biophys Acta (BBA) Protein
  Struct. 405(2):442–51. https://doi.org/10.1016/0005-2795(75)90109-9

## Examples

``` r
if (FALSE) { # \dontrun{
require(dplyr)

set.seed(0)
p <- rnorm(50, mean = 0.7, sd = 0.3) %>% abs()
p[p > 1] <- 1
p[p < 0] <- 0

set.seed(0)
a <- rnorm(50, mean = 0.3, sd = 0.2) %>% abs()
a[a > 1] <- 1
a[a < 0] <- 0

set.seed(0)
backg <- rnorm(1000, mean = 0.4, sd = 0.4) %>% abs()
backg[backg > 1] <- 1
backg[backg < 0] <- 0

# Function use without threshold specification
e <- sdm_eval(p, a)
e

# Function use with threshold specification
sdm_eval(p, a, thr = "max_sorensen")
sdm_eval(p, a, thr = c("lpt", "max_sens_spec", "max_jaccard"))
sdm_eval(p, a, thr = c("lpt", "max_sens_spec", "sensitivity"))
sdm_eval(p, a, thr = c("lpt", "max_sens_spec", "sensitivity", sens = "0.95"))

# Use of bg argument (it will only be used for calculating BOYCE index)
sdm_eval(p, a, thr = "max_sens_spec")
sdm_eval(p, a, thr = c("max_sens_spec"), bg = backg)

# If background will be used to calculate all other metrics
# background values can be used in "a" argument
sdm_eval(p, backg, thr = "max_sens_spec")
} # }
```
