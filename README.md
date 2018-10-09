
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multitestr

`multitestr` implements bootstrap-based multiple testing adjustments for
hypothesis tests conducted on linear regression coefficients. These
adjustments are often used in social science experiments where many
outcomes or subgroups are analyzed and there is a risk of falsely
rejecting the null of a particular hypothesis purely by chance. While
`R`’s built-in function `p.adjust` performs variety of multiple testing
adjustments such as the Bonferroni or Holm corrections, these approaches
can be quite conservative when test statistics are correlated. The
booststrap-based adjustments used in this package allow for correlations
across test-statistics and hence can be more powerful.

Currently, we implement the Westfall and Young (1993) bootstrap based
procedure to control the *familywise error rate* (FWER) for a family of
hypotheses. Users may specify either the “pairs” non-parametric
bootstrap or the Wild parameteric bootstrap.

## Installation

You can install the current version of `multitestr` from
[Github](https://github.com/fdhidalgo/multitestr) with:

``` r
install.packages("devtools")
devtools::install_github("fdhidalgo/multitestr")
```

## Example

To demonstrate `multitestr`‘s functionality, we reproduce (up to
sampling error) some results from Casey et al (2012). The study
estimates the effect of a \`\`community-driven development’’ program on
12 composite outcomes in Sierra Leone.

The user must specify a list of `full_formulas` and `null_formulas`,
where the `full_formulas` are the full regression models being tested
and `null_formulas` are models under the null hypothesis. In this
particular case, the `full_formulas` are formulas with 12 different
outcome variables (`z_score_1`…`z_score_12`), the treatment variable
(`t`), pre-treatment covariates (`tothhs` and `road`), and strata fixed
effects (`ward`). The `null_formulas` are the corresponding models
without a treatment variable, since the null hyptothesis is that the
treatment had no effect.

``` r
library("multitestr")
# Replicate Casey et al 2012, Table II Column 3
set.seed(123)
# Create a list of formulas with different dependent variables
F <- lapply(sprintf("zscore_%d ~ 0 + t + tothhs + road + ward", 1:12),
            as.formula)

#Run the Westfall and Young boostrap step down procedure            
pvals <- boot_stepdown(full_formulas = F,
                       null_formulas = lapply(F, update, . ~ . - t),
                       data = gobifo,
                       coef_list = "t", #this parameter specifies the coefficient of interest
                       nboots = 1000, 
                       parallel = FALSE, 
                       boot_type = "pairs")
#> ===========================================================================
#b
dplyr::mutate_if(pvals, is.numeric, round, 3)
#>       Hypothesis Variable bs_pvalues_unadjusted bs_pvalues_adjusted
#> 1   Hypothesis 1        t                 0.001               0.001
#> 2   Hypothesis 2        t                 0.001               0.001
#> 3   Hypothesis 3        t                 0.001               0.001
#> 4   Hypothesis 4        t                 0.738               0.981
#> 5   Hypothesis 5        t                 0.938               0.981
#> 6   Hypothesis 6        t                 0.134               0.677
#> 7   Hypothesis 7        t                 0.359               0.915
#> 8   Hypothesis 8        t                 0.452               0.915
#> 9   Hypothesis 9        t                 0.312               0.910
#> 10 Hypothesis 10        t                 0.043               0.348
#> 11 Hypothesis 11        t                 0.839               0.981
#> 12 Hypothesis 12        t                 0.357               0.915
```

`bs_pvalues_unadjusted` are the bootstrap-based p-values with no
adjustment for mutltiple testing. `bs_pvalues_adjusted` have been
adjusted for multiple testing using the algorithm described in Westfall
and Young (1993).

## References

Casey, Katherine, Rachel Glennerster, and Edward Miguel. 2012. Reshaping
Institutions: Evidence on Aid Impacts Using a Preanalysis Plan.
*Quarterly Journal of Economics.* 127(4):1755-1812.

Westfall, Peter H., and S. Stanley Young. 1993. *Resampling-Based
Multiple Testing: Examples and Methods for P-Value Adjustment.* New
York: John Wiley & Sons.
