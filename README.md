# fmlogit

**Fractional multinomial logit for shares that add up to one.**

The `fmlogit` package estimates how covariates relate to the allocation of a
total across several categories. It provides fitted shares, marginal and
discrete effects, robust inference, and lightweight plots.

The main methodological reference is
[Mullahy (2015)](https://doi.org/10.1515/jem-2012-0006).
[Wulff’s guide](https://doi.org/10.1177/1094428114560024)
to interpreting multinomial logit results
explains why predicted outcomes and partial effects are useful alongside
coefficient tables.

## What is the fractional multinomial logit model?

Many outcomes describe how a resource is divided: a household’s spending
across goods, a portfolio across asset classes, or a day across activities.
Each observation contains several shares, each between zero and one, with a
total of one. A person may allocate time to several activities on the same
day, producing shares such as 0.5, 0.3, and 0.2.

The model describes the **conditional mean** of each share: the average
allocation among observations with given covariate values. With observation
`i`, category `j`, and `J` categories, its mean function is

$$
p_{ij} = E(y_{ij} \mid x_i) = \frac{\exp(x_i' \beta_j)}{\sum_{h=1}^{J} \exp(x_i' \beta_h)},
\qquad \beta_1 = 0.
$$

Here `x_i` contains the covariates and an intercept. All categories enter the
same denominator, which makes their fitted shares add up to one. The first
category is the reference, or baseline; its coefficients are set to zero
to identify the model. Its fitted share and partial effects are still
estimated and can be substantial.

The same notation is used throughout the package:

| Name | Meaning |
| --- | --- |
| `N`, `J`, `K` | Numbers of retained observations, outcomes, and design columns, including the intercept. |
| `y` | `N`-by-`J` matrix of observed shares; the baseline outcome is first. |
| `x` | `N`-by-`K` design matrix; the first column is the intercept. |
| `beta` | `J`-by-`K` coefficient matrix; the first row is zero. |
| `theta` | The `d = (J - 1) * K` free coefficients, stacked by outcome, with the intercept first in each block. |
| `eta` | `N`-by-`J` matrix of linear indices, with entries `x_i' beta_j`. |
| `p` | `N`-by-`J` matrix of conditional mean shares. |

Mullahy writes the mean share as ξ and uses the last outcome as the
baseline. Here `p` denotes that same mean, and `eta` names the linear
index. Changing the baseline changes the coefficient representation,
while fitted shares and partial effects remain the same.

Observed shares may be exactly zero or one. They can be used directly,
without replacing zeros by small positive numbers or taking logarithms of
the observed shares. For finite coefficients, fitted conditional means are
strictly inside the unit interval. Accepting observed zeros does not mean
the model estimates the probability that a share is exactly zero.

The model follows the fractional-response approach of Papke and Wooldridge,
extended to a system of shares by Mullahy. It specifies conditional means
without requiring a complete distribution for the shares. Its consistency
still relies on a correctly specified mean function and the usual
identification and sampling conditions.

## How to install fmlogit

Install the development version from the current repository:

```r
install.packages("remotes")
remotes::install_github("Fifis/fmlogit")
library(fmlogit)
```

### A small example

The following example generates shares for three activities, with a
continuous income index and a binary urban indicator. Simulated counts
provide a convenient way to construct the example shares; the estimator
also accepts shares without an underlying count denominator.

```r
set.seed(123)
N <- 300
x <- data.frame(income = rnorm(N), urban = rbinom(N, 1, 0.4))

eta <- cbind(
  home = 0,
  travel = 0.4 + 0.8 * x$income - 0.3 * x$urban,
  other = -0.2 - 0.5 * x$income + 0.6 * x$urban
)
p <- exp(eta)
p <- p / rowSums(p)

y <- t(vapply(seq_len(N), function(i) {
  as.numeric(rmultinom(1, size = 20, prob = p[i, ])) / 20
}, numeric(3)))
colnames(y) <- colnames(p)

weights <- ifelse(x$urban == 1, 2, 1)
fit <- fmlogit(y, x, weights = weights, print.level = 0)
summary(fit)
head(fitted(fit))
```

The rows of `y` are observations and its columns are outcome shares.
The rows of `x` describe the same observations. Here `home` is the
baseline because it is the first column of `y`.
The illustrative weights give urban observations twice the weight of other
observations. In a survey application, supply the appropriate sampling
weights for the population you want to describe.

## Why use fmlogit in R?

### 1. Integration with the R workflow

Estimation, prediction, effect summaries, and plots can be kept in the same
R analysis. Standard methods such as `summary()`, `fitted()`,
`residuals()`, `predict()`, and `vcov()` provide access to the results.
Plots use base R graphics and can return their underlying data for further
analysis.

### 2. Post-estimation on the share scale

A coefficient table describes relative-share parameters. Partial effects
translate those parameters into changes in expected shares, which are often
easier to interpret. You can select the covariates of interest, distinguish
derivatives from discrete changes, evaluate effects at a chosen profile,
or average them over the observed sample.

Uncertainty for these effects is available through joint coefficient
simulation. This supports summaries of effect estimates, standard errors,
z statistics, and p-values. Plotting individual effects can also reveal
variation that a single average hides.

### 3. Estimation flexibility

The covariate input can be a matrix or a data frame. Data-frame factors and
character variables are converted to design columns during fitting.
An intercept is added as the first design column, and supplied constant
columns are removed. Other redundant columns cause a rank-deficiency
error, so the design still needs to identify the model.

Optional nonnegative observation weights enter estimation, its robust
covariance calculation, and the mean profiles and average effects used
after fitting. Optimisation is handled by `maxLik`, with several methods
available.

### 4. Aggregates with user-supplied values

Sometimes a change in allocation is useful to express as a weighted total.
For example, categories may have different costs, returns, or external
valuations. The `wtp()` function multiplies their share effects by values
you supply and adds the results. The interpretation and units depend on
those supplied values. Their standard errors retain covariance between
the outcome effects being combined.

## How does the estimator work?

### Step 1. Specify the joint conditional means

The logit mean function respects both the unit-interval bounds and the
requirement that shares add up to one. Estimating separate fractional-logit
equations for each category would not generally impose that common total.

The estimation objective adds up each observed share multiplied by the log
of its fitted share, across categories and observations. Observation weights
multiply these contributions when supplied. This is a multinomial
**quasi-log-likelihood**: the observed shares enter where a categorical
multinomial logit model would use zero-or-one outcome indicators.

Consequently, fitting the model does not require the observed allocation
to be a single categorical choice. It also does not require shares to be
independent of one another. Their adding-up constraint already links them.

### Step 2. Maximise the quasi-log-likelihood

Pairwise fractional-logit fits against the baseline supply starting values;
zero coefficients are used if they give a better initial objective. The
estimator then optimises the full joint objective. With `J` outcomes and `K`
design columns including the intercept, there are `(J - 1) * K` free
coefficients.

The default optimiser is `"BHHH"`. Other supported `maxLik` methods,
such as `"BFGS"` or `"CG"`, can be selected through `method`.
Inspect `fit$convergence` and `fit$conv_code` when assessing the result.
Runtime depends on sample size, the number of coefficients, scaling, and
convergence; there is no general timing guarantee.

### Step 3. Compute robust standard errors

The package computes an HC0 sandwich covariance matrix for the full joint
coefficient vector, using observation-level scores and the weighted
Hessian. This retains covariance between coefficients for different
outcomes, which matters for joint inference.

The Hessian uses each observation's weight once, while the score outer
products use its square. Multiplying every weight by the same positive
constant leaves the estimates and HC0 covariance unchanged. These are
relative observation weights, rather than counts of replicated records.
The optimizer uses weights scaled to mean one, so their absolute units
do not alter its step sizes. Stored weights and the reported
quasi-log-likelihood retain the input scale.

Mullahy explains why ordinary multinomial-logit inverse-Hessian standard
errors are generally unsuitable for fractional shares. Shares can be less
variable than categorical outcomes with the same conditional means, so
robust standard errors need not be larger than conventional ones.

Use `vcov(fit)` to obtain the full joint covariance matrix. The sandwich
calculation allows the conditional variance to differ from the nominal
multinomial variance, while the conditional-mean specification remains
an important assumption.

## How do the post-estimation commands work?

### Reading coefficients

For a covariate entering linearly, its coefficient for a non-baseline
category describes the change per unit in the log ratio of that category’s
expected share to the baseline category’s expected share. The intercept
gives this log ratio when the other design columns are zero. These
coefficients do not directly give percentage-point changes in shares.

In particular, a positive coefficient does not guarantee a positive
marginal effect on that category’s expected share. Other categories also
respond to the covariate and enter the common denominator. Coefficient
tests are useful for hypotheses about these relative-share parameters;
hypotheses about share changes require inference for the corresponding
partial effects.

Changing the reference category changes the coefficient representation,
while fitted shares and their partial effects are invariant to that
normalisation.

### Marginal and discrete effects

For a continuous design variable, the marginal effect is the local slope
of the fitted share:

$$
\frac{\partial p_{ij}}{\partial x_{ik}}
= p_{ij}\left(\beta_{jk} -
  \sum_{h=1}^{J} p_{ih}\beta_{hk}\right).
$$

The formula shows why the effect depends on the whole model and the
covariate profile, rather than on one coefficient alone. An effect of
`0.02` means a local increase of two percentage points per unit of the
covariate. It is a slope; a large finite change may have a different
average slope.

A **discrete effect** compares fitted shares at two covariate values. The
package changes the selected design column from its minimum to its
maximum among observations with positive weight, while holding other
columns fixed. For a zero-or-one variable with both values represented
at positive weight, this is the change from zero to one.

For either kind of effect, the changes across all outcome categories sum
to zero. A gain in one share must be offset by losses elsewhere in the
allocation.

### Effects at the mean and average effects

Two summaries answer different questions:

| Option | Calculation | Interpretation |
| --- | --- | --- |
| `marg.type = "atmean"` | Use the weighted mean design profile, replacing the selected column by its endpoints for a discrete effect. | The effect at one reference profile. |
| `marg.type = "aveacr"` | Evaluate effects for each retained observation, then average using the fitted model's observation weights. | The average effect for the distribution represented by those weights. |

With normalized weights $\alpha_i = w_i / \sum_h w_h$, average effects
are $\sum_i \alpha_i m_{ijk}$, where $m_{ijk}$ is an individual effect.
The default reference profile is $\bar x_w = \sum_i \alpha_i x_i$.
Equal weights give the usual arithmetic means. Zero-weight observations
do not contribute to these summaries or determine automatic endpoints.
This weighted averaging follows Mullahy's Appendix 1; whether the result
represents a population depends on the meaning and validity of the weights.

The default is `"atmean"`. In a nonlinear model, the effect at the mean
generally differs from the mean of the effects. A mean profile can also
contain fractional values of binary variables and need not describe an
actual observation. Wulff discusses why both summaries can be informative
and why either can hide important variation.

For example, calculate average marginal effects of the income index and
average discrete changes in the urban indicator:

```r
set.seed(123)
ame <- effects.fmlogit(
  fit, effect = "marginal", marg.type = "aveacr",
  varlist = "income", se = TRUE, R = 1000
)
summary(ame)

urban_change <- effects.fmlogit(
  fit, effect = "discrete", marg.type = "aveacr", varlist = "urban",
  se = TRUE, R = 1000
)
urban_change$effects
```

With `marg.type = "atmean"`, `at` supplies a custom reference profile:

```r
at <- colSums(fit$x[, -1, drop = FALSE] * fit$weights) / sum(fit$weights)
at["income"] <- 0
effects.fmlogit(fit, varlist = "income", at = at)$effects
```

The values in `at` must follow the non-intercept column order of
`fit$x`. An explicit profile is used as supplied; it is not averaged
again. For multi-level factors or interactions, a change in one design
column may not represent a valid change in the original variable. Use
compatible covariate profiles and update all related design columns when
constructing such comparisons. Derivatives of a raw variable that also
enters through powers or interactions require accounting for those terms.

### Predictions and plots

Predictions describe expected share levels; marginal effects describe
their slopes. Looking at both helps avoid confusing a large share with a
large effect, or a declining positive slope with a declining share.

```r
newdata <- data.frame(income = c(-0.5, 0.5), urban = c(0, 1))
predict(fit, newdata = newdata)
head(residuals(fit))
```

For new observations, `predict()` expects the processed, numeric design
columns in the order of `colnames(fit$x)[-1]`, with the intercept
excluded. Factors must already be encoded into the same design columns.
The simple example above has only numeric covariates.

Following Wulff’s emphasis on graphical interpretation, inspect how
effects vary across observations or meaningful covariate values. The
margins plot uses the observation-level matrices saved by
`marg.type = "aveacr"`:

```r
plot(ame, x_data = fit$x, against.x = "income", group.x = "urban")

plot_data <- plot(
  ame, x_data = fit$x, against.x = "income", plot.show = FALSE
)
head(plot_data$effects)
```

Both plot methods return axis values, an effect matrix, an optional SE
matrix, and grouping information. Observation-level margins plots have
no individual SEs: the SEs in `ame` describe its average effects. Lines
through effects at chosen profiles and plots of individual effects also
answer different questions; choose the display to match the question.
Automatic profile curves use the weighted means for the other covariates
and a grid spanning the positive-weight observations. Fitted shares,
residuals, individual effects, and plotted points remain observation-level
values; the package does not multiply them by their observation weights.

### Standard errors for partial effects

With `se = TRUE`, the package uses the simulation-based **Krinsky–Robb
method**. It draws the entire non-baseline coefficient vector from an
estimated joint normal distribution using `vcov(fit)` and recomputes the
requested effects for each draw using the same fixed observation weights.
The ordinary sample standard deviation across these equally likely draws
is the SE. Observation weights apply within each effect calculation;
they are not weights for the simulation draws.

Keeping the coefficients together preserves their covariance across
outcomes and covariates. `ame$cov.list` stores a `J`-by-`J` covariance
matrix of the outcome effects for each selected covariate, and
`ame$weights` retains the observation weights used for averaging.

These SEs describe coefficient uncertainty while holding the observed
covariates and weights, or the chosen reference profile, fixed. They do
not describe the spread of individual effects across the sample. They
rely on an asymptotic normal approximation. Set a seed for reproducibility and
increase `R` if Monte Carlo variation is material for the application.

The reported z statistics divide each effect estimate by its simulated
SE; p-values use the standard normal distribution. Approximate 95% Wald
intervals can be formed as the estimate plus or minus 1.96 SEs. These are
different from taking percentiles of the simulated effects.

The delta method also supports SEs, confidence intervals, and hypothesis
tests; Wulff uses it in his examples. Mullahy illustrates bootstrap
inference for average partial effects. These are alternative ways to
assess uncertainty. The package currently implements Krinsky–Robb SEs for
partial effects.

### Weighted aggregates of effects

Interpret `wtp()` as a user-valued aggregate of share changes. For
example, assign external values of 0, 10, and 5 to the three categories:

```r
values <- c(home = 0, travel = 10, other = 5)
valued_change <- wtp(urban_change, wtp.vec = values, indv.obs = TRUE)
valued_change$wtp

curve <- plot(
  fit, wtp.vec = values, varlist = "income", against = "income",
  effect = "marginal", t = 40
)
```

Values are used in outcome-column order. The aggregate estimates describe
changes in the supplied value scale, and `wtp.obs` contains the
observation-level changes when requested from average-across margins.
Observation weights have already entered the average effects; the supplied
category values are a separate aggregation across outcomes.

When effect SEs are available, `wtp()` uses the joint outcome covariance
from the same Krinsky–Robb draws: for values $v$ and effect covariance
$C_k$, the aggregate variance is $v'C_kv$. Optional plot bars use the same
calculation. This retains the covariance required by the adding-up
constraint. If all categories have the same value, reallocating shares
changes neither the total value nor its estimated uncertainty, apart
from numerical rounding.

## Practical concerns

**Prepare a meaningful share system.** Supply all mutually exclusive and
exhaustive categories on a zero-to-one scale, with each row summing to one.
Convert percentages before fitting. Shares based on a total of zero are
undefined; decisions about excluding those observations affect the
population to which results apply. Ye and Pendyala’s time-use application
and Mullahy’s portfolio example illustrate why defining the total and the
categories matters.

**Check the estimation sample and design.** Rows with missing responses
or covariates are removed together, and supplied weights are subset using
the same rows. `fit$rowNo` records the retained original positions.
Consider rescaling covariates when their numerical scales differ greatly,
and investigate redundant columns or poorly identified outcomes when
optimisation or covariance estimation fails.

**Keep the weighting target consistent.** The optional `weights` argument
supplies relative observation weights that are retained after fitting and
inherited by average effects, default mean profiles, and profile plots.
Multiplying all weights by a common positive constant changes none of
these summaries. Appropriate survey weights can make the represented
distribution differ from the raw sample; arbitrary weights do not by
themselves establish a population interpretation.

**Assess the mean specification and the sampling assumptions.** Robust
SEs do not repair a poorly specified conditional mean. Logit imposes a
particular relationship between relative shares and the covariates, and
sample averages can conceal heterogeneous effects. Examine predictions
and effects over meaningful profiles, avoiding combinations with little
support in the data. Conditional associations acquire a causal
interpretation only with an appropriate research design and assumptions.

**Match the model to the question.** A model for expected allocation does
not describe the entire distribution of allocations or separately model
the decision to participate in an activity. Ramalho, Ramalho, and
Murteira discuss cases where a separate participation and allocation
model can be useful. Standard beta and Dirichlet models offer other
approaches, usually with stronger distributional assumptions and, in
their basic forms, restrictions on boundary observations. Those models,
two-part specifications, and the specification-test procedures studied
in the papers are not implemented here.

**Use appropriate inference for dependent observations.** The current
sandwich estimator treats observations as independent. Observation
weights alone do not account for survey clusters, strata, serial
correlation, or repeated observations on the same unit. Population-average
effects may also require accounting for uncertainty in the covariate
distribution, beyond the fixed-design simulation provided here.
Ordinary likelihood-ratio procedures for fully specified probability
models should not be applied automatically to the quasi-likelihood.

## References

- **Mullahy, J. (2015).** [Multivariate fractional regression estimation of
  econometric share models](https://doi.org/10.1515/jem-2012-0006).
  *Journal of Econometric Methods*, 4(1), 71–100.
  The main reference for the multivariate mean specification, robust
  inference, and share-scale partial effects.

- **Wulff, J. N. (2014).**
  [Interpreting results from the multinomial logit model: Demonstrated by
  foreign market entry](https://doi.org/10.1177/1094428114560024).
  *Organizational Research Methods*, 18(2), 300–325.
  Guidance on coefficients, marginal effects, predicted outcomes, and plots.

- **Papke, L. E., and Wooldridge, J. M. (1996).** [Econometric methods for
  fractional response variables with an application to 401(k) plan
  participation rates.](https://doi.org/10.1002/(SICI)1099-1255(199611)11:6%3C619::AID-JAE418%3E3.0.CO;2-1)
  *Journal of Applied Econometrics*, 11(6), 619–632.
  The foundational univariate fractional-response approach.

- **Ramalho, E. A., Ramalho, J. J. S., and Murteira, J. M. R. (2011).**
  [Alternative estimating and testing empirical strategies for fractional
  regression models](https://doi.org/10.1111/j.1467-6419.2009.00602.x).
  *Journal of Economic Surveys*, 25(1), 19–68.
  An overview of univariate fractional-regression specifications and tests.

- **Ye, X., and Pendyala, R. M. (2005).** A model of daily time use allocation
  using fractional logit methodology. In H. S. Mahmassani (Ed.),
  *[Transportation and Traffic Theory: Flow, Dynamics and Human Interaction](https://shop.elsevier.com/books/transportation-and-traffic-theory/mahmassani/978-0-08-044680-6)*,
  pp. 507–524. Elsevier.
  An application to daily time allocation.

The original R implementation was developed by James Xinde Ji. Current
development and maintenance are by Andreï V. Kostyrka; see `DESCRIPTION`
for package authorship details.
