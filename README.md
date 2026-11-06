# Willingness-to-Pay for Mixed Logit Models

Simulation-based estimation of willingness-to-pay (WTP) and its uncertainty for discrete choice models estimated with [Apollo](https://www.apollochoicemodelling.com/) in R. The tool handles WTP measures where the numerator and/or denominator include **random coefficients** (normal or uniform), **interactions with individual characteristics**, and **non-linear transformations** such as lognormal cost coefficients. Bootstraping is usually recommended, however it very computationally heavy while this Monte-Carlo simulation only takes seconds. 

## Motivation

WTP is a ratio of coefficients, for example the value of travel time:

$$\text{WTP} = \frac{\beta_{time}}{\beta_{cost}}$$

In a simple MNL model, a standard error can be obtained with the delta method. In mixed logit models this becomes harder:

- the numerator and denominator may themselves be random across individuals,
- the coefficients may be transformed (e.g. $\beta_{cost} = -\exp(\mu + \sigma \eta)$),
- WTP may depend on individual characteristics through interactions,
- ratios of random variables are often skewed or heavy-tailed, so a symmetric interval can be misleading.

This function addresses this by simulation: it propagates both estimation uncertainty and, where relevant, preference heterogeneity through the full ratio, and summarises the resulting distribution.

## Method

The function works in five steps:

1. **Parameter uncertainty (Krinsky–Robb).** Draw `R` parameter vectors from a multivariate normal distribution with mean equal to the estimated coefficients and covariance equal to the (robust) covariance matrix:
   $\theta^{(r)} \sim N(\hat\theta, \hat\Sigma)$, $r = 1, \dots, R$.

2. **Preference heterogeneity.** If the WTP involves random coefficients, each parameter draw is repeated `K` times, and an individual-level coefficient is drawn for each:
   - *Normal:* $\beta^{(r,k)} \sim N(\mu^{(r)}, \sigma^{(r)})$. If a second parameter $\rho$ is supplied (an off-diagonal Cholesky element), the standard deviation becomes $\sqrt{\sigma^2 + \rho^2}$.
   - *Uniform:* $\beta^{(r,k)} \sim U(\mu^{(r)} - |s^{(r)}|,\ \mu^{(r)} + |s^{(r)}|)$.

   This gives `R × K` draws in total.

3. **Transformations.** User-supplied functions can be applied to individual parameters after drawing (e.g. `-exp()` for a negative lognormal coefficient), and to the summed numerator and denominator (e.g. unit conversion).

4. **Interactions.** If parameters interact with individual characteristics, the numerator and denominator are computed per individual as $\sum_j \alpha_j x_{ij}$, using one row per respondent from the data. Within each draw, WTP is averaged across individuals.

5. **Summary.** WTP is computed as numerator / denominator for every draw, and the function returns the mean, median and the 2.5% and 97.5% quantiles.

## What the interval represents

This matters for interpretation:

| Model specification | Number of values | Interpretation of the 2.5–97.5% interval |
|---|---|---|
| Fixed coefficients, no interactions | `R` | Confidence interval for WTP (estimation uncertainty only) |
| Fixed coefficients with interactions | `R` | Confidence interval for **sample-average** WTP |
| Random coefficients | `R × K` | Interval of the **WTP distribution in the population**, including both preference heterogeneity and estimation uncertainty. It is *not* a confidence interval for mean WTP |

When the denominator is a random coefficient whose distribution includes zero (e.g. a normal cost coefficient), the ratio has no finite mean, and `wtp_mean` will be unstable. In that case, report the median, or constrain the sign of the denominator, e.g. with a lognormal specification via `trans_param`.

## Installation

The function is currently provided as a script:

```r
source("WTP_random_coef_function.R")
```

Requires the `MASS` package. Models are assumed to be estimated with `apollo`.

## Usage

### Example 1: Value of travel time with random coefficients

Travel time coefficient normally distributed; cost coefficient negative lognormal. Travel time is measured in minutes, so the result is converted to per hour.

```r
library(apollo)
source("WTP_random_coef_function.R")

# ... standard Apollo workflow ...
model <- apollo_estimate(apollo_beta, apollo_fixed,
                         apollo_probabilities, apollo_inputs)

# Transformation functions are passed by name
neg_exp <- function(x) -exp(x)
per_hour <- function(x) x * 60

vtt <- apollo_WTP(
  alpha_name          = "b_tt",
  random_alpha_normal = list(b_tt = "s_tt"),
  beta_name           = "b_cost_log",
  random_beta_normal  = list(b_cost_log = "s_cost_log"),
  trans_param         = list(neg_exp = "b_cost_log"),
  trans_alpha         = "per_hour",
  estimated_coef      = model$estimate[!names(model$estimate) %in% apollo_fixed],
  estimated_Sigma     = model$robvarcov,
  database            = database,
  R = 1000,
  K = 500
)

vtt$wtp_median
vtt$wtp_95_confidence
```

### Example 2: WTP varying with income (interaction in the denominator)

```r
wtp_income <- apollo_WTP(
  alpha_name       = "b_tt",
  beta_name        = c("b_cost", "b_cost_income"),
  interaction_beta = list(b_cost_income = "income"),
  estimated_coef   = model$estimate[!names(model$estimate) %in% apollo_fixed],
  estimated_Sigma  = model$robvarcov,
  database         = database,
  R = 1000
)
```

Here the denominator for individual *i* is `b_cost + b_cost_income * income_i`. No random coefficients are involved, so `K` is not needed, and the interval is a confidence interval for the sample-average WTP.

## Arguments

| Argument | Description |
|---|---|
| `alpha_name` | Character vector of all parameters in the **numerator**. These are summed (multiplied by their interaction attribute, if any). |
| `random_alpha_normal` | Named list of normally distributed numerator parameters and their standard-deviation parameter(s): `list(b_tt = "s_tt")` or `list(b_tt = c("s_tt", "rho_tt"))`. |
| `random_alpha_uniform` | Named list of uniformly distributed numerator parameters and their spread parameter: `list(b_x = "s_x")`. |
| `interaction_alpha` | Named list mapping numerator parameters to a column in `database`: `list(b_tt_income = "income")`. |
| `beta_name` | Character vector of all parameters in the **denominator**. |
| `random_beta_normal`, `random_beta_uniform`, `interaction_beta` | As above, for the denominator. |
| `trans_param` | Named list applying a function to individual parameters after drawing: `list(neg_exp = c("b_cost_log"))`. Names are function names. |
| `trans_alpha`, `trans_beta` | Character vector of function names applied to the summed numerator/denominator, in the order given. |
| `estimated_coef` | Named numeric vector of estimated coefficients, **excluding fixed parameters**. |
| `estimated_Sigma` | Covariance matrix with row and column names matching `estimated_coef`. The robust covariance matrix (`model$robvarcov`) is recommended. |
| `database` | The estimation data. Only needed when interactions are used. |
| `R` | Number of draws from the parameter distribution. |
| `K` | Number of individual-level draws per parameter draw. Required when random coefficients are included. |

## Value

A list with:

- `wtp_dist()`: function returning the vector of simulated WTP values (stored with class `"hidden_vector"` so the full vector is not printed by accident),
- `wtp_mean`: mean of the simulated values,
- `wtp_median`: median of the simulated values,
- `wtp_95_confidence`: 2.5% and 97.5% quantiles.

## Input validation

The function checks that:

- `estimated_coef` is a named numeric vector, and all names appear in `estimated_Sigma`,
- `estimated_Sigma` is a square matrix with matching row and column names,
- all parameters in `alpha_name` and `beta_name` exist among the estimated coefficients,
- all standard-deviation parameters referenced for random coefficients exist,
- all interaction attributes exist as columns in `database`,
- `K` is supplied when random coefficients are used.

## Practical notes

- **Memory.** With random coefficients, the draw matrix has `R × K` rows. With interactions as well, an `(R × K) × N` matrix is created, where `N` is the number of respondents. Reduce `R` or `K` for large samples.
- **Number of draws.** Check stability by rerunning with a different seed or more draws. Quantiles converge more slowly than the mean.
- **Transformation functions** are looked up by name with `get()`, so they must exist in the environment where the function is called.

## Known limitations

- The respondent identifier column is currently hard-coded as `Respondent_Serial`.
- When a correlation parameter is supplied, it only affects the marginal standard deviation. Individual-level draws for different coefficients are drawn independently, so correlation between numerator and denominator is not carried over.
- Only normal and uniform distributions are supported directly. Other distributions (e.g. lognormal) can be obtained via `trans_param`.

## Other files

- `WTP_simple_function.R`: *[TODO: short description]*
- `test_function.R`: *[TODO: short description of tests/examples]*

## References

- Krinsky, I. & Robb, A. L. (1986). On approximating the statistical properties of elasticities. *Review of Economics and Statistics*, 68(4), 715–719.
- Daly, A., Hess, S. & de Jong, G. (2012). Calculating errors for measures derived from choice modelling estimates. *Transportation Research Part B*, 46(2), 333–341.
- Train, K. (2009). *Discrete Choice Methods with Simulation* (2nd ed.). Cambridge University Press.
- Hess, S. & Palma, D. (2019). Apollo: A flexible, powerful and customisable freeware package for choice model estimation and application. *Journal of Choice Modelling*, 32, 100170.

## Author

Simun Pauli, DTU
