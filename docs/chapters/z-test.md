# Univariate Two-Sample Test

We collect data during or after experiments and the data are used for statistical analysis. A problem of primary interest is to figure out whether the average treatment effect is significant. In this chapter, we introduce popular two-sample tests for this problem.

## Notation and Setup



## Null Hypotesis of Interest

Denote $\delta=\mu_1-\mu_0$. We are interest in testing the null hypothesis $H_0: \delta=0$ versus the alternative hypothesis $H_1: \delta\neq 0$.

## Test Statistic

A test statistic serves as the key to generate a testing procedure. To test the null, we consider the t-statistic 
$$
T = \frac{\bar{y}_1 - \bar{y}_0}{\sqrt{s_1^2/n_1 + s_0^2/n_0}}.
$$
This test statistic is motivated from that $\bar{y}_1 - \bar{y}_0$ is a nutural sample estimator of $\delta$ and its variance is given by
$$
\operatorname{Var}(\bar{y}_1 - \bar{y}_0)=\frac{\sigma_1^2}{n_1} + \frac{\sigma_0^2}{n_0}.
$$
 It can be shown that $s_d^2$ are unbiased and consistent sample estimators of $\sigma^2_d$. By the Central Limit Theorem (CLT) and Slutsky's Theorem , one can verify that $T$ is asymptotically standard normal under $H_0$.

## Critical Region

A critical region, also called rejection region, is a set of values for the sample data, constructed based on the test statistic. The critical region identifies the testing procedure: we reject the null when the observed data fall into the critical region, and have no evidence to reject the null when the observed data do not fall into the critical region.

In this example, suppose we want to control the type 1 error (the chance of rejecting the null when it is actually correct) at level $\alpha$, a two-sided critical region for $T$ can be given as 
$$
R_{\alpha}=(-\infty, -z_{\alpha/2})\cup(z_{\alpha/2},+\infty)
$$
where $z_{\alpha}$ is the upper quantile of the standard normal distribution. Because $T$ is asymptotically standard normal, we observe that the type 1 error of the derived testing procedure satisfies
$$
P(T\in R_{\alpha}\mid H_0) \rightarrow 2(1-\Phi(z_{\alpha/2}))=\alpha.
$$

## P-Value
An alternative way of developing a testing procedure is to compute the so-called p-value, which provides information about the strength of evidence against the null hypothesis. P-value, denoted by $p=p(T)$, is itself a test statistic and is constructed in a way that smaller p-value implies stronger evidence against the null. Therefore, the null is rejected when $p\leq\alpha$. Given a significance level $\alpha$, a p-value is valid if and only if for every $\alpha\in[0,1]$,
$$
P(p(T)\leq \alpha\mid H_0)\leq \alpha.
$$
In other words, the resulting testing procedure controls type 1 error at $\alpha$.
## Pivitol Quantity

