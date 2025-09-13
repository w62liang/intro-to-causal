---
include: math
---

# Regression Analysis of Experimental Data

In the stage of statistical analysis of experimental data, baseline covariates are often available and highly encouraged to be incorporated into analysis to enhance the precision of statistical inference and the credibility of the results. Let $X_i$ denote baseline covariates, collected prior to treatment assignment. In this chapter, we introduce three commonly used regression models for $\mathcal{G}=\\{(Y_i,A_i,X_i),i=1,\cdots,n\\}$, corresponding regression adjustment approaches, and interpret the regression coefficients. We first consider two-arm randomized experiments where $A_i$ is an indicator variable of treatment assignment. 

## Two-Arm Experiments under Covariate-Free Randomization

### Analysis of Variance (ANOVA)

Recall that $\tau=E[Y_i(1)-Y_i(0)]$ is the ATE and $\mu_0=E[Y_i(0)]$. Under covariate-free randomization, one can verify that 
$$
E(Y_i- \tau A_i - \mu_0\mid A_i)=0,
$$
which yields the working model for ANOVA:
$$
E(Y_i\mid A_i)=\tau A_i + \mu_0.
$$
The ANOVA model is correct due to randomization. To estimate $\tau$ using the sample data, we can apply the ordinary least squares (OLS) method and solve
$$
\arg\min_{\tau,\mu_0} \sum_{i=1}^n(Y_i-\tau A_i - \mu_0)^2,
$$
which yields consistent estimators $\hat{\mu}\_0=\bar{Y}\_0$ and $\hat{\tau}\_{AN}=\bar{Y}\_1 - \bar{Y}\_0$ for $\mu_0$ and $\tau$. However, the variance estimator obtained from OLS is slightly different and greater than the one we provided in the former chapter.

### Analysis of Covariance (ANCOVA)

ANCOVA incorporates the baseline covariates $X$ into its working model by assuming
$$
E(Y_i\mid A_i,X_i)=\tau A_i + \beta^\top (X_i-\mu_X) + \beta_0
$$
for some $\beta$ and $\beta_0$. Under covariate-free randomization, simple algebra gives that
$$
\tau(X_i)=E[Y_i(1)-Y_i(0)\mid X_i]=E(Y_i\mid A_i=1,X_i) - E(Y_i\mid A_i=0,X_i) = \tau,
$$
implying the ANCOVA model is incorrect for unless $\tau(X_i)=\tau$, i.e., the ATE is a constant across different subgroups. Let us impose this constant conditional ATE assumption. Then, $\beta_0=\mu_0$. The standard least squares method solves
$$
\arg\min_{\tau,\beta_0,\beta} \sum_{i=1}^n\\{Y_i-\tau A_i - \beta^\top (X_i-\bar{X}) - \beta_0\\}^2
$$
where $\bar{X}=\sum_{i=1}^nX_i/n$ and yields the estimator for $\tau$ as
$$
\hat{\tau}\_{ANC} = \bar{Y}\_1 - \bar{Y}\_0 - \hat{\beta}^\top(\bar{X}\_1 - \bar{X}\_0) 
$$
where 
$$
\hat{\beta} = \left\\{\sum_{i=1}^n(X_i-\bar{X})(X_i-\bar{X})^\top\right\\}^{-1}\sum_{i=1}^n(X_i-\bar{X})Y_i+ o_p(n^{-1/2})
$$
is the OLS estimator for $\beta$. The $\hat{\tau}\_{ANC}$ is an unbiased and consistent estimator of $\tau$. Despite having incoporated the information of baseline covariates, the asymptotic variance of $\hat{\tau}\_{ANC}$ is not necessarily smaller than that of $\hat{\tau}\_{AN}$ unless $n_1=n_0$ or that $\beta_1=\beta_0$ where $\beta_a = \operatorname{cov}(X)^{-1}\operatorname{cov}(X,Y(a))$ for $a=0,1$. In next section, we introduce a regression model which produces an estimator of $\tau$ that is unabised, consistent, and guaranteed to be more efficient than $\hat{\tau}\_{AN}$ and also $\hat{\tau}\_{ANC}$.

### Analysis of Heterogeneous Covariance (ANHECOVA)

ANHECOVA considers the following model:
$$
E(Y_i\mid A_i,X_i)=\tau A_i + \beta^\top (X_i-\mu_X) + \theta^\top A_i(X_i-\mu_X)  + \beta_0,
$$
which, compared to the ANCOVA model, adds the treatment-by-covariate interaction terms. Simple algebra reveals that $\beta_0 + \beta^\top (X_i-\mu_X) = E(Y_i\mid A_i=0,X_i)$, $\beta_0=\mu_0$,
$$
\tau(X_i) = E(Y_i\mid A_i=1,X_i) - E(Y_i\mid A_i=0,X_i) = \tau + \theta^\top (X_i-\mu_X),
$$
and $\tau = E[Y_i(1)-Y_i(0)]$. Therefore, the ANHECOVA model is correct for $\tau$, but incorrect for $\tau(X_i)$ unless $\tau(X_i)$ is truly linear in $X_i$, in which case $\theta$ can be interpreted as the ATE conditional on the covariates $X$. When $\tau(X_i)$ is not linear in $X_i$, $\theta^\top (X_i-\mu_X)$ is the projection of $\tau(X_i)$ to the space spanned by $X_i$. We will explain this later.

The ANHECOVA model has been studied for decades, but the name of "ANHECOVA" was firstly introduced by Ye et al. (2023), who also proved that the ANHECOVA model is applicable to a general class of randomization methods and and the estimator of $\tau$ produced by ANHECOVA is "optimal" among a class of linear-adjusted estimators. Therefore, ANHECOVA is preferred compared to ANCOVA especially when the sample size is sufficiently large or the group sizes are not balanced.

To estimate $\tau$ and $\theta$, we can apply the OLS and solve
$$
\arg\min_{\tau,\beta_0,\beta,\theta} \sum_{i=1}^n\\{Y_i-\tau A_i - \beta^\top (X_i-\bar{X}) - \theta^\top A_i(X_i-\bar{X})  - \beta_0\\}^2,
$$
leading to the estimator
$$
\hat{\tau}\_{ANHC} = \bar{Y}\_1 - \bar{Y}\_0 -(\hat{\beta}+\hat{\theta})^\top (\bar{X}\_1-\bar{X}) + \hat{\beta}^\top (\bar{X}\_0-\bar{X}) + o_p(n^{-1/2})
$$
where
$$
\hat{\beta} = \left\\{\frac{1}{n}\sum_{i=1}^n(X_i-\bar{X})(X_i-\bar{X})^\top\right\\}^{-1}\frac{1}{n_0}\sum_{A_i=0}(X_i-\bar{X})Y_i
$$
and 
$$
\hat{\theta} = \left\\{\frac{1}{n}\sum_{i=1}^n(X_i-\bar{X})(X_i-\bar{X})^\top\right\\}^{-1}\sum_{i=1}^n(X_i-\bar{X})\left\\{\frac{A_i}{n_1}Y_i - \frac{(1-A_i)}{n_0}Y_i  \right\\},
$$
which converges to $\operatorname{cov}^{-1}(X)\operatorname{cov}(X,Y(1)-Y(0))$, i.e., the regression coefficient of $Y(1)-Y(0)$ on $X$. We see that $\hat{\tau}\_{ANHC}=\hat{\tau}\_{AN}-r_n+o_p(n^{-1/2})$ where $r_n=(\hat{\beta}+\hat{\theta})^\top (\bar{X}\_1-\bar{X}) - \hat{\beta}^\top (\bar{X}\_0-\bar{X})$ is the adjusted term. It can be shown that $\operatorname{cov}(\hat{\tau}\_{ANHC},r_n)=o_p(1)$, implying 
$$
\begin{aligned}
\operatorname{var}(\hat{\tau}\_{ANHC}) &= \operatorname{var}(\hat{\tau}\_{AN}-r_n) + o_p(n^{-1})\\\\
&=\operatorname{var}(\hat{\tau}\_{AN}) -\operatorname{var}(r_n) + o_p(n^{-1}) \leq \operatorname{var}(\hat{\tau}\_{AN}) + o_p(n^{-1}).
\end{aligned}
$$
Therefore, ANHECOVA guarantees effciciency gains from covariate adjustment when $\tau(X_i)$ is somewhat correlated with $X_i$. In fact, $\operatorname{var}(\hat{\tau}\_{ANHC})$ has smallest asymptotic variance among the estimators of the form 
$$
\\{\bar{Y}\_1 - \bar{Y}\_0 -\beta_1^\top (\bar{X}\_1-\bar{X}) + \beta_0^\top (\bar{X}\_0-\bar{X}), \beta_1,\beta_2\in\mathbb{R}^p\\}.
$$

## Two-Arm Experiments under Covariate-Adaptive Randomization

We consider regression analysis under covariate-adaptive randomization, where the stratification variables are allowed to share same components with $X_i$. We assume, as we did in Chapter 3, the sample units are divided into $K$ strata indexed by $\\{1,\cdots,K\\}$. We use $Z_i$ to denote a $K$-dimensional indicator variable of strata and let $Z_i=e_k$ if the unit $i$ is assigned to stratum $k$, where $e_k$ is a $K$-dimentional one-hot vector whose $k$-th component is $1$ and others are $0$. 