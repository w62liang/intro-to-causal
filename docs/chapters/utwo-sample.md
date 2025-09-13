---
include: math
---

# Univariate Two-Sample Test

Practitioners collect data during or at the end of an experiment and the data are used for statistical analysis. When the outcome measurement is a scalar variable at a certain time point, a problem of primary interest is whether the average treatment effect (ATE) on the outcome is significant. In this chapter, we introduce two-sample tests for this problem and discuss related issues such as sample size determination.

## Statistical Models and Notation

Suppse we collect data $\mathcal{G}=\\{(Y_i,A_i,X_i),i=1,\cdots,n\\}$ through an experiment where $A_i\in\\{0,1\\}$ are treatment indicators, $Y_i=Y_i(A_i)$ are observed outcomes (survival; cell concentration), and $X_i$ are baseline covariates (age, gender, blood pressure, etc.). We assume

 * **A1**. $(Y_i(1),Y_i(0))$ are i.i.d. copies of $(Y(1),Y(0))$.
 * **A2**. $\\{(Y_i(1),Y_i(0)),i=1,\cdots,n\\}\perp \\!\\!\\!\perp \\{A_1,\cdots,A_n\\}$.
 * **A3**. $E(A_i)=\pi\in(0,1)$ for $i=1,\cdots,n$.

These assumptions are mild and standard for experimental data under covariate-free randomization of treatment assignments. Denote $n_a=\\#\\{i:A_i=a\\}$, $\mu_a=E[Y(a)]$ and $\sigma^2_a=\operatorname{var}\\{Y(a)\\}$ for $a=0,1$. Let $\bar{Y}\_a=n_a^{-1}\sum_{i:A_i=a}Y_i$ be the sample mean in group $A_i=a$ and 
$$
S_a^2 = \frac{1}{n_a-1}\sum_{i:A_i=a}(Y_i-\bar{Y}\_a)^2
$$
the correponding sample variance. The two statistics are respectively unbiased estimators of $\mu_a$ and $\sigma_a^2$. The average treatment effect (ATE), also known as the risk difference (RD) in epidemiology, is defined by $\tau=\mu_1-\mu_0$. 

## Two-Sided Z-Test

We are primarily interest in testing the null hypothesis $H_0: \tau=0$ versus the alternative hypothesis $H_1: \tau\neq 0$.

### T-Statistic

To test the null, we consider the t-statistic $T = T(\mathcal{G})=\hat{\tau}/S\_{\tau}$ where $\hat{\tau}=\bar{Y}\_1 - \bar{Y}\_0$ and $S\_{\tau}=\sqrt{S_1^2/n_1 + S_0^2/n_0}$. This test statistic is motivated from that $\hat{\tau}$ is a natural sample estimator of $\tau$ and $S^2_{\tau}$ is the sample estimator for
$$
\sigma\_{\tau}^2=\operatorname{Var}(\hat{\tau})=\frac{1}{n}\left\\{\frac{\sigma_1^2}{\pi} + \frac{\sigma_0^2}{1-\pi}\right\\}.
$$
By the Central Limit Theorem (CLT) and the Slutsky's Theorem , one can verify that $T$ is asymptotically standard normal under $H_0$.

### The Testing Procedure

Intuitively, $T$ should be close to $0$ when $H_0$ is correct. Let $\mathcal{O}=\\{(y_i,a_i,x_i),i=1,\cdots,n\\}$ denote realized values of the data. The critical region therefore has the form $R(c)=\\{\mathcal{O}:|T(\mathcal{O})|>c \\}$ for some constant $c>0$.

To control Type I Error at $\alpha$, we hope $P(\mathcal{G}\in R(c)\mid H_0)=P(T>c\mid H_0)+P(T<-c\mid H_0)\leq \alpha$. Since $T\mid H_0\overset{d}{\rightarrow}N(0,1)$, we can set $c=z_{\alpha/2}$ as the upper $\alpha/2$-quantile of $N(0,1)$ (i.e., $P(N(0,1)\geq z_{\alpha/2})=\alpha/2$). It follows that
$$
P(\mathcal{G}\in R(c)\mid H_0)=P(T>z_{\alpha/2}\mid H_0)+P(T<-z_{\alpha/2}\mid H_0)\rightarrow \alpha.
$$

### Power and Sample Size Determination

The power of the testing procedure is given by
$$
\begin{aligned}
\beta(\tau,\sigma\_{\tau})=&P_{\tau}(|T|>z_{\alpha/2})\\\\
=& P_{\tau}\left(\frac{\hat{\tau}}{S\_{\tau}}>z_{\alpha/2}\right) + P_{\pi}\left(\frac{\hat{\tau}}{S\_{\tau}}<-z_{\alpha/2}\right)\\\\
=& P_{\tau}\left(\frac{(\hat{\tau}-\tau)}{S\_{\tau}}>z_{\alpha/2}-\frac{\tau}{S\_{\tau}}\right) + P_{\tau}\left(\frac{(\hat{\tau}-\tau)}{S\_{\tau}}<-z_{\alpha/2}-\frac{\tau}{S\_{\tau}}\right)\\\\
= & 1-\Phi\left(z_{\alpha/2}-\frac{\tau}{\sigma\_{\tau}}\right) + \Phi\left(-z_{\alpha/2}-\frac{\tau}{\sigma\_{\tau}}\right) + o(1).
\end{aligned}
$$
Taking partial derivative of the power function with respect to $\sigma\_{\tau}$ gives 
$$
\begin{aligned}
\frac{\partial \beta(\tau,\sigma\_{\tau})}{\partial \sigma\_{\tau}}&=\frac{\tau}{\sigma\_{\tau}^2} \left\\{\phi(z_{\alpha/2}+\sigma\_{\tau}^{-1}\tau) - \phi(z_{\alpha/2}-\sigma\_{\tau}^{-1}\tau)\right\\}+o(1)\\\\
&\approx \frac{|\tau|}{\sigma\_{\tau}^2} \left\\{\phi(z_{\alpha/2}+\sigma\_{\tau}^{-1}|\tau|) - \phi(z_{\alpha/2}-\sigma\_{\tau}^{-1}|\tau|)\right\\}\leq 0,
\end{aligned}
$$
and the equality holds if and only if $\tau=0$. This implies that, asymptotically, the Type II Error of the proposed testing procedure drops down as $n$ increases because $\sigma\_{\tau}$ decreases. When fixing $\tau>0$, $\beta(\tau,\sigma\_{\tau})=\gamma$ has a unique solution which is the standard deviation needed for detecting the ATE of at least $\tau$ with the power $\gamma$. However, the solution has no closed-form and can only be solved from the equation $\beta(\tau,\sigma\_{\tau})=\gamma$ using, for example, the Newton's method.

Let us consider the settings where the outcomes are binary. Then, we can verify that $S_a^2=[n_a/(n_a-1)]\bar{Y}_a(1-\bar{Y}_a)$ and
$$
S\_{\tau}=\left\\{\frac{\bar{Y}_1(1-\bar{Y}_1)}{n_1-1}+\frac{\bar{Y}_0(1-\bar{Y}_0)}{n_0-1}\right\\}^{1/2}.
$$
The Titu's lemma states that $a^2/x+b^2/y\geq (a+b)^2/(x+y)$ when $x,y>0$, followed from which we have
$$
S\_{\tau}^2\geq (n-2)^{-1}\left(\sqrt{\bar{Y}_1(1-\bar{Y}_1)}+\sqrt{\bar{Y}_0(1-\bar{Y}_0)}\right)^2
$$
and thus
$$
n\geq S\_{\tau}^{-2}\left\\{\sqrt{\bar{Y}_0(1-\bar{Y}_0)}+\sqrt{(\bar{Y}_0+\tau)(1-\bar{Y}_0-\tau)}\right\\}^2+2.
$$
Therefore, to detect an ATE of $\tau$ with power at least $\gamma$, roughly we need a sample size of
$$
n \geq \hat{\sigma}\_{\tau}^{-2}(\tau,\gamma)\left\\{\sqrt{\mu_0(1-\mu_0)}+\sqrt{(\mu_0+\tau)(1-\mu_0-\tau)}\right\\}^2+2
$$
where $\hat{\sigma}\_{\tau}(\tau,\gamma)$ is the solution to $\beta(\tau,\sigma\_{\tau})=\gamma$.

When $n_1=n_0$, we have 
$$
n_1=n_0=S\_{\tau}^{-2}\left\\{\bar{Y}_1(1-\bar{Y}_1)+\bar{Y}_0(1-\bar{Y}_0)\right\\} + 1.
$$
In this case, the minimum sample size for each group to detect an ATE of $\tau$ with power at least $\gamma$ is
$$
n_1=n_0= \hat{\sigma}\_{\tau}^{-2}(\tau,\gamma)\left\\{\mu_0(1-\mu_0)+(\mu_0+\tau)(1-\mu_0-\tau)\right\\} + 1.
$$

### Two-Sided P-Value

With the observed data $\mathcal{O}$, a two-sided p-value can be computed by $P(|T|>|T(\mathcal{O})|\mid H_0)=2-2\Phi(|T(\mathcal{O})|)$.

## One-Sided Z-Test

Consider testing the null hypothesis $H_0: \tau>0$ versus the alternative hypothesis $H_1: \tau\leq 0$. The null space of the parameter is $\Theta_0=(0,+\infty)$ and the alternative space is $\Theta_1=(-\infty,0]$.

### T-Statistic

Similarly, we consider the t-statistic $T = T(\mathcal{G})=\hat{\tau}/S\_{\tau}$. Then, for $\tau\in\Theta_0$, we have $T-S\_{\tau}^{-1}\tau\overset{d}{\rightarrow}N(0,1)$.

### The Testing Procedure

Intuitively, $T$ should be significantly greater than $0$ when $H_0$ is correct. The critical region therefore has the form $R(c)=\\{\mathcal{O}: T(\mathcal{O})\leq c \\}$ for some constant $c$.

To control Type I Error at $\alpha$, we can set $c=-z_{\alpha}$. It follows that
$$
\begin{aligned}
P(\mathcal{G}\in R(c)\mid H_0)& \leq \sup_{\tau>0}P_{\tau}(T \leq -z_{\alpha})\\\\
& \leq\sup_{\tau>0}P_{\tau}(T-S\_{\tau}^{-1}\tau +(S\_{\tau}^{-1}-\sigma_{\tau}^{-1})\tau \leq -z_{\alpha}-\sigma_{\tau}^{-1}\tau)\\\\
& \leq\sup_{\tau>0}\Phi(-z_{\alpha}-\sigma_{\tau}^{-1}\tau)+o(1) \leq \alpha+o(1).
\end{aligned}
$$

### Power and Sample Size Determination

The power of the testing procedure is given by
$$
\begin{aligned}
\beta(\tau,\sigma_{\tau})=&P_{\tau}(T \leq -z_{\alpha})\\\\
=& P_{\tau}(T-S\_{\tau}^{-1}\tau +(S\_{\tau}^{-1}-\sigma_{\tau}^{-1})\tau \leq -z_{\alpha}-\sigma_{\tau}^{-1}\tau)\\\\
\approx & \Phi(-z_{\alpha}-\sigma_{\tau}^{-1}\tau)
\end{aligned}
$$
which is trivially a monotonically decreasing function of $\sigma_{\tau}$ when $H_1$ is true. For one-sided z-test, the solution to the equation $\beta(\tau,\sigma_{\tau})=\gamma$ has a closed form:
$$
\sigma_{\tau}(\tau,\gamma) = -\tau/(z_{\alpha}+z_{1-\gamma}) 
$$

Similarly, we consider the settings of binary outcomes. To detect an ATE of $\tau$ with power at least $\gamma$ when $H_1$ is true, roughly we need a sample size of
$$
n \geq (z_{\alpha}+z_{1-\gamma})^2\left\\{\sqrt{\mu_0(1-\mu_0)}+\sqrt{(\mu_0+\tau)(1-\mu_0-\tau)}\right\\}^2/\tau^2+2.
$$
When $n_1=n_0$, the minimum number of samples we need is
$$
n_1 = (z_{\alpha}+z_{1-\gamma})^2\left\\{\mu_0(1-\mu_0)+(\mu_0+\tau)(1-\mu_0-\tau)\right\\}/\tau^2 + 1
$$

### One-Sided P-Value

In one-sided z-test, we can set $W=-T$ which yields the p-value
$$
\sup_{\pi > 0}P_{\pi}(-T\geq -T(\mathcal{O}))=\sup_{\pi > 0}\Phi(T(\mathcal{O})-\tau/\sigma_{\tau})=\Phi(T(\mathcal{O})).
$$

## Z-Test under Covariate-Adaptive Randomization

Under covariate-adaptive randomization, suppose the sample units are divided into $K$ strata indexed by $\\{1,\cdots,K\\}$. We use $Z_i$ to denote a $K$-dimensional indicator variable of strata and assume $Z_i=e_k$ if the unit $i$ is assigned to stratum $k$, where $e_k$ is a $K$-dimentional one-hot vector whose $k$-th component is $1$ and others are $0$. Still, we impose assumptions **A1**--**A3**, except that we modifiy **A2** to that 

* **A2** $\\{(Y_i(1),Y_i(0)),i=1,\cdots,n\\}\perp \\!\\!\\!\perp \\{A_1,\cdots,A_n\\} \mid \\{Z_1,\cdots,Z_n\\}$,

i.e., independence holds for every stratum. Now, we introduce a testing procedure generated by confidence intervals. In other words, we consider interval estimation of $\tau$. 

Let $N_k=\sum_{i=1}^n1\\{Z_i=e_k\\}$ be the sample size of the stratum $k$, $N_{k1}=\sum_{i=1}^n\\{Z_i=e_k\\}A_i$, and $N_{k0}=N_k-N_{k1}$. Let $\bar{Y}\_{ka}=N_{ak}^{-1}\sum_{Z_i=e_k,A_i=a}Y_i$ and $\hat{\tau}\_k=\bar{Y}\_{1k} - \bar{Y}\_{0k}$ be the difference in means estimator of the ATE in the stratum $k$. Then, we can estimate $\tau$ by $\hat{\tau}=\sum_{k=1}^K(N_k/n)\hat{\tau}\_k$, which is unbiased and consistent. The variance of $\hat{\tau}$ is given by
$$
\begin{aligned}
\sigma_{\tau,s}^2=\operatorname{var}(\hat{\tau}) =& E\left\\{\sum_{k=1}^K(N_k^2/n^2)\operatorname{var}(\hat{\tau}\_k\mid Z_1,\cdots,Z_n)\right\\}\\\\
=& E\left\\{\sum_{k=1}^K\frac{N_k^2}{n^2}\left(\frac{\sigma_{1k}^2}{N_{1k}}+\frac{\sigma_{0k}^2}{N_{0k}}\right)\right\\}\\\\
=& \frac{1}{n}\left\\{\frac{E[\operatorname{var}(Y_i(1)\mid Z_i)]}{\pi} + \frac{E[\operatorname{var}(Y_i(0)\mid Z_i)]}{1-\pi}\right\\},
\end{aligned}
$$
where $\sigma_{ak}^2=\operatorname{var}(Y_i(a)\mid Z_i=e_k)$. Note that
$$
\sigma_{\tau}^2 - \sigma_{\tau,s}^2 =\frac{\operatorname{var}[E(Y_i(1)\mid Z_i)]}{n\pi}+\frac{\operatorname{var}[E(Y_i(0)\mid Z_i)]}{n(1-\pi)}\geq 0.
$$
Therefore, owing to stratification, the variance of the estimator for the ATE is reduced and the reduction relies on the between-stratum variance of the outcome variable. We can estimate $\sigma_{\tau,s}^2$ by
$$
S_{\tau,s}^2 = \sum_{k=1}^K\frac{N_k^2}{n^2}\left(\frac{S_{1k}^2}{N_{1k}}+\frac{S_{0k}^2}{N_{0k}}\right)
$$
where $S_{ak}^2=(N_{ak}-1)^{-1}\sum_{A_i=a,Z_i=e_k}(Y_i-\bar{Y}\_{ak})^2$. Under certain regularity conditions, it can be shown that $(\hat{\tau}-\tau)/S_{\tau}$ converges to the standard normal distribution. A two-sided $(1-\alpha)$ confidence interval for $\tau$ is then given as $C_{\alpha}=[\hat{\tau}-S_{\tau,s}z_{\alpha/2},\hat{\tau}+S_{\tau,s}z_{\alpha/2}]$ and we can verify that
$$
P(\tau\in C_{\alpha})=P(|(\hat{\tau}-\tau)/S_{\tau,s}|\leq z_{\alpha/2})=1-\alpha+o(1).
$$
To test the one-sided hypothese $H_0:\tau>0$ versus $H_1:\tau\leq 0$, the $(1-\alpha)$ confidence interval can be given by $C_{\alpha}=(-\infty,\hat{\tau}+S_{\tau,s}z_{\alpha}]$ and a testing procedure can be set as: We reject the null if $0\notin C_{\alpha}(\mathcal{O})$ and accept the null otherwise. The power of the resuling testing procedure under $H_1$ is given by
$$
\begin{aligned}
\beta(\tau,\sigma_{\tau,s}) =& P_{\tau}(\hat{\tau}+S_{\tau,s}z_{\alpha}< 0)\\\\
 =& P_{\tau}((\hat{\tau} - \tau)/S_{\tau,s} < -z_{\alpha} - \tau/S_{\tau,s})\\\\
 =& \Phi(-z_{\alpha}-\sigma_{\tau,s}^{-1}\tau) + o(1),
\end{aligned}
$$
which is a monotonically decreasing function of $\sigma_{\tau,s}$. Because $\sigma_{\tau,s}\leq \sigma_{\tau}$, we have $\beta(\tau,\sigma_{\tau,s})\geq \beta(\tau,\sigma_{\tau})$ under $H_1$.