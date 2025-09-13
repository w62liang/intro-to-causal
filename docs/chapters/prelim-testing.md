---
include: math
---

# Prelimilaries on Hypothesis Testing

> **This section is based on Casella and Berger (2002)**.

Statistical hypothesis testing is a scientific framework for evaluating hypotheses based on observed data, typically made of the following components:

* A statistical model for the data
* Hypothese
* A testing procedure
* Evaluation of the testing procedure
* Reporting testing results

Suppose we have a coin and want to know whether it is fair. In the following, we use this example to illustrate the above concepts.

## Sample Data

>Sample data are a collection of observations with certain measurement variables. 

To assess if the coin is fair, we can flip the coin $n$ times and for each flip we record a measurement as 
$$
X_i = 
\left\\{ 
\begin{array}{cc}
1 & \text{ head} \\\\ 
0 & \text{ tail}
\end{array}
\right. .
$$
Then, the sample data are $\mathbf{X}=\\{X_1,\cdots,X_n\\}$. We use $\mathbf{x}=(x_1,\cdots,x_n)$ to denote a real-valued vector, or say a vector of realizations of $\mathbf{X}$.

## Statistical Models

>A statistical model is a probabilistic model for the data generating process. It is often assumed that the sample data are independently sampled from a hypothetical infinite population. Then, the statistical model is a probabilistic model for the population.

We assume $X_1,\cdots,X_n$ are independent and identically distributed from $\text{Bernoulli}(\pi)$.

## Hypothese

>A hypothesis is a statement about the parameters of the statistical model. The primary interest of hypothesis testing is to figure out whether there is strong statistical evidence against a null hypothesis. An alternative hypothesis is complementary to a null hypothesis and is the one practitioners hope to prove.

We are interested in testing the null $H_0:\pi=0.5$ against the alternative $H_1:\pi\neq 0.5$.

## A Testing Procedure

>A test statistic serves as the foundation for constructing a testing procedure. A critical region, also called rejection region, is a set of values for the sample data constructed upon the test statistic, such that we reject the null when the observed data fall in the critical region, and accept the null when the observed data fall in the acceptance region, complementary to the critical region. A testing procedure is identified by the critical region. 

We consider the test statistic $T=T(\mathbf{X})=n^{-1}\sum_{i=1}^nX_i$. Then, we have $nT\overset{d}{=} B(n,\pi)$ and $nT\mid H_0\overset{d}{=} B(n,0.5)$. If $H_0$ is correct, then $T$ should be close to $0.5$ because the probability mass of $B(n,0.5)$ is concentrated around $0.5n$. Ituitively, therefore, we should reject the null hypothesis when $T$ is far away from $0.5$. Motivated by this reasoning, we define the critical region as
$$
R(c)=\\{\mathbf{x}:|T(\mathbf{x})-0.5|>c\\}
$$ 
for some constant $c\in [0,0.5)$.

## Evaluation of a Testing Procedure

A testing procedure might make two types of errors---Type I Error and Type II Error, which are summarized in the following table:

<table>
  <tr>
    <th></th>
    <th colspan="2" style="text-align:center;">Decision</th>
  </tr>
  <tr>
    <th>Truth</th>
    <td>Accept <i>H</i><sub>0</sub></td>
    <td>Reject <i>H</i><sub>0</sub></td>
  </tr>
  <tr>
    <td><i>H</i><sub>0</sub></td>
    <td>Correct</td>
    <td>Type I Error</td>
  </tr>
  <tr>
    <td><i>H</i><sub>1</sub></td>
    <td>Type II Error</td>
    <td>Correct</td>
  </tr>
</table>

In the coin example, the Type I Error is given by
$$
\begin{aligned}
\alpha(c)&=P(\mathbf{X}\in R(c)\mid H_0)\\\\
 &= P(|T(\mathbf{X})-0.5|>c\mid H_0)\\\\
 &= 1 - P(0.5-c\leq T\leq 0.5+c\mid H_0)\\\\
 &= 1-\sum_{t=\lceil 0.5n-cn \rceil}^{\lfloor 0.5n+cn \rfloor} \binom{n}{t}0.5^t(1-0.5)^{n-t}.
\end{aligned}
$$
The $\alpha(c)$ is complex and not easy to compute. As $n$ goes to infinity, the de Moivre–Laplace theorem states that the binomial distribution $B(n,0.5)$ can be approximated by $N(n/2,n/4)$. In other words, $2\sqrt{n}(T-1/2)\mid H_0$ is asymptotically standard normal which is also a result of the Centre Limit Theorem (CLT). Hence, we have
$$
\alpha(c) = P(2\sqrt{n}|T(\mathbf{X})-0.5|>2\sqrt{n}c\mid H_0) \approx 2 - 2\Phi(2\sqrt{n}c),
$$
which is a monotonically decreasing function of $c$. When $c=0$, we have $R(0)=\\{T(\mathbf{x})\neq 0.5\\}$ and $\alpha(c)\approx 1$. Note that $P(T\neq 0.5)\rightarrow 1$ as $n\rightarrow \infty$. When $c=0.5$, we have $R(0.5)=\empty$ and $\alpha(c)=0$.

The power function is the probability of rejecting $H_0$ under the postulated statistical model, given as
$$
\begin{aligned}
\beta(\pi,c)= & P_{\pi}(\mathbf{X}\in R(c))\\\\
 = & P_{\pi}(T-0.5|>c\mid)\\\\
 = & P_{\pi}\left(\frac{\sqrt{n}(T-\pi)}{\sqrt{\pi(1-\pi)}}<\frac{\sqrt{n}(0.5-c-\pi)}{\sqrt{\pi(1-\pi)}}\right) \\\\
 &+ P_{\pi}\left(\frac{\sqrt{n}(T-\pi)}{\sqrt{\pi(1-\pi)}}>\frac{\sqrt{n}(0.5+c-\pi)}{\sqrt{\pi(1-\pi)}}\right) \\\\
 \approx & \Phi\left(\frac{\sqrt{n}(0.5-\pi-c)}{\sqrt{\pi(1-\pi)}}\right) + 1 - \Phi\left(\frac{\sqrt{n}(0.5-\pi+c)}{\sqrt{\pi(1-\pi)}}\right),
\end{aligned}
$$
which is also a monotonically decreasing function of $c$. Fixing $\pi>0$, the Type II Error is $1-\beta(\pi,c)$ as a monotonically increasing function of $c$. We see that there is a balance between the two types of errors regarding the values of $c$; it is difficult to control the two types of errors simultaneously at a low level by adjusting the testing procedure. The primary goal is to control the Type I Error. The Type II Error can be 
reduced by increasing the sample size.

## Reporting Testing Results

In addition to the decision of rejecting or accepting a null hypothesis, the practitioners always want to know more about the strength of the evidence of making such a decision. We introduce two commonly reported strength measures of the statistical evidence.

### P-Value

>A p-value, denoted by $p(\mathbf{X})$, is a test statistic that lies in $[0,1]$ such that smaller values of $p(\mathbf{X})$ imply stronger evidence against the null hypothesis. A valid p-value should satisfy $P(p(\mathbf{X})\leq \alpha\mid H_0)\leq\alpha$. 

Assume $W=W(\mathbf{X})$ is a test statistic so that large value of $W$ yields stronger evidence against $H_0$. We can construct a p-value by
$$
p(\mathbf{x})=P(W(\mathbf{X})\geq W(\mathbf{x})\mid H_0).
$$
Because large value of $W$ is unlikely to appear under $H_0$, the $p(\mathbf{x})$ reveals the tail probability of observed data under $H_0$. Let $F_{W\mid H_0}(w)$ denote the cumulative distribution function (cdf) of $W\mid H_0$. Then, $F_{W\mid H_0}(W)\mid H_0$ is uniformly distributed on $[0,1]$. One can verify that
$$
\begin{aligned}
P(p(\mathbf{X})\leq \alpha\mid H_0) &= P(1-F_{W\mid H_0}(W)\leq \alpha\mid H_0)\\\\
&= P(F_{W\mid H_0}(W)\geq 1-\alpha\mid H_0)\leq \alpha,
\end{aligned}
$$
implying $p(\mathbf{X})$ is valid. In the coin example, we can set $W=|T-1/2|$ to obtain a valid p-value.

Provided a valid p-value, we can construct the critical region as $R_{\alpha}=\\{\mathbf{x}: p(\mathbf{x})\leq \alpha\\}$. By definition of p-value, it is clear that $P(\mathbf{X}\in R_{\alpha}\mid H_0)=P(p(\mathbf{X})\leq \alpha\mid H_0)\leq\alpha$. A testing procedure generated from a valid p-value always controls the Type I Error.

### Confidence Interval

>Given a significance level $\alpha$, a confidence interval (CI), denoted by $C_{\alpha}=[L(\mathbf{X}),U(\mathbf{X})]$, is a random interval for the population parameter of interest such that $P(\pi\in C_{\alpha})\geq 1-\alpha$, implying that the CI has at least $100(1-\alpha)\%$ chance to cover the true parameter.

The key to find a CI is the so-called pivitol quantity, which is a random variable and also a function of the population parameter such that its distribution is known and independent of all parameters in the statistical model. 

Let $Q=Q(\mathbf{X},\pi)=\sqrt{n}(T-\pi)/\sqrt{\pi(1-\pi)}$ be the scalar pivitol quantity which is asymptotically standard normal. Then, for any $\alpha_1,\alpha_2\in [0,1]$ such that $\alpha_1+\alpha_2=\alpha$, a $(1-\alpha)$ confidence interval for $\pi$ can be constructed as $C_{\alpha}=C_{\alpha}(\mathbf{X}) = \\{ \pi: z_{1-\alpha_1}< Q(\mathbf{X},\pi) \leq z_{\alpha_2} \\}$. A critical region can be set up based on $C_{\alpha}$ as 
$$
R_{\alpha} = \\{\mathbf{x}: 0.5\notin C_{\alpha}(\mathbf{x})\\} = \\{\mathbf{x}: Q(\mathbf{x},0.5)\leq z_{1-\alpha_1}\text{ or }Q(\mathbf{x},0.5) > z_{\alpha_2}\\}.
$$ 
It follows that
$$
P(\mathbf{X}\in R_{\alpha}\mid H_0) = 1-P(Q(\mathbf{X},0.5)\in (z_{1-\alpha_1},z_{\alpha_2}) \mid H_0)=\alpha+o(1).
$$

## References

<cite>Casella, G., & Berger, R. (2002). Statistical inference. Chapman and Hall/CRC.</cite>