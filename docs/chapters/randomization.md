---
include: math
---

# Randomization Methods

In randomized experiments, treatment assignment mechanisms are determined by practitioners at the design stage. Randomization refers to a class of probablisitic treatment assignment schemes which provide a sound statistical basis for causal inference. Regarding theoretical aspects, the randomization methods can be generally classified into two categories: covariate-free methods that do not rely on baseline covariates and covariate-adaptive methods. The first results in independence between the treatment indicators and other factors including the potential outcomes, baseline covariates, and unmeasured confounding variables. The second results in independence across strata. 

Let $A_1,\cdots,A_n$ denote treatment indicators. Under covariate-free randomization, we have
$$
\\{A_1,\cdots,A_n\\} \perp\\!\\!\\!\perp\\{(Y_i(1),Y_i(0)),i=1,\cdots,n\\}.
$$
Therefore, $\hat{\tau} = \bar{Y}\_1 - \bar{Y}\_0$ is an unbiased estimator of $\tau=E[Y(1)-Y(0)]$.
For covariate-adaptive methods such as stratified permuted block, stratified biased coin, and Pocock-Simon’s minimization, certain baseline covariates are used to create strata and covariate-free randomization is carried out in each strata. Then, the independence assumption holds across strata, i.e.,
$$
\\{A_1,\cdots,A_n\\} \perp\\!\\!\\!\perp\\{(Y_i(1),Y_i(0)),i=1,\cdots,n\\}\mid \\{Z_1,\cdots,Z_n\\}
$$
where $Z_1,\cdots,Z_n$ are categorical variables indicating stratum membership. Suppose the observations are split into $K$ strata indexed by $\\{1,\cdots,K\\}$ and the sizes of the strata are $N_1,\cdots,N_K$. For covariate-adaptive randomization, we can estimate $\tau$ by $\hat{\tau}=\sum_{k=1}^K(N_k/n)\hat{\tau}_k$ where $\hat{\tau}_k$ are the stratum-specific ATEs.

## Completely Randomized Experiment
In completely randomized experiment (CRE), the sizes of treated and control, $n_1$ and $n_0$, are pre-determined and each possible treatment assignment sequence occurs with equal probability. More formally, for $(a_1,\cdots,a_n)\in\\{0,1\\}^n$ we have 
$$
P(A_1=a_1,\cdots,A_n=a_n)=\left\\{ 
\begin{array}{cc}
1/\binom{n}{n_1} & \text{if} \sum_{i=1}^na_i=n_1 \\\\ 
0 & \text{o.w.}
\end{array}
\right. .
$$
Under this treatment assignment scheme, one can verify that
$$
\begin{aligned}
P(A_i=1) =& \sum_{\sum_{j\neq i}a_j=n_1-1}P(A_1=a_1,\cdots,A_i=1,\cdots,A_n=a_n)
=\frac{\binom{n-1}{n_1-1}}{\binom{n}{n_1}}
=\frac{n_1}{n}
\end{aligned}
$$
for $i=1,\cdots,n$. The CRE can be realized by randomly shuffling the sample sequence and picking up the first $n_1$ units to treat. Fisher–Yates shuffle is a classical algorithm for shuffling a finite sequence, where for $i$ from $n$ to $1$ we randomly select $j$ from $\\{1,\cdots,i\\}$ and exchange the $i$-the and $j$-th unit. For a CRE with $K>2$ treatment arms, define $e_k$ as a $k$-dimensional vector with the $k$-th component being $1$ and the others being $0$. We let $A_i=e_k$ if the $i$-th unit is assinged to treatment group $k$. Let $n_1,\cdots,n_K$ be pre-specified sample sizes for the $K$ treatment groups. Then, for $(a_1,\cdots,a_n)\in\\{e_1,\cdots,e_K\\}^n$, the treatment assignment mechanism is
$$
P(A_1=a_1,\cdots,A_n=a_n)=\left\\{ 
\begin{array}{cc}
\frac{\prod_{k=1}^Kn_k!}{n!} & \text{if} \sum_{i=1}^n a_i=(n_1,\cdots,n_K) \\\\ 
0 & \text{o.w.}
\end{array}
\right. .
$$
Similarly, to realize CREs with multiple treatment arms, we can shuffle the sample units and assign the $k$-th block of $n_k$ sample units to the treatment group $k$. One can prove that the probability of any unit being selected to treatment group $k$ is $n_k/n$, which is left as an exercise.

The CRE achieves balance in group size and is easy to implement. However, it has a few drawbacks:

* The CRE is not time-adaptive. The treatment assignment scheme relies on the sample size and thus cannot be applied to experiments with consecutive enrollment process. 
* The CRE has a chance to yield covariate imbalance especially when the sample size is small, leading to the loss of estimator's efficiency.

To address the first drawback, we introduce permuted block randomization method and biased coin randomization.

## Permuted Block Randomization

In permuted block randomization (PBR), individuals are randomly filled in consecutive blocks of even size such as $4$, $6$ and $8$, and the CRE is applied to each block. For example, if we hope to balance the treatment group size, we can create consecutive $4$-unit blocks and set a random treatment assignment sequence (TTCC,TCTC,CCTT, CTCT,TCCT,CTTC) for each block. In a three-arm experiment, we can create blocks of $6$ units with treatment ratio such as $2:2:2$ or $1:2:3$.

## Biased Coin Randomization
Biased coin randomization (BCR) is a variation of multinomial trials with dynamic treatment probabilities. The multinomial trials randomization is also called simple randomization, which randomly assigns treatments in a K-arm experiment with fixed probabilities $p_1,\cdots,p_K$. Despite time-adaptive, it may cause serious group imbalance during the study. BCR attempts to overcome imbalance in group size by adjusting the treatment probabilities dynamically based on the proportion of treated units. For example, in a two-arm experiment, the treatment probability is set by $p=0.5$ at the beginning. Suppose $10$ participants are enrolled and $8$ of them are randomly assigned to the treated group. To balance group size, we can adjust the treatment probability to $p=0.2$ for the next $10$ sample units.

The PBR and BCR are both time-adaptive strategies which fit in long-term randmized experiments with consecutive recruiting phase. During the entire study, the sizes of the treatment and control groups are always well balanced. However, they do not resolve the second drawback of a CRE; There might be covariate-imbalance between treatment groups when the sample size is small. This issue can be addressed at either the design stage or the analysis stage. We introduce an experimental design called stratified randomization which mitigates covariate imbalance at the design stage.

## Stratified Randomization

Stratified randomized experiments (SREs) are covariate-adaptive. Baseline covariates are collected prior to treatment assignment and selected risk factors considered associated with the potential outcomes are used to create strata. Covariate-free randomization is then carried out in each stratum. For example, if gender is considered prognostic for the outcome, the practitioners can perform CRE, PBR, or BCR separately for males and females, thereby ensuring that the proportion of males is balanced across treatment groups. For a continuous variable such as age, the practitioners can dichotomize or categorize the variable into intervals, e.g., under 18 years, 18 - 60 years, and 60 years and older. When there are two and more variables for stratification, the number of strata can grow rapidly, which may lead to small sample sizes within certain strata. This potential limitation of SREs highlights the importance of limiting the number of stratification variables.

## Matched Paired Design

## Pocock-Simon’s Minimization

