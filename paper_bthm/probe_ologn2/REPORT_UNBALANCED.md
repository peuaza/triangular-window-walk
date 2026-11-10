# The single-epoch lower bound for unbalanced multisets

**Question.** REPORT_EPOCH_BOUND proves, for a uniform arrangement of a multiset $\mathfrak M\subset\{0,1,2\}^M$
with $\kappa_0=\kappa_2$ and $s^2\ge 1/4$,

$$P(\mathrm{cost}\ge cM)\ \ge\ p(c)-o(1),\qquad
p(c):=P\Big(\sup_{\theta\in[1-c,1]}\mathcal V(\theta)<0\Big),\quad
\mathcal V(\theta)=W(\theta)-2W(\theta/2),$$

and §8 of the paper applies it to a **suffix** $a[i..n-1]$, whose counts satisfy
$\#2-\#0=-X_i$ and are therefore *not* balanced. REFEREE_S8 §4(c) records this as an unproved
extension and calls the patch in REPORT_EPOCH_BOUND §2 ("the drift term $\mu\theta$ is killed by the
contrast, since $\mu(\theta-2\theta/2)=0$") dimensionally wrong.

**Answer: the extension is true and is proved below.** The correct statement is

> **Theorem 1.** For every $c\in(0,1)$ and $\eta>0$ there is $M_0=M_0(c,\eta)$ such that for every
> $M\ge M_0$, **every** multiset $\mathfrak M\subset\{0,1,2\}$ of $M$ cells with $s^2(\mathfrak M)\ge 1/4$,
> and the uniform random arrangement $a$ of $\mathfrak M$,
> $$P(\mathrm{cost}(a)\ge cM)\ \ge\ p(c)-\eta .$$
> $M_0$ does not depend on the composition $(\kappa_0,\kappa_1,\kappa_2)$.

So the second hypothesis of the $\Omega(n^2)$ bound is removable, and the paper's
$c_0=c\,p(c)$ is unchanged. The referee's diagnosis is right in every particular: the drift *is* of
order $\sqrt M$ — $(\mu-1)\sqrt M\,\theta/s$ in the notation below, the referee's $(\mu/s)\sqrt
M\,\theta$ with $\mu$ the mean of $b=a-1$ — and it diverges; but it never enters, because the test is
**exactly** invariant under translation of the cells (Lemma 1), so one centres
the cells *before* normalising and then the limit is a genuine bridge for **every** composition
(Lemma 2). The one genuinely delicate point — uniformity of the CLT over the composition — is
Lemma 3, proved from scratch by Stein's method of exchangeable pairs with an explicit error
$2.89\,(K/s)/\sqrt m+O(1/M)$ that depends only on $K/s\le 4$. No functional CLT, no tightness, no
entropy bound and no combinatorial CLT are cited.

| statement | status |
|---|---|
| cost is invariant under $a\mapsto \alpha a+\beta\mathbf 1$ ($\alpha>0$), exactly, pathwise | **proved** (Lemma 1) |
| hence $\mathcal L(\mathrm{cost}\mid\mathfrak M)=\mathcal L(\mathrm{cost}\mid\mathfrak M-\mu\mathbf 1)$, exactly | **proved** (Cor. 1) |
| exact covariance $\mathrm{Cov}(Y_i,Y_j)=s^2i(M-j)/(M-1)$ for **any** multiset | **proved** (Lemma 2(a)) |
| exact contrast covariance $(2(\theta\wedge\theta')-(\theta\vee\theta'))_+M/(M-1)$ for **any** multiset | **proved** (Lemma 2(b)); exact check, 82 compositions |
| quantitative CLT for linear permutation statistics, error uniform in the composition | **proved** (Lemma 3) |
| discrete scan $\to$ Brownian functional, uniform constants with $K=2$ | **proved** (Lemma 4) |
| Theorem 1 (liminf $\ge p(c)$, uniform over compositions) | **proved** (§6, Prop. 6.1) |
| the form §8 needs: $E[\mathrm{cost}(i,n-1)]\ge c(n-i)(p(c)-o(1))$ | **proved** (Cor. 2) |
| $P(s^2(\text{suffix})<1/4)=e^{-\Theta(M)}$ for the paper's ensembles | **proved** (Lemma 5) |
| Stein-equation bounds $\sup|f'|\le 1.2533,\ \sup|f''|\le 2$ | *cited* (Stein 1986) |
| Hoeffding's hypergeometric tail bound | *cited* (Hoeffding 1963) |

Everything is uniform in $(\kappa_0,\kappa_1,\kappa_2)$; the only numerical inputs are
$|a_i-\mu|\le 2$ and $s^2\ge 1/4$.

---

## 1. Setting and notation

$M$ cells; the multiset $\mathfrak M\subset\{0,1,2\}$ has counts $(\kappa_0,\kappa_1,\kappa_2)$,
$\kappa_0+\kappa_1+\kappa_2=M$, mean and variance

$$\mu=\frac{\kappa_1+2\kappa_2}{M}\in[0,2],\qquad
s^2=\frac1M\sum_{v\in\mathfrak M}(v-\mu)^2=\frac{\kappa_1+4\kappa_2}{M}-\mu^2 .$$

$a=(a_1,\dots,a_M)$ is a uniform random arrangement of $\mathfrak M$ (uniform over the
$M!/\prod\kappa_j!$ distinct words). Put

$$c_i:=a_i-\mu,\qquad Y_m:=\sum_{i\le m}c_i,\qquad A_m:=\frac1m\sum_{i\le m}a_i=\mu+\frac{Y_m}{m},$$

so that $\sum_{i=1}^M c_i=0$ **identically** and $E[Y_m]=0$ for every $m$ and every composition.
With $k(m)=\lfloor m/2\rfloor$,

$$b_2(m)=1\iff A_{k(m)}<A_m\iff
\begin{cases} Y_m-2Y_{k}>0, & m=2k\ \text{even},\\[2pt] Y_m-2Y_k>Y_k/k, & m=2k+1\ \text{odd},\end{cases}
\tag{1.1}$$

and $C=\{m\in[2,M]: b_2(m)=0\}$, $m^*=\max C$, $\mathrm{cost}=M-m^*+1$ ($=M-1$ if $C=\varnothing$).
Finally

$$\mathcal D_M(\theta):=\frac{Y_{\lfloor\theta M\rfloor}-2Y_{\lfloor\theta M/2\rfloor}}{s\sqrt M},
\qquad \mathcal G_M(m):=\frac{Y_m-2Y_{\lfloor m/2\rfloor}}{s\sqrt M}=\mathcal D_M(m/M),
\tag{1.2}$$

the second identity being exact for every $m$ — this is the only "floor" fact used below.
$\mathcal V(\theta)=W(\theta)-2W(\theta/2)$ with $W$ a standard Brownian motion, and
$p(c)=P(\sup_{[1-c,1]}\mathcal V<0)$.

Note $|c_i|\le K:=2$ always, since $a_i,\mu\in[0,2]$. $s\ge 1/2$ by
hypothesis, hence

$$K/s\ \le\ 4\qquad\text{for every admissible multiset.}\tag{1.3}$$

(For a two-point multiset the bound is much better: e.g. $\{1,2\}$ with equal counts has $K=1/2$,
$s=1/2$, $K/s=1$.)

---

## 2. Lemma 1 (exact invariance) — and what the old remark got wrong

**Lemma 1.** For $\alpha>0$, $\beta\in\mathbb R$ and any word $a$,
$\mathrm{cost}(\alpha a+\beta\mathbf 1)=\mathrm{cost}(a)$.

*Proof.* $A_m(\alpha a+\beta\mathbf1)=\alpha A_m(a)+\beta$ for every $m$, so every comparison
$A_{k(m)}\ge A_m$ is unchanged. $\ \blacksquare$

**Corollary 1.** For a uniform arrangement of $\mathfrak M$,
$\mathcal L(\mathrm{cost}\mid\mathfrak M)=\mathcal L(\mathrm{cost}\mid\mathfrak M-\mu\mathbf 1)$,
exactly, at every $M$. Hence in proving Theorem 1 we may assume $\mu=0$, i.e. $\sum_ic_i=0$.

This is an identity of functions, not an asymptotic statement: nothing about the contrast, the CLT
or the normalisation enters. It is the whole content of "the drift does not matter", and it is
strictly stronger than the referee's step (i).

### 2.1 Why the printed justification was wrong, and what the right statement is

Write $X_m=\sum_{i\le m}(a_i-1)=Y_m+(\mu-1)m$. In the *raw* units the contrast is

$$X_m-2X_{\lfloor m/2\rfloor}=Y_m-2Y_{\lfloor m/2\rfloor}+(\mu-1)\bigl(m-2\lfloor m/2\rfloor\bigr)
=Y_m-2Y_{\lfloor m/2\rfloor}\quad(m\text{ even}),$$

because $m-2(m/2)=0$ **exactly**. So the drift cancels pathwise, at finite $M$, in the contrast —
that part of the report's remark is correct. What is *not* correct is the location of the
cancellation. In the normalised prefix process

$$\frac{X_{\lfloor\theta M\rfloor}}{s\sqrt M}=\frac{\mu-1}{s}\sqrt M\,\theta+\frac{Y_{\lfloor\theta
M\rfloor}}{s\sqrt M},$$

the drift is $(\mu-1)\sqrt M\theta/s$ and **diverges**; there is no "$\mu\theta$ term in $\mathfrak B$"
that a computation $\mu(\theta-2\theta/2)=0$ could remove, and the displayed cancellation is indeed
dimensionally meaningless as written. The correct statement is the pair

* *(un-normalised, exact)* $X_m-2X_{\lfloor m/2\rfloor}=Y_m-2Y_{\lfloor m/2\rfloor}$ for even $m$,
  and $X_m-2X_k-\frac{X_k}{k}=Y_m-2Y_k-\frac{Y_k}{k}$ for odd $m=2k+1$;
* *(normalised)* therefore $\mathcal D_M(\theta)=(X_{\lfloor\theta M\rfloor}-2X_{\lfloor\theta
  M/2\rfloor})/(s\sqrt M)$ **is** the contrast of the centred process, and the centred process is a
  bridge, not a drifted Brownian motion, because $\sum_ic_i=0$ exactly.

So the drift is removed by centring the cells, which Lemma 1 licenses exactly — never by
cancelling a term inside a normalised process that contains it at order $\sqrt M$. In the centred
coordinates there is no drift at all, and the balanced case is not special: the "bridge" in
REPORT_EPOCH_BOUND Lemma 1 came from $\sum a=2h$, and after centring $\sum c=0$ holds for every
multiset, so the limit is a bridge for every composition.

---

## 3. Lemma 2 (exact covariance for any multiset)

**Lemma 2.** Let $\mathfrak M$ be **any** multiset of $M$ reals with variance $s^2>0$ (no hypothesis
on the mean, none on $\kappa_0-\kappa_2$), $a$ a uniform arrangement, $Y_m=\sum_{i\le m}(a_i-\mu)$.
Then

**(a)** for $0\le i\le j\le M$,
$$\mathrm{Cov}(Y_i,Y_j)=s^2\,\frac{i(M-j)}{M-1};$$

**(b)** for $\theta,\theta'\in(0,1]$ with $\theta M,\theta M/2,\theta'M,\theta'M/2$ integral,
$$\mathrm{Cov}\bigl(\mathcal D_M(\theta),\mathcal D_M(\theta')\bigr)
=\bigl(2(\theta\wedge\theta')-(\theta\vee\theta')\bigr)_+\cdot\frac{M}{M-1};$$

**(c)** in particular $\mathcal D_M(\theta)\perp\mathcal D_M(\theta')$ **exactly** whenever
$\theta'\ge 2\theta$, and $\mathrm{Var}\,\mathcal D_M(\theta)=\theta M/(M-1)$.

*Proof.* (a) $Y_i$ is a sum over a simple random sample of $i$ cells. Each cell contributes with
probability $i/M$, so
$\mathrm{Var}(Y_i)=\sum_{\text{cells}}c^2\frac{i}{M}\bigl(1-\frac{i}{M}\bigr)
+2\sum_{c<c'}cc'\,\mathrm{Cov}(\mathbf 1_c,\mathbf 1_{c'})$
with $\mathrm{Cov}(\mathbf 1_c,\mathbf 1_{c'})=-\frac{i(M-i)}{M^2(M-1)}$; using $\sum c=0$ this gives
$\mathrm{Var}(Y_i)=s^2i(M-i)/(M-1)$. For $i\le j$, $Y_j-Y_i$ is the sum over the next $j-i$ draws and
$\mathrm{Cov}(Y_i,Y_j-Y_i)=-s^2i(j-i)/(M-1)$ (covariance of disjoint parts of a WOR sample), so
$\mathrm{Cov}(Y_i,Y_j)=s^2[i(M-i)-i(j-i)]/(M-1)=s^2i(M-j)/(M-1)$.

(b) Put $i=\theta M$, $j=\theta'M$ and assume $i\le j$. Then
$$\mathrm{Cov}(\mathcal D(\theta),\mathcal D(\theta'))=\frac1{s^2M}\Bigl[
\mathrm{Cov}(Y_i,Y_j)-2\mathrm{Cov}(Y_i,Y_{j/2})-2\mathrm{Cov}(Y_{i/2},Y_j)+4\mathrm{Cov}(Y_{i/2},Y_{j/2})\Bigr].$$
If $j\ge 2i$ the four terms cancel identically: using (a) with the orders
$i\le j$, $i\le j/2$, $i/2\le j$, $i/2\le j/2$,
$$s^{-2}(M-1)\sum=i(M-j)-2i(M-\tfrac j2)-2\tfrac i2(M-j)+4\tfrac i2(M-\tfrac j2)
=\bigl(i-2i-i+2i\bigr)M+\bigl(-1+1+1-1\bigr)ij=0 .$$
If $i\le j\le 2i$ then $j/2\le i$ and the orders are $i\le j$, $j/2\le i$, $i/2\le j$, $i/2\le j/2$:
$$s^{-2}(M-1)\sum=i(M-j)-j(M-i)-i(M-j)+2i(M-\tfrac j2)=(2i-j)M,$$
i.e. $\mathrm{Cov}=(2\theta-\theta')M/(M-1)$, which is (b) for $\theta\le\theta'$ (and gives
$\theta M/(M-1)$ on the diagonal).

(c) is the case $j\ge2i$ of the display. $\ \blacksquare$

**Remarks.** (i) Nothing in Lemma 2 uses $\kappa_0=\kappa_2$, the value set $\{0,1,2\}$, or
$s^2\ge1/4$ — only $s^2>0$. (ii) Verification is exact and cheap: enumerating all arrangements of all
82 compositions at $M=12$ (and 30 at $M=8$), the largest deviation from the formula in (b) is
$7.6\times10^{-15}$, i.e. round-off (q13 part B). The referee's checks of $(1,3,4)$ and $(0,4,4)$ are
instances. (iii) (b) is what makes the *limit* law composition-free: after dividing by $s^2M$ the
covariance is the Brownian one times $M/(M-1)$, with $s$ scaling out. So the only way the
composition could enter the limit is through the CLT error — which is Lemma 3.

---

## 4. Lemma 3 (quantitative CLT, uniform in the composition) — the delicate step

This is the step the referee flags: "the CLT must be uniform over the possible compositions". The
following gives an explicit error depending only on $K/s$, hence — by (1.3) — bounded by a universal
constant.

**Lemma 3.** Let $c_1,\dots,c_M\in\mathbb R$ with $\sum_ic_i=0$, $K=\max_i|c_i|$,
$s^2=M^{-1}\sum_ic_i^2>0$. Let $u_1,\dots,u_M\in\mathbb R$ and put
$A_r=\sum_iu_i^r$ $(r=1,2,4)$, $U=\max_i|u_i|$. Assume $A_1^2\le A_2/4$. Let $\pi$ be a uniform
random permutation of $\{1,\dots,M\}$ and
$$Z=\sum_{i=1}^M u_i\,c_{\pi(i)},\qquad
\sigma^2=\mathrm{Var}\,Z=s^2\,\frac{MA_2-A_1^2}{M-1}>0 .$$
Then for every $z\in\mathbb R$,
$$\Bigl|P(Z\le z)-\Phi(z/\sigma)\Bigr|\ \le\
\frac{1.2533}{\sqrt2}\,\frac{K}{s}\,\frac{\sqrt{A_4}}{A_2}
\;+\;\frac{2KU}{s\sqrt{A_2}}\;+\;\frac1M .
\tag{4.1}$$

*Proof.* Write $W=Z$.

**Step 1 (mean and variance).** $E[W]=\frac{A_1}{M}\sum_ic_i=0$. For $i\ne j$,
$E[c_{\pi(i)}c_{\pi(j)}]=\frac{1}{M(M-1)}\sum_{k\ne l}c_kc_l=-\frac{s^2}{M-1}$; hence
$\mathrm{Var}(W)=s^2A_2-\frac{s^2}{M-1}(A_1^2-A_2)=s^2\frac{MA_2-A_1^2}{M-1}$.

**Step 2 (exchangeable pair).** Let $I$ be uniform on $\{1,\dots,M\}$ and let $\pi'$ be $\pi$ with
position $I$ re-randomised: $\pi'(J)=\pi(J)$ for $J\ne I$ and $\pi'(I)$ uniform on the $M-1$ values
not used by $\pi$ elsewhere. The transition kernel is symmetric — from either $\pi$ or $\pi'$ the
move has probability $\frac1M\cdot\frac1{M-1}$ — so with $W'=\sum_iu_ic_{\pi'(i)}$ the pair
$(W,W')$ is exchangeable. Since
$E[c_{\pi'(I)}\mid\pi,I=i]=\frac{1}{M-1}\sum_{j\ne\pi(i)}c_j=-\frac{c_{\pi(i)}}{M-1}$,
$$E[W'\mid\pi]=W-\frac1M\sum_i u_i\frac{c_{\pi(i)}}{M-1}=W\Bigl(1-\frac1{M-1}\Bigr),\qquad
\lambda:=\frac1{M-1}.\tag{4.2}$$
Moreover, conditionally on $\pi$ and on $\pi(i)=v$,
$E[(c_{\pi'(i)}-v)^2]=\frac{1}{M-1}\sum_{j\ne\pi(i)}(c_j-v)^2=\frac{M(s^2+v^2)}{M-1}$, so with
$D:=\sum_iu_i^2c_{\pi(i)}^2$,
$$E[(W'-W)^2\mid\pi]=\frac{s^2A_2+D}{M-1},\qquad
\alpha:=E[(W'-W)^2]=\frac{2s^2A_2}{M-1},\qquad
\sigma_\lambda^2:=\frac{\alpha}{2\lambda}=s^2A_2.\tag{4.3}$$

**Step 3 (Stein).** Let $|h|\le1$ and let $f$ solve the Stein equation
$f'(w)-\sigma_\lambda^{-2}wf(w)=h(w)-E[h(\sigma_\lambda G)]$, $G\sim N(0,1)$. Then
$E[h(W)]-E[h(\sigma_\lambda G)]=E[f'(W)-\sigma_\lambda^{-2}Wf(W)]$. By exchangeability
$E[(W'-W)(f(W')+f(W))]=0$, whence $E[(W'-W)(f(W')-f(W))]=-2E[(W'-W)f(W)]$, and by (4.2)
$E[(W'-W)f(W)]=-\lambda E[Wf(W)]$. Therefore
$E[Wf(W)]=\frac1{2\lambda}E[(W'-W)(f(W')-f(W))]$. Taylor's theorem gives
$f(W')-f(W)=(W'-W)f'(W)+\frac12(W'-W)^2f''(W_*)$, so
$$E[h(W)]-E[h(\sigma_\lambda G)]
=E\Bigl[f'(W)\Bigl(1-\frac{E[(W'-W)^2\mid W]}{\alpha}\Bigr)\Bigr]
-\frac1{2\alpha}E\bigl[(W'-W)^3f''(W_*)\bigr],$$
using $2\lambda\sigma_\lambda^2=\alpha$. Hence
$$\bigl|E[h(W)]-E[h(\sigma_\lambda G)]\bigr|\le
\sup|f'|\;E\Bigl|1-\frac{E[(W'-W)^2\mid W]}{\alpha}\Bigr|
+\frac{\sup|f''|}{2\alpha}E|W'-W|^3 .\tag{4.4}$$

**Step 4 (the two terms).** Since $|W'-W|\le 2KU$ and $E[(W'-W)^2]=\alpha$,
$E|W'-W|^3\le 2KU\alpha$, and the Stein bounds ($\sup|f'|\le1.2533$, $\sup|f''|\le2/\sigma_\lambda$)
give
$$\text{second term}\ \le\ \frac{2}{\sigma_\lambda}\cdot\frac{2KU\alpha}{2\alpha}
=\frac{2KU}{s\sqrt{A_2}} .$$
For the first term, $E[(W'-W)^2\mid W]-\alpha=(D-E D)/(M-1)$, so it equals
$E|D-ED|/((M-1)\alpha)=E|D-ED|/(2s^2A_2)$. Now $D=\sum_iv_ix_{\pi(i)}$ with $v_i=u_i^2$,
$x_i=c_i^2$, so by the variance formula used in Step 1,
$\mathrm{Var}(D)=\mathrm{Var}(x)\,(MA_4-A_2^2)/(M-1)$ with
$\mathrm{Var}(x)\le E[x^2]=M^{-1}\sum c_i^4\le K^2s^2$. Hence
$\mathrm{sd}(D)\le Ks\sqrt{2A_4}$ and
$$\text{first term}\ \le\ \frac{Ks\sqrt{2A_4}}{2s^2A_2}=\frac{K}{\sqrt2\,s}\frac{\sqrt{A_4}}{A_2}.$$

**Step 5 (replace $\sigma_\lambda$ by $\sigma$).** From (4.3) and Step 1,
$\bigl|\sigma_\lambda^2/\sigma^2-1\bigr|=|A_1^2-A_2|/(MA_2-A_1^2)\le \frac{A_2/4+A_2}{MA_2-A_2/4}\le\frac2M$
(using the hypothesis $A_1^2\le A_2/4$), and since
$\sup_z|\Phi(z/a)-\Phi(z/b)|\le 0.4\,|a/b-1|$ for $a,b>0$, the substitution costs at most $0.4/M$.
Adding, and using $1.2533/\sqrt2=0.886$, gives (4.1). $\ \blacksquare$

**Remarks.**
1. *(Rate.)* For a single WOR sum ($u_i=\mathbf 1\{i\le m\}$) this is $A_2=A_4=m$, $U=1$, and (4.1)
   reads $0.886(K/s)/\sqrt m+2(K/s)/\sqrt m+O(1/M)$ — the classical $O(1/\sqrt m)$ Berry–Esseen rate
   for sampling without replacement, with an explicit constant.
2. *(Uniformity.)* The right side of (4.1) involves the multiset only through $K/s\le4$; $A_2,A_4,U$
   are functions of $u$ alone. **This is the uniformity the referee asks for**: no composition, and
   no lower bound on $s$ beyond $s\ge1/2$, enters.
3. *(Vectors.)* For fixed $\theta_1<\dots<\theta_\ell$ and $\alpha\in\mathbb R^\ell$, apply (4.1) to
   $u_i=\sum_r\alpha_r\bigl(\mathbf 1\{i\le\lfloor\theta_rM\rfloor\}-2\cdot\mathbf 1\{i\le\lfloor\theta_rM/2\rfloor\}\bigr)$,
   for which $Z/(s\sqrt M)=\sum_r\alpha_r\mathcal D_M(\theta_r)$. Then $A_2=M\,v_\alpha(1+O(1/M))$ with
   $v_\alpha=\int_0^1\bigl(\sum_r\alpha_r\varphi_{\theta_r}\bigr)^2dt>0$,
   $\varphi_\theta=\mathbf 1\{t\le\theta\}-2\cdot\mathbf 1\{t\le\theta/2\}$, since distinct $\theta$'s
   give linearly independent $\varphi$'s; $U\le3\sum|\alpha_r|$; $A_1=O(1)$ so $A_1^2\le A_2/4$ for
   $M$ large. So the error is $O_{\alpha}(K/s)/\sqrt M\to0$ **uniformly in the multiset**. Cramér–Wold
   then gives joint convergence of $(\mathcal D_M(\theta_r))_r$ to the centred Gaussian vector with
   covariance $\Sigma_{\theta\theta'}=(2(\theta\wedge\theta')-(\theta\vee\theta'))_+$, which by
   Lemma 2(b) is the covariance of $(\mathcal V(\theta_r))_r$.
4. *(No tightness needed.)* Everything below uses only these finite-dimensional statements. There is
   no functional CLT, no Serfling inequality, no entropy bound. This is why the uniformity is easy to
   certify: the only place a composition could have entered the asymptotics is (4.1), and it does not.

---

## 5. Lemma 4 (from the Brownian functional to the discrete scan)

**Lemma 4.** Let $a$ be a word of length $M$ with variance $s^2>0$, and let $\mathcal G_M$,
$\mathcal D_M$ be as in (1.2). Then for every $m\in[2,M]$:

**(a)** if $m$ is even, $b_2(m)=1\iff\mathcal G_M(m)>0$;
if $m=2k+1$ is odd, $b_2(m)=1\iff\mathcal G_M(m)>Y_k/(ks\sqrt M)$, and $|Y_k/k|\le K$, so
$$b_2(m)=1\ \Longleftarrow\ \mathcal G_M(m)>K/(s\sqrt M),\qquad
b_2(m)=0\ \Longleftarrow\ \mathcal G_M(m)\le -K/(s\sqrt M).$$

**(b)** $\mathcal G_M(m)=\mathcal D_M(m/M)$ exactly.

**(c)** Consequently, if $\inf_{\theta\in[1-c,1]}\mathcal D_M(\theta)>3K/(s\sqrt M)$ then $b_2(m)=1$
for every integer $m\in[(1-c)M,M]$, hence $\mathrm{cost}\ge cM$.

*Proof.* (a) is (1.1) together with $|Y_k|=\bigl|\sum_{i\le k}c_i\bigr|\le kK$. (b) is (1.2). For (c):
for such $m$, $m/M\in[1-c,1]$, so
$\mathcal G_M(m)=\mathcal D_M(m/M)\ge\inf_{[1-c,1]}\mathcal D_M>3K/(s\sqrt M)$, and the threshold in
(a) is at most $K/(s\sqrt M)$; moreover $3K/(s\sqrt M)>K/(s\sqrt M)$, so $b_2(m)=1$. Hence
$C\cap[(1-c)M,M]=\varnothing$; then either $C=\varnothing$ and $\mathrm{cost}=M-1\ge cM$, or
$m^*\le(1-c)M$ and $\mathrm{cost}=M-m^*+1\ge cM+1$. $\ \blacksquare$

The constant is $3$ rather than $1$ only to leave room between "$>$" and "$\ge$"; with $K=2$,
$s\ge1/2$, the threshold is at most $12/\sqrt M$.

---

## 6. Theorem 1, proof

The uniformity over compositions is not a by-product of the proof: it is what makes Corollary 2
usable, so it is proved explicitly, by reducing the uniform statement to a statement about an
arbitrary sequence of compositions. This is where the *uniformity* of (4.1) — as opposed to the mere
fact of convergence for each fixed multiset — is used, and it is the reason no multivariate
Berry–Esseen theorem is needed.

**Proposition 6.1.** Fix $c\in(0,1)$. For every sequence of multisets $\mathfrak M_k\subset\{0,1,2\}$
with $M_k\to\infty$ and $s^2(\mathfrak M_k)\ge1/4$, and every $a_k$ uniform arrangement of
$\mathfrak M_k$,
$$\liminf_{k\to\infty}P\bigl(\mathrm{cost}(a_k)\ge c\,M_k\bigr)\ \ge\ p(c).$$

*Proposition 6.1 implies Theorem 1.* Suppose Theorem 1 fails: then for some $\eta>0$ and every $M$
there are $M'\ge M$ and $\mathfrak M'$ with $P(\mathrm{cost}\ge cM')<p(c)-\eta$. Taking $M=1,2,\dots$
produces a sequence $(M_k,\mathfrak M_k)$ with $M_k\to\infty$ violating Proposition 6.1. $\square$

*Proof of Proposition 6.1.* Fix a sequence as in the statement. By Corollary 1 we may replace each
$\mathfrak M_k$ by its centred version, so $\sum_ic_i=0$ for every $k$; then $K_k\le2$, $s_k\ge1/2$
and $K_k/s_k\le4$ for all $k$.

*Step 1 (a grid).* $\mathcal V$ is a.s. continuous, so for a nested sequence of finite grids
$G\uparrow\mathbb Q\cap[1-c,1]$ the events $\{\inf_G\mathcal V>0\}$ increase to
$\{\inf_{[1-c,1]}\mathcal V>0\}$. Choose $G$ finite, containing $1$, with
$$P\bigl(\inf_G\mathcal V>0\bigr)\ \ge\ p(c)-\eta/4 .$$
Since $\{\inf_G\mathcal V>\varepsilon\}\uparrow\{\inf_G\mathcal V>0\}$ as $\varepsilon\downarrow0$,
choose $\varepsilon>0$ with $P(\inf_G\mathcal V>\varepsilon)\ge p(c)-\eta/2$.

*Step 2 (finite-dimensional limit along the sequence).* By Lemma 3 and remark 3, applied with
$M=M_k$, the vector $(\mathcal D_{M_k}(\theta))_{\theta\in G}$ converges in distribution to
$(\mathcal V(\theta))_{\theta\in G}$: the bound (4.1) is at most $C_{\alpha}(K_k/s_k)/\sqrt{M_k}\to0$
for each Cramér–Wold direction $\alpha$, with $C_\alpha$ depending only on $G$ and $\alpha$, **not**
on $k$ or on the composition. The set $U=\{x\in\mathbb R^{|G|}:\min_rx_r>\varepsilon\}$ is open, so
Portmanteau gives
$$\liminf_{k\to\infty}P\Bigl(\inf_G\mathcal D_{M_k}>\varepsilon\Bigr)\ \ge\
P\Bigl(\inf_G\mathcal V>\varepsilon\Bigr)\ \ge\ p(c)-\eta/2 .$$
(The last probability is the one in Theorem 1 because $\mathcal V=_d-\mathcal V$, so
$P(\inf\mathcal V>0)=P(\sup\mathcal V<0)=p(c)$; this is the same sign symmetry as in the balanced
case, and it is the only place it is used.)

*Step 3 (transfer).* Since $3K_k/(s_k\sqrt{M_k})\le12/\sqrt{M_k}<\varepsilon$ for $k\ge k_1$, and
$\inf_{[1-c,1]}\mathcal D_{M_k}\le\inf_G\mathcal D_{M_k}$,
$$\{\inf_G\mathcal D_{M_k}>\varepsilon\}\ \subseteq\ \{\inf_{[1-c,1]}\mathcal D_{M_k}>\varepsilon\}
\ \subseteq\ \{\mathrm{cost}\ge cM_k\},$$
the last inclusion by Lemma 4(c). Hence $P(\mathrm{cost}\ge cM_k)\ge P(\inf_G\mathcal D_{M_k}>\varepsilon)$
for $k\ge k_1$.

*Step 4 (assembly).* Combining Steps 2 and 3,
$\liminf_kP(\mathrm{cost}\ge cM_k)\ge p(c)-\eta/2$. As $\eta$ was arbitrary,
$\liminf_kP(\mathrm{cost}\ge cM_k)\ge p(c)$. $\ \blacksquare$

**Remark 6.2 (why this is the uniformity the referee asks for).** Nothing in Steps 1–4 is required to
hold simultaneously over compositions: Step 1 chooses $(G,\varepsilon)$ once, Step 2 is a statement
about the single sequence at hand (the error bound depends on $k$ only through $M_k$), and Step 3 is
deterministic. So the argument proves the *uniform* statement by contradiction, without needing a
multivariate Berry–Esseen theorem or a rate for the event probability
$P(\inf_G\mathcal D_M>\varepsilon)$ — a rate for the one-dimensional projections suffices. This is
exactly the step at which a composition-dependent CLT error would have destroyed the conclusion, and
(4.1) is where it is ruled out.

**Where each hypothesis is used, and where it is not.**

| ingredient | used? |
|---|---|
| $\kappa_0=\kappa_2$ | **nowhere** (Lemma 2 needs $s^2>0$ only) |
| $\sum a=2h$ (bridge) | replaced by $\sum c=0$ after centring, which is an identity for every multiset |
| $a_i\in\{0,1,2\}$ | only through $K\le2$ (Lemma 3, 4) and $s\ge1/2$ |
| $s^2\ge1/4$ | only through $K/s\le4$ (1.3); any $s^2\ge s^2_{\min}>0$ works with $K/s_{\min}$ |
| the mean $\mu$ | only through Lemma 1 (exact invariance); no drift term anywhere |

---

## 7. Corollary 2 — the form in which §8 uses the bound

§8 needs $E[\mathrm{cost}(i,n-1)]\ge c_0(n-i)(1-o(1))$, where $\mathrm{cost}(i,n-1)$ is the cost of
the suffix word $a[i..n-1]$, whose multiset is random and correlated with the history.

**Corollary 2.** Let $a$ be a random tape of an exchangeable ensemble whose conditional law given the
counts is the uniform arrangement (both ensembles of the paper qualify), and let $M=n-i$. Then
conditionally on the suffix multiset $\mathfrak M_i=\{a_i,\dots,a_{n-1}\}$ the suffix word is a
uniform arrangement of $\mathfrak M_i$, so by Theorem 1, for $M\ge M_0(c,\eta)$,
$$P\bigl(\mathrm{cost}(i,n-1)\ge cM\ \big|\ \mathfrak M_i\bigr)\ \ge\ p(c)-\eta
\qquad\text{whenever } s^2(\mathfrak M_i)\ge1/4,$$
the bound being **uniform over $\mathfrak M_i$**. Consequently
$$E\bigl[\mathrm{cost}(i,n-1)\bigr]\ \ge\ cM\,\bigl(p(c)-\eta\bigr)\ -\ M\,P\bigl(s^2(\mathfrak M_i)<1/4\bigr).
\tag{7.1}$$

*Proof.* Exchangeability makes the tape a uniform arrangement of its own multiset; conditionally on
the suffix multiset, the ordering of the suffix is uniform and independent of everything else. The
bound of Theorem 1 is uniform in the composition, so it holds simultaneously for all $\mathfrak M_i$.
Finally $\mathrm{cost}\le M$ always, so the contribution of the exceptional event is at most $M$ times
its probability. $\ \blacksquare$

**Lemma 5 (the exceptional event).** Let $\mathfrak M$ have $M$ cells with $p_0=\kappa_0/M$,
$p_2=\kappa_2/M$. If $p_0\ge\rho$ and $p_2\ge\rho$ with $\rho>1/8$, then for a uniform $L$-subset $A$
($L\le M$),
$$P\bigl(s^2(A)<1/4\bigr)\ \le\ 2\exp\bigl(-2L(\rho-\tfrac18)^2\bigr).$$

*Proof.* For a multiset with proportions $(q_0,q_1,q_2)$, $s^2=q_1+4q_2-(q_1+2q_2)^2$. Suppose
$q_0\ge1/8$ and $q_2\ge1/8$. With $u=q_1$, $w=q_2$, $q_0=1-u-w\ge1/8$ gives $u\in[0,\frac78-w]$ and
$w\in[\frac18,\frac78]$; $s^2=u+4w-(u+2w)^2$ is concave in $u$, so its minimum over the interval is at
an endpoint: at $u=0$ it is $4w(1-w)\ge7/16$; at $u=\frac78-w$ it is $\frac7{64}+\frac54w-w^2$, concave
in $w$, whose endpoint values are $\frac14$ (at $w=\frac18$) and $\frac7{16}$ (at $w=\frac78$). So
$q_0,q_2\ge\frac18\Rightarrow s^2\ge\frac14$, with equality only at $(q_0,q_1,q_2)=(\frac18,\frac34,\frac18)$.
Therefore $s^2(A)<1/4$ forces $\#0(A)<L/8$ or $\#2(A)<L/8$. Each count is hypergeometric with mean
$\ge\rho L$, and Hoeffding's inequality for sampling without replacement gives
$P(\#0(A)<L/8)\le\exp(-2L(\rho-1/8)^2)$; add the two. $\ \blacksquare$

For the two-column ensemble $\kappa_0=\kappa_2=M_{\text{tape}}/2$, so $\rho=1/2$ and the bound is
$2e^{-0.28L}$; in the "uniform slice" ensemble the counts concentrate at $(n/3,n/3,n/3)$ with
exponential tails, so $\rho\ge0.3$ w.h.p. and the bound is $2e^{-0.06L}$. Either way
$M\,P(s^2<1/4)\to0$ geometrically, and (7.1) gives
$$E[\mathrm{cost}(i,n-1)]\ \ge\ c\,(n-i)\,(p(c)-o(1))\qquad (n-i\to\infty),$$
uniformly in $i$, which is exactly the input of §8 with $c_0=c\,p(c)$.

**On the two objections to (7.1) raised in REFEREE_S8 §4(b).** The prose "the suffix is
unconditioned" is still wrong (the event that epoch $i$ exists constrains the suffix), and the
argument still does not use the conditional law given the history: it bounds the *marginal*
$E[\mathrm{cost}(i,n-1)]$, and the conditioning is paid for by the $(n-i)P(A_i^c)$ term. What is new
here is that (7.1) is now valid — the comparison variable's law is covered by Theorem 1 for **every**
composition, so the random and history-correlated composition of the suffix is no longer an obstacle.
A by-product: because the bound is uniform, one also gets the conditional statement
$E[\mathrm{cost}(i,n-1)\mid\mathfrak M_i]\ge cM(p(c)-\eta)$ for every realisable $\mathfrak M_i$ with
$s^2\ge1/4$, which is stronger than needed.

**Consequence for the paper.** With Theorem 1 in place, the theorem of §8 reads: for every
$c\in(0,1)$, with $c_0=c\,p(c)$,
$$E[T]\ \ge\ \Bigl[\frac{c_0^2}{2(1+c_0)}-O(q_N)\Bigr]n^2(1-o(1)),$$
and the **only** remaining hypothesis is $q_N=P(\text{meet before epoch }N)\to0$. The sentence "the
extension of the single-epoch theorem to unbalanced multisets" can be deleted from the list of open
inputs.

---

## 8. What is proved, what is cited, what is not claimed

**Cited (two items, both standard and used only with explicit constants).**
1. The Stein-equation bounds $\sup|f'|\le1.2533$, $\sup|f''|\le2$ for $|h|\le1$ and the standard
   normal target (Stein 1986; Chen–Goldstein–Shao, Lemma 2.4). Used in Lemma 3 Step 4.
2. Hoeffding's tail inequality for sampling without replacement (Hoeffding 1963). Used in Lemma 5.
   (Any hypergeometric tail bound would do.)
Plus the Portmanteau theorem.

**Not cited, and not needed:** Hoeffding's combinatorial CLT, Serfling's inequality, the
Chentsov–Kolmogorov tightness criterion, Dudley's entropy bound, Borell–Tsirelson. The proof of
Theorem 1 uses only finite-dimensional convergence (Lemma 3) plus the deterministic transfer
(Lemma 4). This is a strict reduction in the citation load relative to REPORT_EPOCH_BOUND, where the
functional CLT was needed.

**Not claimed.**
* Convergence in distribution of $\mathrm{cost}/M$ is *not* asserted. Only the liminf lower bound,
  which is all the paper uses. (Convergence is presumably true and is visible in the numbers, but it
  would need a continuity/tightness argument for the argmax functional that is not given here.)
* No rate in $M$ for Theorem 1 beyond "for $M\ge M_0(c,\eta)$": $M_0$ is effective in principle
  (through Lemma 3 and the grid choice) but not computed.
* The constant $K/s\le4$ is what $s^2\ge1/4$ buys. Compositions with $s^2<1/4$ are *not* covered;
  empirically they behave the same (the referee's $(0,64,192)$ at $s^2=.188$ gives $.135$), but no
  proof is offered, and Lemma 5 shows they do not arise for $s^2(\text{suffix})$ in the paper.
* The odd-$m$ thresholds are handled with the crude bound $|Y_k/k|\le K=2$. This costs a factor 2 in
  the constants and nothing else; it is not claimed to be sharp.

---

## 9. Numerical verification

All runs in `probe_ologn2/` with
`$P = C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`,
`$P q13_unbalanced.py all > q13.log` (35 s). The script re-derives everything from the definition of
`cost` and never uses any of the lemmas.

**A. Translation/affine invariance, exact.** `cost(a) == cost(a+shift)` for 2000 random words and
shifts $\{+7,-5,+\tfrac14\}$: **0 violations in 6000 tests**. (Affine maps with $\alpha\ne1$ checked
implicitly by C2 below.) The drift table shows the referee's point and its resolution:

    comp (0, 8, 8)  mu=1.500  sqrt(M)(mu-1)theta/s =  1.00  2.00  3.00  4.00   (theta=.25,.5,.75,1)
    comp (8, 8, 0)  mu=0.500                         -1.00 -2.00 -3.00 -4.00
    comp (6, 2, 6)  mu=1.000                          0.00  0.00  0.00  0.00

and at $M=1024$, composition $(0,512,512)$, over 2000 tapes:

    mean of        X_M/(s sqrt M)      =  32.00   (predicted (mu-1)sqrt(M)/s = 32.00)
    mean of (X_M - 2X_{M/2})/(s sqrt M) =  0.002   (predicted 0),  sd = 0.995

i.e. the uncentred normalised mean diverges like $\sqrt M$ exactly as the referee says ($32$ at
$M=1024$, and it grows), while the normalised *contrast* has mean $0$ to three decimals and standard
deviation $1$. The cost is unchanged throughout.

**B. Exact covariance identity, every composition.** At $M=12$, 82 compositions with $s^2\ge0.05$ and
at least 30 arrangements; at $M=8$, 30. Max
$|\mathrm{Cov}_{\text{emp}}-(2(\theta\wedge\theta')-(\theta\vee\theta'))_+M/(M-1)|$:
$3.3\times10^{-15}$ ($M=8$), $7.6\times10^{-15}$ ($M=12$). Includes $s^2<1/4$ cases such as
$(0,3,9)$ at $s^2=.188$, and compositions with $\mu\ne1$. This confirms Lemma 2 has no balance
hypothesis.

**C2. Exact cost law along fixed composition families (full enumeration).** `(M/2,0,M/2)`,
`(0,M/2,M/2)`, `(M/2,M/2,0)` give **identical** laws at every $M$ — they are related by the increasing
affine maps $a\mapsto a/2+1$ and $a\mapsto a+1$, which Lemma 1 permits — a clean confirmation that
only $s^2$ (through the normalisation) can matter, and here it does not even do that. Families not
affinely equivalent to these differ at finite $M$ by an amount consistent with $O(M^{-1/2})$:

| $M$ | `(M/4,M/2,M/4)` $s^2=.5$: KS vs ref, KS$\cdot\sqrt M$ | `(3M/8,M/4,3M/8)` $s^2=.75$ | `(M/8,3M/4,M/8)` $s^2=.25$ |
|---|---|---|---|
| 8  | .071, .202 | .129, .364 | .057, .162 |
| 12 | .069, .240 | .217, .750 | — |
| 16 | .058, .230 | .080, .319 | .036, .144 |

with $E[\mathrm{cost}]/M$ moving from .17–.28 at $M=8\ldots20$ towards the Brownian .132. The parity
effect is visible and is understood exactly: at even $M$ the top test can tie, and a tie counts as
$b_2=0$; for $(M/2,0,M/2)$ at $M$ even,
$P(\mathrm{cost}=1)=\frac12\bigl(1+P(X_{M/2}=0)\bigr)$, which gives
$\frac12(1+400/924)=.7165$ at $M=12$ — the tabulated value, exactly. So the finite-$M$ bias in
$P(\mathrm{cost}=1)$ is half the tie probability at the top test, an $O(M^{-1/2})$ local-CLT quantity.

**D. Berry–Esseen, exact and Monte Carlo.** Exact Kolmogorov distance of the standardised
$\mathcal D_M(1)$ from $N(0,1)$ is .08–.19 at $M=10\ldots13$ (bound (4.1) is vacuous there, as
expected). Monte Carlo ($2\times10^4$ arrangements, $\theta=1$, i.e. the top test):

| $M$ | 64 | 128 | 256 | 512 | 1024 | 2048 |
|---|---|---|---|---|---|---|
| KS, `(M/8,3M/4,M/8)` $\mu=1,s^2=.25$ | .0991 | .0717 | .0506 | .0391 | .0275 | .0199 |
| KS, `(0,M/2,M/2)` $\mu=1.5,s^2=.25$ | .0992 | .0716 | .0525 | .0404 | .0279 | .0213 |
| KS, `(M/2,0,M/2)` $\mu=1,s^2=1$ | .0980 | .0706 | .0509 | .0412 | .0295 | .0208 |
| bound (4.1), `(M/8,3M/4,M/8)` ($K/s=2$) | .720 | .509 | .360 | .255 | .180 | .127 |
| bound (4.1), the other two ($K/s=1$) | .360 | .255 | .180 | .127 | .090 | .064 |

The three compositions — means $1.0$, $1.5$, $1.0$, variances $.25$, $.25$, $1.0$, two-point and
three-point supports — give the **same** distance to within the Monte-Carlo resolution at every $M$
(the largest discrepancy between columns is $.002$), and the bound of Lemma 3 holds in every cell.
This is the numerical form of "uniform in the composition". (The measured decay is slightly slower
than $M^{-1/2}$: KS$\cdot\sqrt M$ drifts from .56 to .97 over $M=64\ldots2048$; any $o(1)$ rate
suffices for Theorem 1.)

**E. Monte Carlo cost law, $M=256\ldots2048$, seven compositions, 4000 tapes each.** $M=2048$:

| composition | $s^2$ | $E[\mathrm{cost}]/M$ | $P(\ge.05M)$ | $P(\ge.1M)$ | $P(\ge.2M)$ | $P(\ge.3M)$ | $P(\ge.5M)$ |
|---|---|---|---|---|---|---|---|
| (256,1536,256) | .250 | .1328 | .385 | .334 | .260 | .193 | .091 |
| (512,1024,512) | .500 | .1326 | .384 | .329 | .250 | .194 | .093 |
| (1024,0,1024) | 1.000 | .1316 | .390 | .336 | .251 | .188 | .092 |
| (0,1024,1024) | .250 | .1345 | .394 | .339 | .251 | .193 | .095 |
| (1024,1024,0) | .250 | .1341 | .382 | .337 | .258 | .193 | .093 |
| (768,512,768) | .750 | .1336 | .379 | .329 | .259 | .195 | .091 |
| (682,682,684) | .667 | .1253 | .369 | .320 | .242 | .177 | .081 |
| **Brownian** | — | **.1324** | **.3833** | **.3301** | **.2534** | **.1934** | **.0924** |

Monte-Carlo standard errors: .008 on a probability, .005 on the mean. At $M=2048$ all seven
compositions agree with the Brownian values and with each other within $2$ s.e.; at $M=256$ the same
seven span .131–.146 and $P(\ge.1M)\in[.327,.353]$, the expected $O(M^{-1/2})$ spread. The Brownian
row is computed at grid $N=4096$ with $2\times10^4$ paths (its own Monte-Carlo s.e. is .0035); the
grid bias is real and negative — the same run gives .3880/.3368/.2607/.1978/.0965 at $N=1024$ and
.3848/.3328/.2569/.1966/.0955 at $N=2048$ — so the continuum $p(c)$ is a few thousandths below the
tabulated value, and this grid effect, not the composition, is the largest systematic in this
comparison. (REPORT_EPOCH_BOUND §7D quotes .3810/.3282/.2504/.1914/.0913 from $4\times10^4$ paths.)

**F. Suffix variance.** $P(s^2(\text{suffix of length }L)<1/4)$ computed exactly (multivariate
hypergeometric, $n=4096$):

| $L$ | 32 | 64 | 128 | 256 | 512 | 1024 |
|---|---|---|---|---|---|---|
| two columns $(\tfrac n2,0,\tfrac n2)$ | 2.4e-6 | 5.6e-11 | 2.7e-20 | 1.3e-39 | 7.4e-82 | 1.2e-185 |
| uniform slice $(1365,1366,1365)$ | 3.4e-3 | 6.1e-5 | 2.1e-8 | 1.6e-15 | 7.6e-31 | 1.1e-67 |

and the algebraic criterion behind Lemma 5 was checked on a grid: $\min s^2$ over
$\{p_0,p_2\ge1/8\}$ is exactly $0.2500$, attained at $(1/8,3/4,1/8)$.

**G. Suffix cost on real tapes (the exact object of Corollary 2).** $n=1024$, 1500 tapes per ensemble,
suffixes $a[i..]$:

| $i$ | 0 | 64 | 128 | 256 | 512 |
|---|---|---|---|---|---|
| two columns: $E[\mathrm{cost}]/M$ | .1408 | .1384 | .1424 | .1437 | .1371 |
| uniform slice: $E[\mathrm{cost}]/M$ | .1338 | .1386 | .1380 | .1365 | .1435 |
| two columns: $P(\ge.1M)$ | .350 | .348 | .358 | .348 | .336 |

with $P(s^2(\text{suffix})<1/4)=0$ in every cell, against a Brownian $.132$ and $.330$ (Monte-Carlo
s.e. .013 on a probability, .008 on the mean). This is the end-to-end check: the bound survives the
random, history-correlated composition of the suffix.

---

## 10. Files

```
probe_ologn2/q13_unbalanced.py    A invariance/drift  B exact covariance  C all compositions
                                  C2 fixed families   D/D2 Berry-Esseen   E0/E Brownian vs MC
                                  F suffix variance   G suffix on real tapes
probe_ologn2/q13.log              full output (35 s), quoted above
probe_ologn2/REPORT_UNBALANCED.md this file
```

No `.tex` file was touched.
