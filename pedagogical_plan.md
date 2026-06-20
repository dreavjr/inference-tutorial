# Step 1: Classical statistical test

## Pedagogical plan

Classical inference doesn't accept parameters as random variables because that would clash with the frequentist notion of probability. To bypass that blocking point, it uses the experiment as the source of variation. Each experiment is associated to a random sample, which allow statistics or estimators to have a sampling distribution even though the true parameters are still consided fixed. Instead of the *parameters*, whose uncertainty classical inference dosn't even have language to define, it's the *estimators* that will be the object of significance or confidence analyses.

### Concepts
- Experimental sample, sampling distribution, test statistic
- Tail areas and p-values
- Statistical tests and confidence intervals

### Try this
1. Reset the position of the sliders to their initial values by refreshing the page or clicking on the "reset" button.
2. Drag Δ to 0: you just made H₀ true! Resample a few dozen times or press "simulate": D̂ and p bounce; ~1 in 20 dips below 0.05 — that is α, the false-positive rate.
3. Δ = 0.5, a real but modest effect: H₀ is now false, but will rarely be rejected. The test is underpowered
); raise n to 100 and it's significant almost always — the band tightens like 1/√n while D̂ stays put. Significance measures evidence, not importance.
4. Crank σ² up: null band widens, significance melts away — a tug of war between signal (D̂) and noise (band width se = s_p√(2/n)).
5. Drop n to 3–4: band turns squat and wide, ±t·se critical ticks slide far out; raise n and the band grows so tall and tight it fades out at the top.

### Take-homes
- The p-value answers one question: "If there were truly no difference, how surprising would my test statistic be?"
- It is *not* the (unconditional) probability that H₀ is true — classical inference doesn't have language for that
- A classical hypothesis is either true or false, our ignorance doesn't make it probabilistic

 "What's the chance that these samples come from different means?" you ask.

 "That question is nonsense," classical inference says. "Either the means are the same or they are different, there's no chance involved."

 "But how do I know which is which?"

 "You don't, but tell you what: if they're the same, I'll only say they're different a small fraction (alpha) of the times you ask me."

 "But right now, for *these samples* right here, are the means different?""

"Yes!"

"What is the chance you're pulling my leg?""

"This is again nonsense. Either I'm lying or I'm saying the truth, it's not a matter of chance."

"But how do I know if, right now, you're lying?""

"You don't."

## Inferences

### Inference 1 — Pooled two-sample t-test

- **Short description:** testing whether the two means of two Gaussian populations with the same unknown variance are equal (classical pooled two-sample t-test).
- **Mathematical model:** frequentist, so *no prior or posterior*.
  - *Likelihood:* each group i.i.d. Gaussian, `x_A ~ Normal(μ_A, σ²)`, `x_B ~ Normal(μ_B, σ²)`, common unknown σ². H₀: Δ = μ_B − μ_A = 0.
  - *Sampling distribution / test statistic:* `t = (x̄_B − x̄_A) / (s_p·√(1/n + 1/m))` which, under H₀, follows a Student-t with `df = n + m − 2`. The p-value is the two-sided tail area; reject if p ≤ α.
- **Critical code** (`assets/stats.js`, `tTest`):

```js
const sp2 = ((n - 1) * svar(xa) + (m - 1) * svar(xb)) / (n + m - 2);
const dObs = mean(xb) - mean(xa);                  // observed difference of means
const se = Math.sqrt(sp2 * (1 / n + 1 / m));       // its standard error
const df = n + m - 2, t = dObs / se;
return { t, df, p: 2 * (1 - tCdf(Math.abs(t), df)), dObs, se, sp: Math.sqrt(sp2) };
```

### Inference 2 — Power / Type-II error at the true Δ

- **Short description:** computing the power (1 − β) of the above test at a hypothesized true effect Δ, via the exact noncentral-t distribution of the statistic under H₁.
- **Mathematical model:** frequentist; no prior/posterior. Under the true Δ the statistic is a noncentral t, `t = Z/√(W/df)` with `Z ~ Normal(ncp, 1)`, noncentrality `ncp = Δ/se_true`, and `W ~ χ²_df` the independent scale from the pooled-variance estimate. We reject when `|t| > t_crit`, so power is the conditional normal tail `Φ(ncp − t_c·s) + Φ(−ncp − t_c·s)` (with `s = √(W/df)`) **averaged over the χ² spread of the estimated SD** — not pretending the SD is known, which is the error a plain normal approximation makes (badly so at small n).
- **How computed:** numerically — Simpson's rule over the χ² bulk of W (`~1e-4` accurate). *(Upgraded from the earlier plain-normal approximation `Φ(ncp − t_c) + Φ(−ncp − t_c)`.)*
- **Critical code** (`assets/stats.js`, `tTestPower`):

```js
const seTrue = sigma * Math.sqrt(2 / n);
const tc = tCrit(alpha, df), ncp = delta / seTrue;
const sd = Math.sqrt(2 * df);                 // χ²_df sd; integrate over its bulk
const lo = Math.max(1e-6, df - 10 * sd), hi = df + 12 * sd;
const M = 400, h = (hi - lo) / M;
const logK = -(df / 2) * Math.log(2) - logGamma(df / 2);
const chi2pdf = (w) => Math.exp(logK + (df / 2 - 1) * Math.log(w) - w / 2);
let power = 0;
for (let i = 0; i <= M; i++) {
  const w = lo + h * i, s = Math.sqrt(w / df);
  const tail = normCdf(ncp - tc * s) + normCdf(-ncp - tc * s);
  const wt = (i === 0 || i === M) ? 1 : (i % 2 ? 4 : 2);   // Simpson weights
  power += wt * tail * chi2pdf(w);
}
power = Math.min(1, power * h / 3);
return { power, beta: 1 - power };
```

# Step 2: What we want vs. can control

## Pedagogical plan

The ideas of power and significance appear with different names: sensitivity and specificity, true positive rate and true negative rate, etc. They are formal properties of the test but do not correspond to the challenges we will face in reality. PPV and NPV correspond much more to real-life concerns ("the test tagged this as positive/negative: how often can I expect it to be right?") while power and significance answer "the sample is positive/negative: how often can I expect the test to tag it right" (which is not as useful, because the true labels are, of course, unknown: if I had them I wouldn't need a test!).

### Take-homes
- Test-bound knobs are easy (change the test).
- Reality-bound knobs are hard (changing the reality is much harder).
- Sensitivity (1 - beta) and specificity (1 − α) are properties of the test, conditioned on what is true in the world.
- PPV and NPV are properties of the real world, conditioned on the test answer. They is useful because a test is always observable, while latent reality seldom is.
- Sensitivity and specificity don't bound false positives and negatives, not even asymptotically, unless I make assumptions on the base rate.

### Try this
1. Reset the sliders to default. This is the comfortable zone textbooks describe. PPV is not exactly specificity and NPV is not exactly sensitivity, but the match is close enough to be useful.
2. Base rate down to 10%: α and power unchanged, yet PPV falls to ≈ 64% (one positive in three is false).
3. Down to 2%: PPV ≈ 25% — three "discoveries" in four are wrong, same respectable test.
4. Try to rescue PPV without touching the base rate. Does it respond more to sensitivity or specificity? Do you get why? (Hint: look at the change in the arrows!)
5. Crank total cases up to 1,000,000: every count grows 1000×, none of the four metrics budge. Sample size will not save you.


## Inferences

### Inference — Predictive values from operating characteristics + base rate (diagnostic Bayes)

- **Short description:** inferring PPV/NPV (how often a positive/negative *verdict* is correct) from sensitivity (power), specificity (1 − α), and the base rate — Bayes' theorem in its diagnostic-test form. This page has no continuous parameter and no sampling; the "inference" is the deterministic inversion on a 2×2 table.
- **Mathematical model:** Bayes on counts, not a density.
  - *Prior:* the base rate, `P(condition +) = prev`.
  - *Likelihood:* sensitivity `P(predicted + | condition +)` and specificity `P(predicted − | condition −)`.
  - *Posterior:* PPV `= TP/(TP+FP) = P(condition + | predicted +)` and NPV `= TN/(TN+FN)`.
  - *Computation:* exact / analytic (closed-form Bayes); the page realizes it as integer counts over a finite population.
- **Critical code** (`02-pvalues-and-power.html`):

```js
function confusion(prev, sens, spec, total) {
  const [P, N] = split(total, prev);     // prior: base-rate split into condition ±
  const [TP, FN] = split(P, sens);       // likelihood: sensitivity within positives
  const [TN, FP] = split(N, spec);       // likelihood: specificity within negatives
  return { P, N, TP, FN, FP, TN, PP: TP + FP, PN: FN + TN };
}
// posterior (predictive values): PPV = pct(c.TP, c.PP), NPV = pct(c.TN, c.PN)
```

# Step 3: Bayesian inference's pipeline: prior -> lihelihood -> posterior -> da capo

## Pedagogical plan

Full Bayeasian treatment: representation of uncertainty, marginalization of nuisance variables (the variance, in this case). Conclusion is not a single value or hard interval: the posterior _is_ the conclusion.

### Concepts
- Prior, likelihood, posterior

### Take-homes
- Likelihoods are not normalized, nor the product prior x likelihood. Normalizing is _the_ main challenge of Bayesian inference
- Priors should reflect domain knowledge and be lax enough to allow the model to "change opinions", but making them _too_ lax leads to missed opportunities

## Current material

### Intro text
Names step 3's sleight of hand: the flat prior, which declared every effect size (Δ = 0.1, 3, a million) equally believable — not humility but a strong claim. Makes the prior explicit via Bayes as multiplication: posterior ~ prior × likelihood. The unknowns come as a pair — effect Δ and noise σ² — so the prior is a distribution over that plane; the conjugate choice is the Normal-inverse-gamma: Δ | σ²_d ~ Normal(m₀, σ²_d/κ₀) and σ²_d ~ Inv-Gamma(a₀, b₀). Three panels draw prior (purple), likelihood (terracotta), posterior (dark) as filled equal-volume contours over the (Δ, σ) plane (σ as standard deviation, both axes same scale), with a green cross at the true (Δ, σ). Below: the data and the marginal over Δ as the three 1-D shadows. A simplification keeps the algebra on one screen: the model is fed the n per-pair differences dᵢ = x_{B,i} − x_{A,i}, each Normal(Δ, σ²_d) with σ²_d = 2σ².

### Concept table
- m₀ → where you expect Δ to be, before data → a slider
- κ₀ → how strongly, in pseudo-observations: your hunch is worth κ₀ imaginary data points → a slider (try the extremes)
- a₀, b₀ → the same two roles, for the variance σ² → b₀ is a slider (where you expect σ²); a₀ held fixed and gentle

### Try this
1. Defaults (m₀ = 0, κ₀ = 1, n = 12): posterior teardrop sits between prior and likelihood, much closer to the likelihood (12 real obs outvote 1 imaginary); the marginal's dark curve peaks just short of the terracotta — a whisper of skepticism.
2. Go lax: κ₀ → 0.1 — prior flattens to a vast sheet, posterior becomes the likelihood (you've rebuilt step 3); resample at n = 12 and the posterior chases noise across the plane.
3. Skeptical but fair: m₀ = 0, κ₀ ≈ 10 — posterior pulled toward 0 (slightly biased when truth is 0.8) but resampling barely moves it; with small n, trading bias for stability is a bargain — what a prior buys.
4. Tight and wrong: m₀ = −1, κ₀ = 100 — 100 imaginary obs outvote 12 real, posterior parks near −1, green cross far outside the teardrop; raise n to 100 and the data claw back. Tight priors are loud.
5. Crank n to 100 and sweep κ₀ from 0.1 to 10: posterior barely budges — enough data washes out any reasonable prior. Priors matter when data are scarce, exactly when "non-informative" is most dangerous.

### Take-homes
- A flat prior is not absence of opinion; it is the opinion that preposterous effect sizes are as credible as modest ones — and with little data that leaks into the posterior as width and wobble.
- κ₀ prices your conviction in the only honest currency: data. Ask "how many observations is my hunch worth?" — both 0.1 (abandons to noise) and 100 (deaf to evidence) can be wrong.
- When n is small, weakly-informative beats non-informative; when n is large the argument dissolves on its own — so the cases where a lax prior "plays it safe" are precisely the cases where it doesn't.
- The prior is a dial to disclose, not a sin to hide: state it, vary it, show the posterior barely moves — that is what a robust claim looks like.

## Inferences

### Inference — Bayesian joint posterior over (Δ, σ²_d) via a Normal-inverse-gamma conjugate prior

- **Short description:** jointly estimating the effect Δ and the noise σ² from the per-pair differences, with a conjugate Normal-inverse-gamma prior, then reporting the marginal posterior over Δ (the variance integrated out as a nuisance).
- **Mathematical model:**
  - *Prior (Normal-inverse-gamma):* `Δ | σ²_d ~ Normal(m₀, σ²_d/κ₀)` and `σ²_d ~ Inv-Gamma(a₀, b₀)`. Here `σ²_d = 2σ²` (the variance of a difference of two groups). a₀ is held fixed and gentle; b₀ tracks the σ²₀ slider; m₀, κ₀ are sliders.
  - *Likelihood:* the n per-pair differences `dᵢ ~ Normal(Δ, σ²_d)`.
  - *Posterior:* Normal-inverse-gamma with updated `(mₙ, κₙ, aₙ, bₙ)`; the *marginal over Δ* is a located-scaled Student-t with `df = 2aₙ`, `loc = mₙ`, `scale = √(bₙ/(aₙ·κₙ))`.
- **How the posterior was computed:** *conjugate prior–likelihood → closed-form (analytic) update* (`nigUpdate`). The 2-D (Δ, σ) surfaces are evaluated on a grid only for the equal-volume heatmap drawing (with a Jacobian for the σ² → σ reparametrization); the inference itself is the analytic conjugate update, and the Δ marginal is the analytic Student-t.
- **Critical code** (`03-bayesian-approach.html`, `render`):

```js
const { dbar, ssd } = drawSample(pool, { delta, sigma, n });
const { kn, mn, an, bn } = nigUpdate(m0, k0, A0, B0, n, dbar, ssd);   // conjugate update
...
const pri = nigMarginalDelta(m0, k0, A0, B0);
const pos = nigMarginalDelta(mn, kn, an, bn);                          // Student-t marginal over Δ
```

  with the closed-form update and marginal in `assets/stats.js`:

```js
function nigUpdate(m0, k0, a0, b0, n, dbar, ssd) {
  const kn = k0 + n;
  const mn = (k0 * m0 + n * dbar) / kn;
  const an = a0 + n / 2;
  const bn = b0 + ssd / 2 + k0 * n * (dbar - m0) ** 2 / (2 * kn);
  return { kn, mn, an, bn };
}
function nigMarginalDelta(m, k, a, b) {
  return { loc: m, scale: Math.sqrt(b / (a * k)), df: 2 * a };   // df = 2a
}
```

# Step 4: Everyone is entitled to their opinion (but they'll lose money if they disagree with me)

## Pedagogical plan

The Bayesian treatment is the logical way to think under uncertainty. In other words, given what you know, the probabilities it outputs are the  unique fully-coherent conclusion. Even with the cost model known by both, a decision based on classical inference has a free parameter (alpha) that must be chosen externally to the decision apparatus, while the Bayesian decision apparatus can marginalize over all nuisances to make a decision that maximizes expected reward. The consequence is that under a system of bets, the Bayesian decision procedure will inevitably make better bets (in expectation), provided that you don't start from too wrong _and_ too strong assumptions.

## Current material

### Intro text
Step 4 handed us a posterior — a whole distribution over Δ — but a distribution is not yet a decision; sooner or later you must act (ship the drug, set the dose, choose B over A). Bayesian decision theory closes the gap with one rule: act to maximise expected utility under your posterior. It needs the posterior plus a reward model U(x) = 1 − |Δ − x|^p − c·|Δ − x|: land on the truth and earn a full 1; miss and two penalties bite, a polynomial |Δ − x|^p and a linear c·|Δ − x|. For positive c the reward crosses zero at a break-even error ε (the largest profitable miss); raising p oddly pushes ε outward for small misses; with c = 0, ε = 1. The bottom row poses the sharp question — given my belief, should I bet at all, and on what value? — both readings betting on the peak Δ̂ but deciding whether to bet differently. Left: the honest posterior, where the Bayesian rule averages reward against the whole posterior, E[U] = ∫ reward(Δ − Δ̂)·P(Δ) dΔ, and bets only when that is positive (the gold lens is the overlap). Right: the classical reading, a (1 − α) confidence interval drawn as a flat error bar (MLE Δ̂ centred, margin E either way); with no distribution to average over the rule turns worst-case — bet iff E ≤ ε, so even the interval's worst boundary still pays. Each panel shows an estimated pill (the rule's forecast — posterior's expected reward, or CI's min…max) turning green to bet / grey to fold, and an actual pill (the payoff once hidden Δ is revealed, or zero for a refused bet) green for a win / red for a loss / grey for a no-play. The gap between the pills is the price of not knowing the truth.

### Concept table
- p → the polynomial degree, how the penalty curves: a gentle bowl at p = 1, a flat top with steep walls as p climbs → a slider (1…10)
- c → the linear penalty, a flat charge per unit of error; raising it tilts the whole reward down → a slider (−10…10)
- Δ → where the reward is centred (the truth itself), drawn only because we simulate → the green line

### Try this
1. Defaults: with c = 0 the only penalty is polynomial, so break-even ε = 1 and the CI fits inside it; both rules bet and win — the paradigms agree.
2. Drop n to 3–4: the CI balloons past ε = 1, its worst-case boundary turns negative, the classical rule refuses (banks 0); the compact posterior still expects positive return — it bets and wins. The Bayesian edge is sharpest when data are scarce.
3. Crank n to 100: both beliefs sharpen, both bet, actual payoffs converge near peak reward 1 — with enough data the rules agree again.
4. Raise c (linear penalty): break-even ε pulls in below 1; first the wide-interval classical rule folds, push further and even the posterior's expectation drops below zero so the Bayesian folds too. (With c > 0, a higher p widens ε, nudging both back toward betting.)
5. Tighten the prior (κ₀ ↑) and shift m₀ off the truth: the posterior's bet pulls away from the interval's midpoint — watch a confident-but-wrong prior keep its estimated pill green and betting while the actual pill slides red: overconfidence made quantitative.

### Take-homes
- A posterior is an input to a decision, not the decision itself; pair it with a reward model and act on expected utility — the whole of Bayesian decision theory.
- A confidence interval is not a probability distribution; reading it as "uniform over plausible values" is a specific, usually worse, belief — and a decision made on it earns less.
- The shape the posterior carries (mass near the estimate, thinning into the tails) is what a uniform block discards and what a reward centred on the truth pays you to keep.
- Same data, same interval width — yet a different decision, because of how honestly the uncertainty inside the interval is shaped.
- Two numbers, never one: what a method thinks it will earn and what it actually earns. A trustworthy rule keeps the pair close; a tight wrong prior, or a mis-shaped interval, pries them apart.

=> Bibliography
- Accessible resources
- Practical materials and tools
- Academic material

## Inferences

### Inference — Bayesian decision by expected utility under the posterior (vs. the classical CI worst-case rule)

- **Short description:** turning the step-4 posterior over Δ into an *action* (bet / fold) by maximizing expected utility; the inferential object is the expected reward `E[U] = ∫ reward(Δ − Δ̂)·P(Δ) dΔ` averaged against the posterior, contrasted with a classical (1−α) confidence-interval rule read worst-case.
- **Mathematical model:**
  - *Prior / Likelihood / Posterior:* identical to step 4 — the same Normal-inverse-gamma conjugate posterior, marginal over Δ a Student-t.
  - *Reward (utility):* `U(x) = 1 − |Δ − x|^p − c·|Δ − x|`, clamped at −1; sliders for degree p and linear penalty c.
  - *Bayesian rule:* bet iff `E[U] > 0`, betting on the posterior mode mₙ.
  - *Classical rule:* take the (1−α) CI `d̄ ± t·se`; with no distribution to average over, act worst-case — bet iff the reward at the interval edge stays ≥ 0.
- **How computed:** the posterior is the same analytic conjugate update (`nigUpdate` / `nigMarginalDelta`); the *decision* expectation is then obtained by **numerical integration** (trapezoid rule, M = 4000 over [−12, 12]) of posterior × reward — not a closed form, since the reward is non-conjugate.
- **Critical code** (`04-the-decision.html`, `computeState` + helpers):

```js
const { kn, mn, an, bn } = nigUpdate(m0, k0, A0, B0, n, dbar, ssd);
const pos = nigMarginalDelta(mn, kn, an, bn);
...
const rew = (x, ctr) => rewardClamped(x, ctr, p, c);   // U = max(-1, 1 - |x-ctr|^p - c|x-ctr|)
// Bayesian rule: expected reward of betting on the posterior mode
const erPost = integrate(x => post(x) * rew(x, modePost));
const betPost = erPost > 0;
// CI rule: worst-case reward at the interval edge
const hw = tCrit(alpha, dfC) * se;
const ciLo = dbar - hw, ciHi = dbar + hw;
const ciEdgeR = rew(ciLo, modeCI);
const erCI = Math.min(1, ciEdgeR), betCI = erCI >= 0;
```

  with the numerical integrator:

```js
function integrate(f) {                 // ∫f dx by the trapezoid rule
  const A = -12, B = 12, M = 4000, h = (B - A) / M;
  let s = 0;
  for (let i = 0; i <= M; i++) s += ((i === 0 || i === M) ? 0.5 : 1) * f(A + h * i);
  return s * h;
}
```

# Step 5: A difference in style?

## Pedagogical plan

Bayesian probabilities are not tied to frequencies: they represent what you must conclude when you are reasoning under uncertainty. A probability distribution represent how much you can locate the value of a random variable given all that you know. This small change changes everything, because now parameters don't have to be fixed: their distribution represent what you know about them, not frequencies that depend on repetition.

### Concepts
- Bayesian probabilities break the asymmetry between parameters and data: everything is a random variable with a distribution
- Inference by Bayes' rule: posterior is proportional to prior x likelihood. The posterior represents how much I expect the parameter's value to be in certain regions of the parameter space
- Priors don't have to be informative

Myths to dismiss:
- Bayesian probability is all about priors => you can often use non-informative or little-informative priors (although that is not necessarily a good idea for most problems)
- I can use Bayesian inference to mimick hard-threshold classical tests and still reap the benefits => if you don't add informative priors ("let the data speak for itself" zealousness) and force sharp decisions instead of exploit the uncertainty, Bayesian inference can give exactly the same results as classical inference

### Take-homes
- This is not the Bayesian treatment: (1) I'm using maximum a posteriori with an uninformative prior = same as maximum likelihood. Credible intervals force a binarization of the decision that goes against the Bayesian ethos.

## Current material

## Pedagogical plan

### Intro text
Step 1 asked only "if H₀ were true, how surprising is our t?" and answers with a yes/no stamp — never telling us how large the effect is or how sure we are. The Bayesian move needs no new machinery, just a change of reading: take the same t-distribution and re-centre it on the observed difference instead of on 0, letting it describe the whole range of plausible effect sizes (with flat priors this is exactly the posterior for the mean difference — a shifted, scaled t). Both panels share step 1's square scale; the right panel is the distribution over Δ. Because it is a genuine density, a very certain posterior becomes a tall thin spike that fades out at the top of the box. Hovering the Δ curve maps any candidate effect size to the two population means it implies.

### Concept table
- the curve's peak → the observed difference Δ̂ = x̄₂ − x̄₁, our best estimate
- the curve's width → our uncertainty, wide when data are few or noisy
- the curve's height → it's a real density (area = 1), so concentrating ⇒ taller
- the dashed line at 0 → the null; if it sits out in the tail, the data doubt it
- the green dashed Δ → the real effect, unknown in life, shown here because we simulate

### Try this
1. Hover slowly across the Δ curve: the two implied means slide apart and back — the curve shows which separations the data consider plausible.
2. Is 0 under the curve? With defaults, barely (the Bayesian echo of "significant"); drag Δ toward 0 and the curve drifts until it straddles 0 and "no effect" becomes most plausible.
3. Raise n: the curve concentrates (taller, narrower) around the truth — certainty about size, not just a smaller p; far enough and the spike fades out the top. This is power from the inside.
4. Raise σ²: the curve fattens and sinks — noise widens the band of indistinguishable effect sizes.
5. "New sample" a few times: the curve jumps (its centre Δ̂ is the noisy observed effect) but usually keeps true Δ in its bulk; each study hands you one curve and you never see the green line.

### Take-homes
- A p-value returns a verdict; this returns an estimate with honest error bars — same t-distribution read the other way round.
- The shape is the message: "significant" can hide a curve spanning trivial to enormous effects; "not significant" can hide a curve tight around zero. The stamp throws that away; the distribution keeps it.
- One sleight of hand remains — "with flat priors" was itself a choice, made explicit in the next step.

## Inferences

### Inference — Bayesian posterior over the mean difference Δ, with flat priors

- **Short description:** estimating the effect size Δ = μ_B − μ_A as a full posterior distribution under flat (non-informative) priors — which yields a shifted, scaled Student-t, i.e. step 1's t-distribution re-centred on the observed difference.
- **Mathematical model:**
  - *Likelihood:* Gaussian, as in step 1 (`dᵢ ~ Normal(Δ, σ²_d)`).
  - *Prior:* flat (improper) on Δ and on log σ — "let the data speak."
  - *Posterior over Δ* (σ² marginalized out): a located-scaled Student-t centred at `Δ̂ = x̄_B − x̄_A`, scale = the t-test's standard error `se`, `df = n + m − 2`.
- **How the posterior was computed:** *analytic* — the flat-prior result coincides with the classical sampling object, so the page re-uses `tTest`'s `(dObs, se, df)`; no numerical integration. The (1−α) HDI is exactly the equal-tailed band `Δ̂ ± t·se`.
- **Critical code** (`05-bayesian-t-test.html`, `render` + `posteriorSvg`):

```js
const { df, dObs, se } = tTest(xa, xb);
...
const dens = ds.map(d => tPdf((d - dObs) / se, df) / se);   // posterior density over Δ
const tc = tCrit(alpha, df);
const hdiLo = Math.max(dObs - tc * se, -HALF), hdiHi = Math.min(dObs + tc * se, HALF);
```

# Step 6 (do-it-yourself appendix): Bayes on binary data

## Role and framing

Steps 1–5 are the live session (≈1 hour, adult software developers). Step 6 is **not** a
sixth live step — it is a take-home coding exercise, deliberately left a little less
finished than the rest. The reasons we settled on this:

- Step 5 carries the thesis of the whole tutorial (Making a decision apparatus
  marginalizes nuisances and maximizes expected reward; classical inference leaves α as a
  free parameter chosen outside the apparatus). It must stay live and must land on the
  continuous model the learner has lived in for four steps — one new idea (decision) on
  established ground, not two new ideas at once. So we do **not** cut Step 5 to make room.
- The continuous thread (1→3→4→5) is a complete spine. Step 2 is the only binary step and
  the only one that never gets a Bayesian answer; its closing line even promises "Bayes'
  theorem — where we go next" and then the thread is dropped. Step 6 is that owed answer.
- For this audience the strongest capstone is an **inversion**: for five steps they dragged
  sliders on a machine someone else built; here they build the machine behind the slider.
  "You've been playing with the posterior — now compute one."
- The deliberate incompleteness *is* the assignment: the scaffold and the questions are
###   fully built, but the morals and part of the concept table
  fills in. Transfer (apply the pipeline to new data) is exactly the move that consolidates
  it.

## Pedagogical plan

The exercise re-runs the entire pipeline — prior → likelihood → posterior → expected
utility → act — on a binary classifier with a single unknown: its success rate *p*. It is a
Jupyter-style notebook (an in-browser, client-side Python engine to be wired in later;
the current page is a non-functional mock-up so we can settle the design first).

Why binary data is the right place to *code* this rather than slider it:

- *p* is a probability, so its entire universe is [0, 1] — the prior, likelihood and
  posterior are just arrays over a fine grid. No tails to truncate, no axis bounds to guess
  (the thing that made Step 4's (Δ, σ) plane fiddly). The finite fixed support makes this
  the *most* concrete of the six, not the least.
- Beta is conjugate to the Binomial, so the exact posterior `Beta(α+k, β+n−k)` is known on
  paper. The grid is a *numerical* answer to a problem that also has an *analytic* one —
  which is not busywork but the only honest way to learn a numerical method: practise it
  where the answer key exists, then carry it where it doesn't (cf. learning numerical
  integration on ∫x² dx first).

Concepts (each a callback, made tactile by code):
- The prior **is** Step 2's base rate, written as `Beta(α, β)` — a few imaginary cases,
  priced in pseudo-observations (the pseudo-count α + β echoes Step 4's κ₀).
- Normalization is "the main challenge of Bayesian inference" (Step 4) shrunk to one
  Riemann sum: divide the curve by its own area.
- The decision is Step 5's rule on new data: pair the posterior with a reward, act on
  expected utility.

### Take-homes
- The same four boxes (prior, likelihood, posterior, decision) are, on binary data, ten
  lines of array arithmetic.
- Marginalization vs expectation: there is only one unknown, so what the learner codes in
  the decision step is an **expectation** over *p*, not a marginalization of a nuisance.
- The grid that solves this in one line is hopeless in higher dimensions (N cells per axis →
  Nᵈ cells); that wall is the whole reason MCMC and variational inference exist — the
  forward-pointer past the course.

Decisions about the model (where we deviated from first sketches):
- **Notation:** `Beta(alpha, beta)` for the prior; `n, p` for the Binomial (so the rate is
  `p`, gridded over [0, 1]); data is `(n, k)`.
- **No frequentist component on this page** (no Wald interval, no CI-vs-posterior contrast).
  The only backward link to Step 2 is the conceptual base-rate→prior bridge.
- **Reward simplified to a bet/fold** with a `GAIN` on a true positive and a `COST` on a
  false positive; folding is a guaranteed 0. We chose an FP cost rather than the originally
  suggested FN penalty because with FP = TN = 0 the decision degenerates ("always act,"
  no interior threshold). Consequence to keep honest in the text: this reward is linear in
  *p*, so the decision *also* has a closed form (`GAIN·E[p] − COST·(1−E[p])`, with
  `E[p] = (α+k)/(α+β+n)`). That reverses the earlier "the decision needs the grid" beat — so
  the page now frames the grid sum as the *general* recipe (works for any reward) and uses
  the closed form as a *second* answer-key check. The "hard part" message rests on
  normalization (cell 3) and the curse of dimensionality (cell 7).
  - *Open question still on the table:* keep the FP-cost version (clean interior threshold,
    closed-form decision), honor the literal FN penalty (decision degenerates to "always
    bet"), or go to the full three-cell model (TP +, FP −, FN −) for a richer threshold.

## Current material

A notebook mock-up (page-local CSS for the cell vocabulary, to be promoted to shared CSS
once the design settles and the engine is wired). Cells alternate **provided** (read-only
scaffolding / answer keys) and **your code** (signature + docstring contract + `# TODO`);
the **Run** buttons are intentionally inert. A callout at the top states plainly that the
page is deliberately unfinished and the engine is not connected yet.

### Intro text
Frames Step 6 as the Bayesian answer Step 2 never got, and as the
capstone inversion (build the machine, don't drag its sliders). States the pipeline on new
data and the two reasons binary data is the right place to code it (finite support [0, 1];
Beta-Binomial conjugacy as an answer key).

The seven cells:
0. **Setup** (provided): the grid `p = linspace(0,1,N)`, `dp`; the data `n, k`; the prior
   `alpha, beta`.
1. **Prior** (you code): `prior_kernel(p, alpha, beta)` → un-normalized Beta kernel.
2. **Likelihood** (you code): `likelihood(p, n, k)` → `p**k * (1-p)**(n-k)`, the
   n-choose-k constant dropped (washes out in normalization).
3. **Normalize** (you code): `normalize(density, dp)` — the one hard line — then
   `posterior = normalize(prior_kernel * likelihood, dp)`.
4. **Conjugate check** (provided): overlay `Beta(α+k, β+n−k)`; numeric must match analytic
   (if not, the bug is in `normalize`).
5. **Belief → act** (you code): `GAIN`/`COST`, `expected_utility_act(p, posterior, dp)` as a
   weighted sum, `decide(...)` returning bet/fold; output checks against the closed form.
6. **Why this was a toy** (read-me): the curse-of-dimensionality thought experiment.

Closing prose:
- "Questions to answer (once it runs)" — morals-as-questions: the flat prior `(1,1)` and its
  hidden bet; a tight wrong prior `(40,4)` vs more data; a rare-success low-mean prior
  `(2,40)`; raising `COST` until `decide` folds and tying break-even to `COST/(GAIN+COST)`
  and Step 5.
- "The moral (write your own)" — three `TODO` prompts for the learner.
### - Concept table
  `decide`→Step 5, with a deliberate `TODO` row.

Not yet done: wire the in-browser Python engine and the plot helpers (`plot_three`,
`plot_overlay`, the E[U]/bet-fold output); promote the notebook CSS to shared; decide whether
to link the page from `index.html` (currently unlinked so it doesn't read as a sixth live
step); resolve the reward-model open question above.

## Inferences

(Reference implementation: `bayesian-decision-do-it-yourself.ipynb` — the notebook with answers
filled in; the page `do-it-yourself.html` is the learner-facing mock-up.)

### Inference — Beta-Binomial posterior for a classifier's success rate p, computed numerically on a grid

- **Short description:** estimating a binary classifier's success rate p from `(n, k)` (k successes in n cases) using a Beta prior conjugate to the Binomial likelihood; the posterior is computed *numerically* on a [0, 1] grid (normalization by a Riemann sum) and cross-checked against the analytic `Beta(α+k, β+n−k)`.
- **Mathematical model:**
  - *Prior:* `p ~ Beta(α, β)`, kernel `p^(α−1)·(1−p)^(β−1)` — step 2's base rate, now as a density; pseudo-count α+β echoes step 4's κ₀.
  - *Likelihood:* `k ~ Binomial(n, p)`, kernel `p^k·(1−p)^(n−k)` (the n-choose-k constant dropped — it washes out in normalization).
  - *Posterior:* `Beta(α+k, β+n−k)`.
- **How computed:** *numerical* — grid `p = linspace(0,1,N)`, multiply prior kernel × likelihood, then normalize by dividing by its own Riemann-sum area (`density.sum()·dp`). The conjugacy supplies an *analytic* answer key (`scipy.stats.beta`) that the numerical posterior is asserted to match to 1e-3.
- **Critical code** (`bayesian-decision-do-it-yourself.ipynb`, cells 1–4):

```python
def prior_kernel(p, alpha, beta):
    return p**(alpha - 1) * (1 - p)**(beta - 1)        # un-normalised Beta

def likelihood(p, n, k):
    return p**k * (1 - p)**(n - k)                     # Binomial kernel in p

def normalize(density, dp):
    return density / (density.sum() * dp)              # the one hard line

posterior  = normalize(prior_kernel(p, alpha, beta) * likelihood(p, n, k), dp)
post_exact = beta_dist(alpha + k, beta + n - k).pdf(p) # conjugate answer key
assert np.allclose(posterior, post_exact, atol=1e-3)
```

### Inference — Bet/fold decision by expected utility over the posterior

- **Short description:** the step-5 decision rule on binary data — pair the posterior over p with a reward (GAIN on a true positive, COST on a false positive) and act iff the expected utility of acting is positive.
- **Mathematical model:** per-case payoff of acting is `GAIN·p − COST·(1−p)`; `E[U(act)] = ∫ (GAIN·p − COST·(1−p))·posterior(p) dp`, a single weighted sum over the grid. Bet iff `E[U] > 0`. Because the reward is *linear* in p, there is also a closed form using only the posterior mean `E[p] = (α+k)/(α+β+n)`, with break-even at `COST/(GAIN+COST)`.
- **How computed:** numerical expectation over the grid posterior (a Riemann sum), with the linear closed form as a second answer-key check. Note this is an *expectation* over the single unknown p, not a *marginalization* of a nuisance (there is only one parameter).
- **Critical code** (`bayesian-decision-do-it-yourself.ipynb`, cell 5):

```python
def expected_utility_act(p, posterior, dp):
    payoff = GAIN * p - COST * (1 - p)
    return np.sum(payoff * posterior) * dp

def decide(p, posterior, dp):
    value = expected_utility_act(p, posterior, dp)
    return value, ("bet" if value > 0 else "fold")

E_p         = (alpha + k) / (alpha + beta + n)         # closed-form check
closed_form = GAIN * E_p - COST * (1 - E_p)
```
