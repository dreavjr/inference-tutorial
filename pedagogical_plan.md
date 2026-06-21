# Step 1: How classical tests work

## Main argument

Classical inference doesn't accept parameters as random variables because that would clash with the frequentist notion of probability. To bypass that blocking point, it uses the experiment as the source of variation. Each experiment is associated to a random sample, which allow statistics or estimators to have a sampling distribution even though the true parameters are still considered fixed. Instead of the *parameters*, whose uncertainty classical inference doesn't even have language to define, it's the *estimators* that will be the object of significance or confidence analyses.

### Concepts
- Experimental sample, sampling distribution, test statistic
- Tail areas and p-values
- Statistical tests and confidence intervals

### Try this
1. Play with the sliders. Note what happens with reality (plot 1) and the test (plot 2) in different scenarios.
2. For a fixed reality, how does the test respond to different sample sizes and significances (α)?
3. Run a simulation. Do the observed mismatched significances (dotted lines) correspond to the expected rates α and β shown in plot 2?
4. Reset the position of the sliders to their initial values by refreshing the page or clicking on the "reset" button.
5. Drag Δ to 0: you just made H₀ true! Resample a few dozen times or simulate: D̂ and p bounce; ~1 in 20 dips below 0.05 — that is α, the false-positive rate.
6. Set Δ = 0.5, a real but very modest effect: H₀ is now false, but will rarely be rejected. The test is underpowered.
7. Raise n to 100 and it's significant almost always. Significance measures evidence for the difference, not its magnitude.
7. Crank σ² up: the null-hypothesis band widens, significance melts away — a tug of war between signal and noise.
8. Drop n to 3–4: band turns squat and wide, ±t·se critical ticks slide far out; raise n and the band grows so tall and tight it fades out at the top.

### Take-homes
- The p-value answers one question: "If there were truly no difference, how surprising would my test statistic be?"
- It is *not* the (unconditional) probability that H₀ is true — classical inference doesn't have language for that
- A classical hypothesis is either true or false, our ignorance doesn't make it probabilistic

### Dialog
 "What's the chance that these samples come from different means?" you ask.

 "That question is nonsense," classical inference says. "Either the means are the same or they are different, there's no chance involved."

 "But how do I know which is which?"

 "You don't, but I'll tell you what: *if* they're the same, I'll only say they're different a small fraction α of the times you ask me."

 "But right now, for *these samples* right here, are the means different?""

"Yes!"

"What is the chance you're pulling my leg?""

"That is again nonsense. Either I'm lying or I'm saying the truth, it's not a matter of chance."

"But how do I know if, right now, you're lying?""

"You don't."

### Links
https://en.wikipedia.org/wiki/Null_hypothesis_significance_testing
https://en.wikipedia.org/wiki/Statistical_significance
https://en.wikipedia.org/wiki/Power_(statistics)
https://en.wikipedia.org/wiki/P-value
https://en.wikipedia.org/wiki/Student%27s_t-test

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

## TODO

- Add the MLE estimator on the plot?
- change \hat(D) -> \hat(Δ)

# Step 2: The world vs. the test

## Main argument

The ideas of power and significance appear with different names: sensitivity and specificity, true positive rate and true negative rate, etc. They are formal properties of the test but do not necessarily correspond to the challenges we will face in reality. PPV and NPV correspond much more to real-life concerns.

Frequentist and classical analysis very often converge to the same answers for a large number of samples, but on this particular issue, increasing sampling size won't work: the four metrics are proportions and insensitive (in expectation) to sample size.

### Try this
1. Play with the sliders and observe how the contingency table (plot 1) and the fork arrows (plot 2) change.
2. Reset the sliders to default. This is the "textbook zone" (significance=95%, power=80%). PPV and NPV do not match sensitivity and specificity exactly, but the difference is moderate enough to be still useful.
3. Slide base rate down to 10%. Significance and power are unchanged, but PPV falls to ≈ 64% (one positive in three is false).
4. Slide base rate up to 90%. PPV is great, but now NPV falls to ≈ 34% (two negatives in three are false).
5. Play with the base rates. Can you find one that makes PPV ≈ Significance and NPV ≈ Specificity? What about vice-versa? Is there something special about those base rates?
6. Reset the sliders and bring base rate down to 2%. Try to rescue PPV without touching the base rate. Does it respond more to sensitivity or specificity? Do you get why? (Hint: look at the change in the arrows!)
7. Pick any base rate and crank total cases up to 1,000,000. Every count grows 1000×, but the four metrics don't budge. Sample size will not save you this time.

### Take-homes
- Sensitivity and specificity are properties of the test, while PPV and NPV are properties of both the test and the reality.
- Sensitivity (1 - beta) and specificity (1 − α) are conditioned on what is true in the world.
- PPV and NPV are conditioned on the test answer. They are more useful because test results are observable, while the latent truth often isn't.
- Sensitivity and specificity don't bound false positives and negatives, not even asymptotically, unless you make assumptions about the base rate.

### Dialog

"I can't afford too many false positives right now," you say. "Can you keep those reasonable?"

"Well," hesitates the test, "I can promise that if the sample is negative I'll rarely say it's positive."

"What's the catch? Will you give me lots of false negatives?"

"Not necessarily. I can also promise that if the sample is positive I'll rarely say it's negative."

"Good! That means I can rely on you, right?"

"That depends. What's the base rate? I mean, how many positive cases are there for each negative?"

"I have no idea, that's exactly why I need you."

"If you can't give me the base rate, I can't guarantee being reliable. But at least you can rely on me when I say something is positive. Or maybe when I say something is negative."

"How do I know which is which?""

"You don't, that depends on the base rate."

### Links
https://en.wikipedia.org/wiki/Confusion_matrix#Table_of_confusion
https://en.wikipedia.org/wiki/Sensitivity_and_specificity
https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values
https://en.wikipedia.org/wiki/Base_rate_fallacy

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

# Step 3: The Bayesian approach

## Main argument

Bayesian inference brings a profound conceptual change: probabilities are no longer just about repeatable measurements, they are about degrees of ignorance. Now, both data *and parameters* can be random variables with probability distributions.

This brings a huge conceptual simplification. I no longer need to conceive a different test statistic for every kind of hypothesis I want to test. I only need a generative model of the data from the parameters (that's the likelihood) and whatever prior reflects my opinion on the parameters before seeing the data. In the domain of the parameters I multiply the prior by the likelihood to get the unnormalized posterior. Then I marginalize (integrate) all the parameters I don't care about.


### Concepts
- Prior
- Likelihood
- Posterior

### Take-homes
- Likelihoods are not normalized, nor the product prior x likelihood. Normalizing is _the_ main challenge of Bayesian inference
- Priors should reflect domain knowledge and be lax enough to allow the model to "change opinions", but making them _too_ lax leads to missed opportunities

### Try this
1. Play with the sliders. Which plots respond to changes on the samples and which respond to changes on the prior.
2. Reset the sliders to default values.
3. Slide the prior weight all the way down. How does the posterior respond to changes on reality? Or about resampling?
3. Slide the prior weight all the way up. How does the posterior now respond to reality changes and resampling?
4. Set the prior center very wrong (Δ and σ² very different than the prior center for them).
5. Slide the prior weight all the way down and play with resampling and sample sizes. How does the posterior respond? Set the prior weight all the way up and try the same.
6. Suppose you had to bet with a colleague doing this same exercise in this same page on the parameters of their populations. They will show you their samples, but not the true distributions. How should you set your prior to maximize your winning chances?
7. Can you see marginalization in action? What are the parameters used in the estimation? Which do you and which don't you care about.

### Take-homes
- Bayesian probabilities break the asymmetry between parameters and data: everything is a random variable with a distribution
- Inference by Bayes' rule: posterior is proportional to prior x likelihood. The posterior represents how much I expect the parameter's value to be in certain regions of the parameter space
- Conceptually, the prior flatness should reflect the knowledge you have about the parameters. If you look at a flat prior and think "well, the parameter could _never_ be in this range", you are not exploiting all the knowledge you have.
- In practice, the influence of the prior is small, provided you have enough data or your prior is not at once assertive and wrong.
- Conjugate priors measure the strength of the prior as pseudo-counts, i.e., as the equivalent "convincing power" of that number of data samples.
- The likelihood is the part where usually everyone agrees (including classical inference), because it's relatively easier to establish how data comes from the parameters than vice-versa.
- Bayesian inference is often called subjective, but that is true only in the sense that different actors may have different knowledge about the problem, and be exposed to different data. Two Bayesian actors starting from the same "state of ignorance" and seeing exactly the same data must reach _exactly_ the same conclusions, no matter how they organize the inference computation.

### Dialog

 "What's the chance that these samples come from different means?" you ask.

"Essentially 100%," says Bayesian inference.

"Really? Why's that?"

"Well, the event of the means being _exactly_ the same has probability 0. So the means are different almost surely."

"I see what you mean, even if that answer doesn't really help me. Wait. Doesn't the question bother you?"

"Why would it? It's a perfectly valid question."

"Someone told me that either the means are true or they aren't, there's no chance involved."

"I disagree. Chance is not about what things _are_: it's about what _we know_ about them."

"But what do we know about the means, then?"

"That is _exactly_ encoded by this posterior distribution."

"But how do I know that the means come from that distribution?"

"That question makes no sense. As I said, the distribution is not about the means, it's about your knowledge about them."

(...)  [TODO: revise and complete]

### Links

https://en.wikipedia.org/wiki/Bayes%27_theorem
https://en.wikipedia.org/wiki/Bayesian_inference
https://en.wikipedia.org/wiki/Cox%27s_theorem
https://en.wikipedia.org/wiki/Cromwell%27s_rule


## Inferences

### Inference — Bayesian joint posterior over (Δ, σ²_d) via a Normal-inverse-gamma conjugate prior

- **Short description:** jointly estimating the effect Δ and the noise σ² from both groups (between-subjects, pooled), with a conjugate Normal-inverse-gamma prior, then reporting the marginal posterior over Δ (the variance integrated out as a nuisance).
- **Mathematical model:**
  - *Prior (Normal-inverse-gamma):* `Δ | σ²_d ~ Normal(m₀, σ²_d/κ₀)` and `σ²_d ~ Inv-Gamma(a₀, b₀)`. Here `σ²_d = 2σ²` (the difference-of-means variance scale). m₀, κ₀, σ²₀ are sliders. κ₀ is the prior weight on **Δ only** (its pseudo-observations). The variance prior is **decoupled** and held fixed: a₀ = ν₀/2 with ν₀ a gentle, fixed pseudo-count (ν₀ = 4 ⇒ a₀ = 2), kept above 2 so E[σ²] is finite. b₀ then pins the **expected** variance to σ²₀: `b₀ = (a₀−1)·2σ²₀`, so E[σ²] = b₀/(a₀−1)/2 = σ²₀ — no a₀ = 1 singularity. (Earlier drafts coupled `a₀ = κ₀/2` and pinned the *mode*; that made the mean diverge for κ₀ ≤ 2 and produced a discontinuity in the prior at κ₀ = 2.)
  - *Likelihood (between-subjects, pooled):* the difference of group means `Δ̂ = x̄_B − x̄_A ~ Normal(Δ, σ²_d/n)` together with the pooled within-group sum of squares (`2n − 2` df, the grand mean marginalised out under a flat prior). The two groups are **independent — there is no pairing** — so σ²_d is estimated from all 2n observations, not from n formed differences. The integrated likelihood over (Δ, σ²_d) carries variance-power (2n−1)/2.
  - *Posterior:* Normal-inverse-gamma with updated `(mₙ, κₙ, aₙ, bₙ)`, where `aₙ = a₀ + n − ½` (the pooled 2n−2 df plus the ½ borrowed by completing the square on the mean). The *marginal over Δ* is a located-scaled Student-t with `df = 2aₙ`, `loc = mₙ`, `scale = √(bₙ/(aₙ·κₙ))`. In the flat-prior limit (a₀=−½, b₀=0, κ₀=0) it reduces **exactly** to step 1/5's pooled t: `df = 2n−2`, `loc = Δ̂`, `scale = se`.
- **How the posterior was computed:** *conjugate prior–likelihood → closed-form (analytic) update* (`nigUpdate`). The 2-D (Δ, σ²) surfaces are evaluated on a grid only for the equal-volume heatmap drawing (in the variance's own scale, so no Jacobian); the inference itself is the analytic conjugate update, and the Δ marginal is the analytic Student-t.
- **Critical code** (`03-bayesian-approach.html`, `render`):

```js
const { a0: A0, b0: B0 } = variancePrior(s0S.value);                  // shared prior (steps 3 & 4)
const { dbar, ss2 } = drawSample(pool, { delta, sigma, n });          // Δ̂ and the pooled SS (df 2n−2)
const { kn, mn, an, bn } = nigUpdate(m0, k0, A0, B0, n, dbar, ss2);   // conjugate update
...
const pri = nigMarginalDelta(m0, k0, A0, B0);
const pos = nigMarginalDelta(mn, kn, an, bn);                          // Student-t marginal over Δ
```

  with the closed-form update and marginal in `assets/stats.js`:

```js
function nigUpdate(m0, k0, a0, b0, n, dbar, ss2) {   // ss2 = pooled SS on the σ²_d scale (df 2n−2)
  const kn = k0 + n;
  const mn = (k0 * m0 + n * dbar) / kn;
  const an = a0 + n - 0.5;                            // pooled 2n−2 df + ½ from the mean term
  const bn = b0 + ss2 / 2 + k0 * n * (dbar - m0) ** 2 / (2 * kn);
  return { kn, mn, an, bn };
}
function nigMarginalDelta(m, k, a, b) {
  return { loc: m, scale: Math.sqrt(b / (a * k)), df: 2 * a };   // df = 2a
}
```

# Step 4: Making a decision

## Main argument

The Bayesian treatment is the logical way to think under uncertainty. In other words, given what you know, the probabilities it outputs are the unique fully-coherent conclusion. Even with the cost model known by both, a decision based on classical inference has a free parameter (alpha) that must be chosen externally to the decision apparatus, while the Bayesian decision apparatus can marginalize (integrate) over all nuisances to make a decision that maximizes expected reward. The consequence is that under a system of bets, the Bayesian decision procedure will inevitably make better bets (in expectation), provided that you don't start from too wrong _and_ too strong assumptions.

### Try this

1. Play with the sliders and watch which panels respond to what: the *truth & sample* move all four panels, the *prior* moves the marginal (2) and the Bayesian decision (3), and the *reward* moves both decision panels (3, 4). Notice that α moves only the confidence-interval panel (4) — the Bayesian rule never looks at it.
2. Reset to defaults and read the two decision panels side by side. Panel 3 weighs the *whole* posterior against the reward — the gold lens is the expected reward E[U], and the rule bets iff E[U] > 0. Panel 4 keeps only the interval's two edges and acts on the *worst* one.
3. Double-click panel 3. The reward model sweeps across and the integrand fills in behind it; the number it accumulates is E[U]. That single integral — posterior × reward, summed — is the entire Bayesian decision.
4. Run a simulation. Compare "earned" for the two rules over 40 samples: the Bayesian rule earns more on average. Check the other tallies too — "earned if always bet" shows what folding the bad bets buys you, and "earning ceiling" is the unreachable best (betting exactly on the truth every time).
5. Raise the reward penalty c (or the degree p): mistakes hurt more and both rules turn cautious — but the CI rule, judging by its worst edge, folds first and leaves reward on the table the Bayesian rule still collects.
6. Drag α down toward 0.005. The CI widens, its worst edge slides further out, and the CI rule folds almost everywhere. The Bayesian rule, blind to α, keeps betting wherever E[U] > 0. Ask yourself: why should a *false-positive rate* gate a decision that's really about *reward*?
7. Push n up to 100. Both rules sharpen and converge — a tight posterior leaves the CI's worst case no longer pessimistic. The gap between the rules is a small-sample phenomenon.
8. Now set the prior at once *wrong and strong*: m₀ far from the true Δ and κ₀ near 100. Simulate. The Bayesian edge can vanish or reverse — the guarantee holds only if you don't start from assumptions that are both confident and mistaken.

### Take-homes
- A posterior is an input to a decision, not the decision itself; pair it with a reward model and act on expected utility — the whole of Bayesian decision theory.
- A confidence interval allows bounding the likely rewards but using only the best and worst case to guide the decision underperforms in comparison with Bayesian decision, which exploits all information available.


### Dialog

TODO!


### Links
https://en.wikipedia.org/wiki/Dutch_book_arguments
https://en.wikipedia.org/wiki/Decision_theory
https://en.wikipedia.org/wiki/Expected_utility_hypothesis

## Inferences

### Inference — Bayesian decision by expected utility under the posterior (vs. the classical CI worst-case rule)

- **Short description:** turning the step-4 posterior over Δ into an *action* (bet / fold) by maximizing expected utility; the inferential object is the expected reward `E[U] = ∫ reward(Δ − Δ̂)·P(Δ) dΔ` averaged against the posterior, contrasted with a classical (1−α) confidence-interval rule read worst-case.
- **Mathematical model:**
  - *Prior / Likelihood / Posterior:* identical to step 3 — the same between-subjects pooled Normal-inverse-gamma conjugate posterior (κ₀ weighting the mean only, the variance prior decoupled with fixed a₀ = ν₀/2 and `b₀ = (a₀−1)·2σ²₀` pinning the expected variance), marginal over Δ a Student-t. Both rules now read the same pooled variance estimate: the Bayesian posterior and the classical CI (`tTest`, df = 2n−2) coincide in the flat-prior limit, so the head-to-head is a clean prior-vs-α contrast rather than a model mismatch.
  - *Reward (utility):* `U(x) = 1 − |Δ − x|^p − c·|Δ − x|`, clamped at −1; sliders for degree p and linear penalty c.
  - *Bayesian rule:* bet iff `E[U] > 0`, betting on the posterior mode mₙ.
  - *Classical rule:* take the (1−α) CI `d̄ ± t·se`; with no distribution to average over, act worst-case — bet iff the reward at the interval edge stays ≥ 0.
- **How computed:** the posterior is the same analytic conjugate update (`nigUpdate` / `nigMarginalDelta`); the *decision* expectation is then obtained by **numerical integration** (trapezoid rule, M = 4000 over [−12, 12]) of posterior × reward — not a closed form, since the reward is non-conjugate.
- **Critical code** (`04-the-decision.html`, `computeState` + helpers):

```js
const { kn, mn, an, bn } = nigUpdate(m0, k0, A0, B0, n, dbar, ss2);
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

# Step 5: A Bayesian t-test?

## Main argument

It's tempting to ask for a "Bayesian t-test": Step 1's familiar question, answered with Step 3's machinery. We can certainly build one — put a flat (improper) prior on the effect Δ and on the noise, multiply by the Gaussian likelihood, and read off a posterior over Δ. But that posterior is not new. When the prior is flat, prior × likelihood *is* the likelihood, and the result is exactly Step 1's t-distribution — merely re-centred on the observed difference Δ̂ and re-read as a degree of belief rather than a sampling distribution. The 95% credible interval is, number for number, the classical confidence interval Δ̂ ± t·se.

So what did we gain? Nothing — and that is precisely the lesson. The power of Steps 3 and 4 came from two moves this page deliberately undoes. First, an *informative* prior: real knowledge, priced in pseudo-observations, pulling the posterior toward what we already know — a flat prior throws that away, declaring every effect from 0.1 to a million equally believable. Second, a *decision* that integrates the whole posterior against a reward — here we instead collapse the posterior to a single yes/no (does the interval exclude zero?), which is just the significance test wearing Bayesian clothes, hard threshold α and all. Maximum a posteriori under a flat prior is maximum likelihood; binarizing the posterior is the classical test. Use Bayesian inference to imitate a classical procedure and you should expect classical results: throw the information away and you don't get to keep its benefits. (The credible interval still earns its keep as an honest visual summary — useful when an audience would balk at a point estimate with no sense of its spread.)

(This is a bit of a tangent: the idea is to show that we cheat by using Bayesian inference while throwing away information, we should not expect to reap the benefits. This is a bit of a side-lesson: Steps 1-4 and the exercise of Step 6 should not depend on it.)

### Concepts

TODO!

### Try this

1. Play with the sliders. Panel 2 is read as a *posterior* over Δ — but it's the very curve from Step 1's t-test, only re-centred on the observed Δ̂ and relabelled "degrees of belief" instead of "sampling distribution."
2. Reset to defaults. Hover panel 2 to read off the plausibility of any effect size; the shaded band is the 95% credible interval (the HDI).
3. Open Step 1 in another tab with the same Δ, σ², N, α. The interval Δ̂ ± t·se is *identical* — with flat priors the Bayesian credible interval and the classical confidence interval coincide number for number.
4. Watch the verdict pill. "Does the HDI exclude 0?" is exactly Step 1's "is it significant?" in Bayesian clothing — we've quietly smuggled the hard threshold back in.
5. Drag α from 0.05 toward 0.005: the band widens and "tolerable risk" can flip to "intolerable." A genuinely *credible* interval shouldn't depend on a *false-positive rate* — α's presence here is the tell that we're still running the classical test.
6. Set Δ small (≈0.2–0.3) with small N: the HDI straddles 0 — "intolerable." Raise N until it clears 0. It's Step 1's evidence-vs-magnitude tug of war all over again, just renamed.
7. Run a simulation. The 40 repeated HDIs behave exactly like Step 1's confidence intervals — about 1 in 20 misses the true Δ. The picture is identical because the arithmetic is identical.
8. Step back: we used Bayes' theorem (flat prior × likelihood) and then collapsed the whole posterior to one yes/no, ignoring any reward. We did all the Bayesian work and kept none of the Step-4 payoff.

### Take-homes
- Priors don't have to be informative and can even be improper, but almost always using such a flat prior is a missed opportunity
- You can dress-up hard-threshold classical tests in Bayesian concepts, but you shouldn't expect to reap the benefits of Bayesian inference
- This is not the Bayesian treatment: (1) I'm using maximum a posteriori with an uninformative prior = same as maximum likelihood. (2) Credible intervals force a binarization of the decision that goes against the Bayesian ethos of using all information available
- Still, credible intervals (the Bayesian answer to classical confidence intervals) are a useful visual summary for communicating results, especially when the public would be disturbed by the lack of interval estimation

### Dialog

TODO!

### Links

TODO!
https://en.wikipedia.org/wiki/Maximum_a_posteriori_estimation
https://en.wikipedia.org/wiki/Credible_interval


## Inferences

### Inference — Bayesian posterior over the mean difference Δ, with flat priors

- **Short description:** estimating the effect size Δ = μ_B − μ_A as a full posterior distribution under flat (non-informative) priors — which yields a shifted, scaled Student-t, i.e. step 1's t-distribution re-centred on the observed difference.
- **Mathematical model:**
  - *Likelihood:* Gaussian between-subjects, as in step 1 (both groups, common σ²).
  - *Prior:* flat (improper) on Δ and on log σ — "let the data speak."
  - *Posterior over Δ* (σ² marginalized out): a located-scaled Student-t centred at `Δ̂ = x̄_B − x̄_A`, scale = the t-test's standard error `se`, `df = n + m − 2`.
- **How the posterior was computed:** *analytic* — the flat-prior result coincides with the classical sampling object, so the page re-uses `tTest`'s `(dObs, se, df)`; no numerical integration. The (1−α) HDI is exactly the equal-tailed band `Δ̂ ± t·se`. This is also **exactly the flat-prior limit of the step 3–4 conjugate model** (a₀=−½, b₀=0, κ₀=0 ⇒ df = 2n−2, scale = se), so steps 3–5 share one between-subjects pooled variance model — flattening the step-3 prior lands precisely here.
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

- Step 4 carries the thesis of the whole tutorial (Making a decision apparatus
  marginalizes nuisances and maximizes expected reward; classical inference leaves α as a
  free parameter chosen outside the apparatus). It must stay live and must land on the
  continuous model the learner has lived in for four steps — one new idea (decision) on
  established ground, not two new ideas at once. So we do **not** cut Step 5 to make room.
- The continuous thread (1→3→4) is a complete spine. Step 2 is the only binary step and
  the only one that never gets a Bayesian answer; its closing line even promises "Bayes'
  theorem — where we go next" and then the thread is dropped. Step 6 is that owed answer.
- For this audience the strongest capstone is an **inversion**: for five steps they dragged
  sliders on a machine someone else built; here they build the machine behind the slider.
  "You've been playing with the posterior — now compute one."
- The deliberate incompleteness *is* the assignment: the scaffold and the questions are
###   fully built, but the morals and part of the concept table
  fills in. Transfer (apply the pipeline to new data) is exactly the move that consolidates
  it.

## Main argument

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

### concepts
- The prior **is** Step 2's base rate, written as `Beta(α, β)` — a few imaginary cases,
  priced in pseudo-observations (the pseudo-count α + β echoes Step 3's κ₀).
- Normalization is "the main challenge of Bayesian inference" (Step 3) shrunk to one
  Riemann sum: divide the curve by its own area.
- The decision is Step 4's rule on new data: pair the posterior with a reward, act on
  expected utility.

### Take-homes
- The same four boxes (prior, likelihood, posterior, decision) are, on binary data, ten
  lines of array arithmetic.
- Integrating over the things you want to ignore to get the ones you care about is the bread and butter of Bayesian inference. Marginalization and expected values are different facets of this strategy.
- Integration becomes a problem as the number of parameters grows => curse of dimensionality. Things like MCMC and, in particular, HMC have greatly broadened the universe of problems for which we can apply Bayesian inference, but this is still a problematic or even completely intractable step for many problems.

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

# Further reading

TODO!

- Doing Bayesian Data Analysis, Kruschke => easier, accessible, includes a practical cookbook

- Statistical rethinking, McElreath  => rigorous but accessible, great explanation for conceptual issues, companion course available on the web

- Bayesian Data Analysis, Gelman et al. => _the_ rigorous textbook for modern Bayesian inference, very comprehensive, steep in difficulty
