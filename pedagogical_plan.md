# Step 1: Classical statistical test

## Pedagogical plan

Classical inference doesn't accept parameters as random variables because that would clash with the frequentist notion of probability. To bypass that blocking point, we define the "experiment" as the source of variation. Each experiment is associated to a random sample, thus the importance of sampling distributions and test statistcs. Instead of conceptualizing the uncertainty of a parameter, we will instead conceptualize significance: how often our experiment leads to the right decision (test) or to the right range of possibilities (confidence interval).

We cannot control directly for power, but, under the assumptions of the statistical model, we can compute how much power our current sample size gives us, or vice versa.

Concepts:
- Experimental sample, sampling distribution, test statistic
- Tail areas and p-values
- Statistical tests and confidence intervals

Take-homes:
- Classical inference is easy to apply (at least for traditional problems such as comparing population means) but hard to interpret. Classical tests have lots of assumptions. Significance and power apply across multiple experiments, not one particular experiment
- The bread and butter of classical inference for parameters is maximum-likelihood to locate a parameter and use sampling probabilites and tail-area analysis for determining uncertainty

## Current material

Page: `01-the-t-test.html` — "How classical tests work: the t-test"

Intro text (summary): Sets up two unobserved populations (treatment A vs. B), each a full distribution with hidden knobs — the true difference of means Δ and the true variance σ². What we observe is a sample of n points per population. Walks through the classical recipe as a 5-step list: (1) state H₀ (Δ = 0) vs. H₁; (2) condense the sample into the test statistic t = (x̄_B − x̄_A) / (s_p·√(2/n)), the observed difference in units of its own estimation noise; (3) under H₀ the difference wobbles as a Student's t scaled by the standard error, drawn on the same axis as the data so D̂ appears in both panels — the band is far tighter than the populations and tightens like 1/√n; (4) the p-value is the two-sided tail area at least as extreme as D̂; (5) reject H₀ if p ≤ α (fixed at 0.05).

Concept table ("a map of the picture" — concept → where to see it below):
- population → the two smooth curves, never observed directly
- sample → the dots above (blue) and below (orange) the line
- H₀ → "blue and orange coincide"; drag Δ to 0 to make it true
- observed difference → the marker D̂ (with its t value) in the right panel
- its null distribution → the dark scaled-t curve, how the difference scatters if H₀ holds
- p-value → the shaded tail area beyond ±D̂, also in numbers
- α → the bar p must clear; ticks mark the critical difference ±t·se

Try this:
1. Drag Δ to 0 (H₀ true), resample a dozen times: D̂ and p bounce; ~1 in 20 dips below 0.05 — that is α, the false-positive rate.
2. Δ = 0.5, n = 12: real but modest effect, rarely significant (underpowered); raise n to 100 and it's significant almost always — the band tightens like 1/√n while D̂ stays put. Significance measures evidence, not importance.
3. Crank σ² up: null band widens, significance melts away — a tug of war between signal (D̂) and noise (band width se = s_p√(2/n)).
4. Drop n to 3–4: band turns squat and wide, ±t·se critical ticks slide far out; raise n and the band grows so tall and tight it fades out at the top.

The moral:
- The p-value answers one question: "If there were truly no difference, how surprising would my t be?" — a tail probability under H₀, nothing more.
- It is not the probability that H₀ is true, and a rejection is not a measure of effect size or importance.
- "Reject H₀" just stamps positive on the study; how much belief that deserves is the subject of the next step.

# Step 2: What we want vs. can control

## Pedagogical plan

The ideas of power and significance appear throughout the literature with different names: sensitivity and specificity, true positive rate and true negative rate, etc. They are formal properties of the test but do not correspond to what we will observe in the field, because they disconsider base rates. PPV and NPV correspond much more to what we will see in the field ("the test tagged this as positive/negative: how often can I expect it to be right?") while power and significance answer "the sample is positive/negative: how often can I expect the test to tag it right" (which is not as useful, because the true labels are, of course, unknown: if I had them I wouldn't need a test!).

Take-homes:
- Test-bound knobs are easy (change the test)
- Reality-bound knobs are hard (changing the reality is much harder)
- Sensitivity and specificity don't bound Type I and II errors, not even asymptotically. However, if I make assumptions about the base rate, then I can find useful bounds.

## Current material

Page: `02-pvalues-and-power.html` — "What p-values and power can do for you — and what they can't"

Intro text (summary): Takes the step-1 machine that stamps significant / not significant and runs it at scale — a thousand studies (or screened patients, or hypotheses). Each case is either condition + (effect real) or condition − (no effect), and the machine returns predicted + ("significant!") or predicted −. Two famous numbers describe the procedure, but both condition on an "if about the truth" — exactly what you don't know. The question you actually face runs the other way: "My result came out positive — how likely is the effect real?" That is the PPV (precision); its mirror is the NPV. Neither α nor power answers it alone.

Concept table (testing language → diagnostic language → the question it answers):
- 1 − α (significance level) → specificity → "If there is nothing, how often do we correctly stay quiet?"
- power (1 − β) → sensitivity, recall → "If there is something, how often do we catch it?"

Try this:
1. Defaults (α = 0.05, power = 0.80, base rate 50%): PPV ≈ 94% — the comfortable textbook world.
2. Base rate down to 10%: α and power unchanged, yet PPV falls to ≈ 64% (one positive in three is false).
3. Down to 2%: PPV ≈ 25% — three "discoveries" in four are wrong, same respectable test.
4. Try to rescue it without touching the base rate: pushing power to 100% barely helps (PPV ≈ 29%); what works is specificity — push to 99.5% (α = 0.005) and false positives melt away.
5. At base rate 2%, crank total cases from 1,000 to 1,000,000: every count grows 1000×, no percentage moves. PPV is a property of proportions — sample size will not save you.

The moral:
- Sensitivity (power) and specificity (1 − α) are properties of the test; they condition on the truth.
- PPV and NPV are properties of your result; they condition on what you observed and depend inescapably on the base rate.
- Same test, different population → a positive means something different; p-value and power are blind to this, and more data won't fix it.
- Getting from (sensitivity, specificity, base rate) to (PPV, NPV) is exactly Bayes' theorem — where we go next.

# Step 3: A difference in style?

## Pedagogical plan

Bayesian probabilities are not tied to frequencies: they represent what you must conclude when you are reasoning under uncertainty. A probability distribution represent how much you can locate the value of a random variable given all that you know. This small change changes everything, because now parameters don't have to be fixed: their distribution represent what you know about them, not frequencies that depend on repetition.

Concepts:
- Bayesian probabilities break the asymmetry between parameters and data: everything is a random variable with a distribution
- Inference by Bayes' rule: posterior is proportional to prior x likelihood. The posterior represents how much I expect the parameter's value to be in certain regions of the parameter space
- Priors don't have to be informative

Myths to dismiss:
- Bayesian probability is all about priors => you can often use non-informative or little-informative priors (although that is not necessarily a good idea for most problems)
- I can use Bayesian inference to mimick hard-threshold classical tests and still reap the benefits => if you don't add informative priors ("let the data speak for itself" zealousness) and force sharp decisions instead of exploit the uncertainty, Bayesian inference can give exactly the same results as classical inference

Take-homes:
- This is not the full Bayesian treatment: (1) I'm using maximum a posteriori with an uninformative prior = same as maximum likelihood. Credible intervals force a binarization of the decision that goes against the Bayesian ethos.

## Current material

## Pedagogical plan

Page: `03-bayesian-t-test.html` — "A Bayesian t-test?: from is there an effect? to how big is it?"

Intro text (summary): Step 1 asked only "if H₀ were true, how surprising is our t?" and answers with a yes/no stamp — never telling us how large the effect is or how sure we are. The Bayesian move needs no new machinery, just a change of reading: take the same t-distribution and re-centre it on the observed difference instead of on 0, letting it describe the whole range of plausible effect sizes (with flat priors this is exactly the posterior for the mean difference — a shifted, scaled t). Both panels share step 1's square scale; the right panel is the distribution over Δ. Because it is a genuine density, a very certain posterior becomes a tall thin spike that fades out at the top of the box. Hovering the Δ curve maps any candidate effect size to the two population means it implies.

Concept table (you see → it means):
- the curve's peak → the observed difference Δ̂ = x̄₂ − x̄₁, our best estimate
- the curve's width → our uncertainty, wide when data are few or noisy
- the curve's height → it's a real density (area = 1), so concentrating ⇒ taller
- the dashed line at 0 → the null; if it sits out in the tail, the data doubt it
- the green dashed Δ → the real effect, unknown in life, shown here because we simulate

Try this:
1. Hover slowly across the Δ curve: the two implied means slide apart and back — the curve shows which separations the data consider plausible.
2. Is 0 under the curve? With defaults, barely (the Bayesian echo of "significant"); drag Δ toward 0 and the curve drifts until it straddles 0 and "no effect" becomes most plausible.
3. Raise n: the curve concentrates (taller, narrower) around the truth — certainty about size, not just a smaller p; far enough and the spike fades out the top. This is power from the inside.
4. Raise σ²: the curve fattens and sinks — noise widens the band of indistinguishable effect sizes.
5. "New sample" a few times: the curve jumps (its centre Δ̂ is the noisy observed effect) but usually keeps true Δ in its bulk; each study hands you one curve and you never see the green line.

The moral:
- A p-value returns a verdict; this returns an estimate with honest error bars — same t-distribution read the other way round.
- The shape is the message: "significant" can hide a curve spanning trivial to enormous effects; "not significant" can hide a curve tight around zero. The stamp throws that away; the distribution keeps it.
- One sleight of hand remains — "with flat priors" was itself a choice, made explicit in the next step.

# Step 4: Bayesian inference's pipeline: prior -> lihelihood -> posterior -> da capo

## Pedagogical plan

Full Bayeasian treatment: representation of uncertainty, marginalization of nuisance variables (the variance, in this case). Conclusion is not a single value or hard interval: the posterior _is_ the conclusion.

Concepts:
- Prior, likelihood, posterior

Take-homes:
- Likelihoods are not normalized, nor the product prior x likelihood. Normalizing is _the_ main challenge of Bayesian inference
- Priors should reflect domain knowledge and be lax enough to allow the model to "change opinions", but making them _too_ lax leads to missed opportunities

## Current material

Page: `04-full-bayesian.html` — "The full Bayesian approach: tight, lax, and honest"

Intro text (summary): Names step 3's sleight of hand: the flat prior, which declared every effect size (Δ = 0.1, 3, a million) equally believable — not humility but a strong claim. Makes the prior explicit via Bayes as multiplication: posterior ∝ prior × likelihood. The unknowns come as a pair — effect Δ and noise σ² — so the prior is a distribution over that plane; the conjugate choice is the Normal-inverse-gamma: Δ | σ²_d ~ Normal(m₀, σ²_d/κ₀) and σ²_d ~ Inv-Gamma(a₀, b₀). Three panels draw prior (purple), likelihood (terracotta), posterior (dark) as filled equal-volume contours over the (Δ, σ) plane (σ as standard deviation, both axes same scale), with a green cross at the true (Δ, σ). Below: the data and the marginal over Δ as the three 1-D shadows. A simplification keeps the algebra on one screen: the model is fed the n per-pair differences dᵢ = x_{B,i} − x_{A,i}, each Normal(Δ, σ²_d) with σ²_d = 2σ².

Concept table (knob → what it says → here):
- m₀ → where you expect Δ to be, before data → a slider
- κ₀ → how strongly, in pseudo-observations: your hunch is worth κ₀ imaginary data points → a slider (try the extremes)
- a₀, b₀ → the same two roles, for the variance σ² → b₀ is a slider (where you expect σ²); a₀ held fixed and gentle

Try this:
1. Defaults (m₀ = 0, κ₀ = 1, n = 12): posterior teardrop sits between prior and likelihood, much closer to the likelihood (12 real obs outvote 1 imaginary); the marginal's dark curve peaks just short of the terracotta — a whisper of skepticism.
2. Go lax: κ₀ → 0.1 — prior flattens to a vast sheet, posterior becomes the likelihood (you've rebuilt step 3); resample at n = 12 and the posterior chases noise across the plane.
3. Skeptical but fair: m₀ = 0, κ₀ ≈ 10 — posterior pulled toward 0 (slightly biased when truth is 0.8) but resampling barely moves it; with small n, trading bias for stability is a bargain — what a prior buys.
4. Tight and wrong: m₀ = −1, κ₀ = 100 — 100 imaginary obs outvote 12 real, posterior parks near −1, green cross far outside the teardrop; raise n to 100 and the data claw back. Tight priors are loud.
5. Crank n to 100 and sweep κ₀ from 0.1 to 10: posterior barely budges — enough data washes out any reasonable prior. Priors matter when data are scarce, exactly when "non-informative" is most dangerous.

The moral:
- A flat prior is not absence of opinion; it is the opinion that preposterous effect sizes are as credible as modest ones — and with little data that leaks into the posterior as width and wobble.
- κ₀ prices your conviction in the only honest currency: data. Ask "how many observations is my hunch worth?" — both 0.1 (abandons to noise) and 100 (deaf to evidence) can be wrong.
- When n is small, weakly-informative beats non-informative; when n is large the argument dissolves on its own — so the cases where a lax prior "plays it safe" are precisely the cases where it doesn't.
- The prior is a dial to disclose, not a sin to hide: state it, vary it, show the posterior barely moves — that is what a robust claim looks like.

# Step 5: Everyone is entitled to their opinion (but they'll lose money if they disagree with me)

## Pedagogical plan

The full Bayesian treatment is the logical way to think under uncertainty. In other words, given what you know, the probabilities it outputs are the  unique fully-coherent conclusion. Even with the cost model known by both, a decision based on classical inference has a free parameter (alpha) that must be chosen externally to the decision apparatus, while the Bayesian decision apparatus can marginalize over all nuisances to make a decision that maximizes expected reward. The consequence is that under a system of bets, the Bayesian decision procedure will inevitably make better bets (in expectation), provided that you don't start from too wrong _and_ too strong assumptions.

## Current material

Page: `05-the-decision.html` — "Making a decision"

Intro text (summary): Step 4 handed us a posterior — a whole distribution over Δ — but a distribution is not yet a decision; sooner or later you must act (ship the drug, set the dose, choose B over A). Bayesian decision theory closes the gap with one rule: act to maximise expected utility under your posterior. It needs the posterior plus a reward model U(x) = 1 − |Δ − x|^p − c·|Δ − x|: land on the truth and earn a full 1; miss and two penalties bite, a polynomial |Δ − x|^p and a linear c·|Δ − x|. For positive c the reward crosses zero at a break-even error ε (the largest profitable miss); raising p oddly pushes ε outward for small misses; with c = 0, ε = 1. The bottom row poses the sharp question — given my belief, should I bet at all, and on what value? — both readings betting on the peak Δ̂ but deciding whether to bet differently. Left: the honest posterior, where the Bayesian rule averages reward against the whole posterior, E[U] = ∫ reward(Δ − Δ̂)·P(Δ) dΔ, and bets only when that is positive (the gold lens is the overlap). Right: the classical reading, a (1 − α) confidence interval drawn as a flat error bar (MLE Δ̂ centred, margin E either way); with no distribution to average over the rule turns worst-case — bet iff E ≤ ε, so even the interval's worst boundary still pays. Each panel shows an estimated pill (the rule's forecast — posterior's expected reward, or CI's min…max) turning green to bet / grey to fold, and an actual pill (the payoff once hidden Δ is revealed, or zero for a refused bet) green for a win / red for a loss / grey for a no-play. The gap between the pills is the price of not knowing the truth.

Concept table (knob → what it does → here):
- p → the polynomial degree, how the penalty curves: a gentle bowl at p = 1, a flat top with steep walls as p climbs → a slider (1…10)
- c → the linear penalty, a flat charge per unit of error; raising it tilts the whole reward down → a slider (−10…10)
- Δ → where the reward is centred (the truth itself), drawn only because we simulate → the green line

Try this:
1. Defaults: with c = 0 the only penalty is polynomial, so break-even ε = 1 and the CI fits inside it; both rules bet and win — the paradigms agree.
2. Drop n to 3–4: the CI balloons past ε = 1, its worst-case boundary turns negative, the classical rule refuses (banks 0); the compact posterior still expects positive return — it bets and wins. The Bayesian edge is sharpest when data are scarce.
3. Crank n to 100: both beliefs sharpen, both bet, actual payoffs converge near peak reward 1 — with enough data the rules agree again.
4. Raise c (linear penalty): break-even ε pulls in below 1; first the wide-interval classical rule folds, push further and even the posterior's expectation drops below zero so the Bayesian folds too. (With c > 0, a higher p widens ε, nudging both back toward betting.)
5. Tighten the prior (κ₀ ↑) and shift m₀ off the truth: the posterior's bet pulls away from the interval's midpoint — watch a confident-but-wrong prior keep its estimated pill green and betting while the actual pill slides red: overconfidence made quantitative.

The moral:
- A posterior is an input to a decision, not the decision itself; pair it with a reward model and act on expected utility — the whole of Bayesian decision theory.
- A confidence interval is not a probability distribution; reading it as "uniform over plausible values" is a specific, usually worse, belief — and a decision made on it earns less.
- The shape the posterior carries (mass near the estimate, thinning into the tails) is what a uniform block discards and what a reward centred on the truth pays you to keep.
- Same data, same interval width — yet a different decision, because of how honestly the uncertainty inside the interval is shaped.
- Two numbers, never one: what a method thinks it will earn and what it actually earns. A trustworthy rule keeps the pair close; a tight wrong prior, or a mis-shaped interval, pries them apart.

=> Bibliography
- Accessible resources
- Practical materials and tools
- Academic material

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
  fully built, but the morals and part of the concept table are left as prompts the learner
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

Take-homes:
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

Page: `do-it-yourself.html` — "Do it yourself: Bayes on binary data"

A notebook mock-up (page-local CSS for the cell vocabulary, to be promoted to shared CSS
once the design settles and the engine is wired). Cells alternate **provided** (read-only
scaffolding / answer keys) and **your code** (signature + docstring contract + `# TODO`);
the **Run** buttons are intentionally inert. A callout at the top states plainly that the
page is deliberately unfinished and the engine is not connected yet.

Intro text (summary): Frames Step 6 as the Bayesian answer Step 2 never got, and as the
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
- Concept table mapping `prior_kernel`→Step 4 (κ₀), `normalize`→Step 4, `expected_utility_act`/
  `decide`→Step 5, with a deliberate `TODO` row.

Not yet done: wire the in-browser Python engine and the plot helpers (`plot_three`,
`plot_overlay`, the E[U]/bet-fold output); promote the notebook CSS to shared; decide whether
to link the page from `index.html` (currently unlinked so it doesn't read as a sixth live
step); resolve the reward-model open question above.
