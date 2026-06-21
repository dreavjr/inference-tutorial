"use strict";
/* stats.js — the statistical model, built on math.js and shared by every page:
   the Student-t family (native, no library) and normal CDF; sampling pools and
   drawSample (turn a pool into one experiment); the two-sample t-test and its
   power; and the Normal-inverse-gamma conjugate model with helpers that read it
   back as 1-D curves over Δ. Pages compose these; they hold no inference of their
   own — only how to draw it. */

const POOL = 100; // fixed pool of random draws — dots don't jitter while you drag

/* ---- sampling ---------------------------------------------------------- */

function makePool(seed) {
  const r = mulberry32(seed === undefined ? (Math.random() * 4294967296) >>> 0 : seed);
  const norm = () => {
    let u = 0;
    while (u === 0) u = r();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * r());
  };
  const fill = (f) => Array.from({ length: POOL }, f);
  return { z1: fill(norm), z2: fill(norm), j1: fill(r), j2: fill(r) };
}

/* Turn a pool into one experiment: two groups of n drawn from independent
   populations centred at ∓Δ/2 with a shared spread σ. The two groups are NOT
   paired (there are no subjects to pair), so we summarise them the between-subjects
   way: the difference of group means Δ̂ = x̄_B − x̄_A, and the pooled within-group
   sum of squares ss2 — expressed on the difference-variance scale σ²_d = 2σ² so the
   conjugate model can work in Δ. ss2 = 2·[Σ(x_{A,i}−x̄_A)² + Σ(x_{B,i}−x̄_B)²] carries
   2n−2 degrees of freedom and E[ss2] = (2n−2)·σ²_d. Δ̂ has variance σ²_d/n = 2σ²/n.
   Every page draws its sample this way, so the same pool + knobs always yield the
   same data across the series. */
function drawSample(pool, { delta, sigma, n }) {
  const xa = pool.z1.slice(0, n).map(z => -delta / 2 + sigma * z);
  const xb = pool.z2.slice(0, n).map(z =>  delta / 2 + sigma * z);
  const dbar = mean(xb) - mean(xa);                          // Δ̂ = x̄_B − x̄_A
  const ss2  = 2 * (n - 1) * (svar(xa) + svar(xb));          // pooled SS on the σ²_d scale, df = 2n−2
  return { xa, xb, dbar, ss2 };
}

/* ---- Student's t distribution ------------------------------------------ */

function tPdf(x, df) {
  return Math.exp(logGamma((df + 1) / 2) - logGamma(df / 2))
         / Math.sqrt(df * Math.PI) * Math.pow(1 + x * x / df, -(df + 1) / 2);
}

function tCdf(x, df) {
  const p = 0.5 * ibeta(df / 2, 0.5, df / (df + x * x));
  return x >= 0 ? 1 - p : p;
}

function tInv(p, df) { // quantile, by bisection
  if (p < 0.5) return -tInv(1 - p, df);
  let lo = 0, hi = 150;
  for (let i = 0; i < 100; i++) {
    const mid = (lo + hi) / 2;
    if (tCdf(mid, df) < p) lo = mid; else hi = mid;
  }
  return (lo + hi) / 2;
}

const tCrit = (alpha, df) => tInv(1 - alpha / 2, df);   // two-sided critical t

/* ---- normal distribution ----------------------------------------------- */

const normPdf = (x) => Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);
const normCdf = (x) => 0.5 * (1 + erf(x / Math.SQRT2));

/* ---- two-sample t-test, pooled variance (= scipy.stats.ttest_ind) ------- */

function tTest(xa, xb) {
  const n = xa.length, m = xb.length;
  const sp2 = ((n - 1) * svar(xa) + (m - 1) * svar(xb)) / (n + m - 2);
  const dObs = mean(xb) - mean(xa);                  // observed difference of means
  const se = Math.sqrt(sp2 * (1 / n + 1 / m));       // its standard error
  const df = n + m - 2, t = dObs / se;
  return { t, df, p: 2 * (1 - tCdf(Math.abs(t), df)), dObs, se, sp: Math.sqrt(sp2) };
}

/* Power of the test at the true Δ. The statistic is a noncentral t(df, ncp):
   t = Z/√(W/df) with Z ~ Normal(ncp, 1), ncp = Δ/se_true, and W ~ χ²_df the
   independent scale from the pooled-variance estimate. So power is the normal-tail
   formula Φ(ncp − t_c·s) + Φ(−ncp − t_c·s) — with s = √(W/df) — averaged over the
   χ² spread of the estimated SD, not pretending it is known (the error a plain
   normal makes, badly so at small n). Simpson over the χ² bulk; ~1e-4 accurate.
   β = 1 − power is the chance this n misses a real effect. */
function tTestPower({ delta, sigma, n, alpha, df }) {
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
}

/* ---- Normal-inverse-gamma conjugate model ------------------------------
   joint log density over (Δ, σ²_d), the integrated Gaussian log-likelihood of the
   two-sample experiment (grand mean marginalised out under a flat prior), and the
   conjugate prior→posterior update used by steps 3 & 4. The likelihood is the
   between-subjects pooled model: Δ̂ ~ Normal(Δ, σ²_d/n) and the pooled within-group
   SS carries 2n−2 df, so the integrated likelihood has v-power (2n−1)/2. */

function logNIG(d, s2, m, k, a, b) {
  return 0.5 * Math.log(k / (2 * Math.PI * s2)) - k * (d - m) ** 2 / (2 * s2)
       + a * Math.log(b) - logGamma(a) - (a + 1) * Math.log(s2) - b / s2;
}

function logLik(d, s2, n, dbar, ss2) {
  return -((2 * n - 1) / 2) * Math.log(2 * Math.PI * s2) - (ss2 + n * (dbar - d) ** 2) / (2 * s2);
}

/* Conjugate update for the pooled two-sample model. The mean Δ carries weight n
   (Δ̂ has variance σ²_d/n), while the variance picks up the pooled 2n−2 df: that and
   the +½ borrowed by completing the square on the mean give aₙ = a₀ + n − ½. In the
   flat limit (a₀=−½, b₀=0, κ₀=0) the Δ marginal is t with df = 2n−2 and
   scale² = (ss2/2)/((n−1)n) = se² — exactly step 1/5's pooled t-test. */
function nigUpdate(m0, k0, a0, b0, n, dbar, ss2) {
  const kn = k0 + n;
  const mn = (k0 * m0 + n * dbar) / kn;
  const an = a0 + n - 0.5;
  const bn = b0 + ss2 / 2 + k0 * n * (dbar - m0) ** 2 / (2 * kn);
  return { kn, mn, an, bn };
}

/* ---- reading the model as 1-D curves over Δ ----------------------------- */

/* a located, scaled Student-t density: x ↦ t((x−loc)/scale; df)/scale. The flat-prior
   posterior of step 3 and both NIG marginals of steps 4-5 are exactly this shape. */
const locScaleTpdf = (loc, scale, df) => (x) => tPdf((x - loc) / scale, df) / scale;

/* the marginal of a Normal-inverse-gamma over Δ (σ²_d integrated out) is a
   Student-t: df = 2a, centre m, scale √(b/(a·κ)). */
function nigMarginalDelta(m, k, a, b) {
  return { loc: m, scale: Math.sqrt(b / (a * k)), df: 2 * a };
}

/* the Gaussian profile likelihood over Δ, peak-normalised to 1 — a relative
   support weighting, not a density: centre Δ̂, sd = the pooled standard error
   se = √((ss2/2)/((n−1)·n)) = sp·√(2/n). */
function likelihoodProfile(dbar, ss2, n) {
  const sd = Math.sqrt(ss2 / 2 / (n - 1) / n);
  return (x) => Math.exp(-0.5 * ((x - dbar) / sd) ** 2);
}

/* ---- where the (Δ, σ²) surface peaks, as the group variance σ² (step 3) -------
   The heatmap is drawn in (Δ, σ²) with NO Jacobian, so the dot marks the true joint
   NIG mode: at Δ = m the slice is p(σ²_d) ∝ (σ²_d)^−(a+3/2)·exp(−b/σ²_d), peaking at
   σ²_d = b/(a+3/2); in the group variance that is σ² = b/(a+3/2)/2. */
const nigVarMode = (a, b) => b / (a + 1.5) / 2;

/* the mean of the surface in the group variance σ² — E[σ²_d] = b/(a−1), so
   E[σ²] = b/(a−1)/2. Exists only for a > 1; callers hide the ✕ otherwise. */
const nigVarMean = (a, b) => b / (a - 1) / 2;

/* where the integrated-likelihood surface peaks, as the group variance σ². At Δ=Δ̂
   it maximises at σ²_d = ss2/(2n−1), so in the group variance σ² = ss2/(2n−1)/2. */
const likVarMode = (ss2, n) => ss2 / (2 * n - 1) / 2;
