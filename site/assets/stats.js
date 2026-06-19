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

/* Turn a pool into one experiment: two groups of n drawn from populations
   centred at ∓Δ/2 with spread σ, plus the n per-pair differences
   dᵢ = x_{B,i} − x_{A,i} ~ Normal(Δ, 2σ²). Every page draws its sample this way,
   so the same pool + knobs always yield the same data across the series. */
function drawSample(pool, { delta, sigma, n }) {
  const xa = pool.z1.slice(0, n).map(z => -delta / 2 + sigma * z);
  const xb = pool.z2.slice(0, n).map(z =>  delta / 2 + sigma * z);
  const d  = xb.map((v, i) => v - xa[i]);
  const dbar = mean(d);
  const ssd  = d.reduce((s, v) => s + (v - dbar) ** 2, 0);   // Σ(dᵢ − d̄)²
  return { xa, xb, d, dbar, ssd };
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
   joint log density over (Δ, σ²), the Gaussian log-likelihood of n observed
   differences, and the conjugate prior→posterior update used by steps 4 & 5. */

function logNIG(d, s2, m, k, a, b) {
  return 0.5 * Math.log(k / (2 * Math.PI * s2)) - k * (d - m) ** 2 / (2 * s2)
       + a * Math.log(b) - logGamma(a) - (a + 1) * Math.log(s2) - b / s2;
}

function logLik(d, s2, n, dbar, ssd) {
  return -(n / 2) * Math.log(2 * Math.PI * s2) - (ssd + n * (dbar - d) ** 2) / (2 * s2);
}

function nigUpdate(m0, k0, a0, b0, n, dbar, ssd) {
  const kn = k0 + n;
  const mn = (k0 * m0 + n * dbar) / kn;
  const an = a0 + n / 2;
  const bn = b0 + ssd / 2 + k0 * n * (dbar - m0) ** 2 / (2 * kn);
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
   support weighting, not a density: centre d̄, sd √(ssd/((n−1)·n)). */
function likelihoodProfile(dbar, ssd, n) {
  const sd = Math.sqrt(ssd / (n - 1) / n);
  return (x) => Math.exp(-0.5 * ((x - dbar) / sd) ** 2);
}

/* ---- where the (Δ, σ) surfaces peak, as the group std dev σ (steps 4-5) --
   At Δ = m the σ²_d slice of the NIG is Inv-Gamma(a+½, b) — the extra ½ comes from
   the Normal factor's σ⁻¹ — with mode b/(a+3/2); then σ = √(mode/2). */
const nigSigmaMode = (a, b) => Math.sqrt(b / (a + 1.5) / 2);

/* the MLE of the group std dev from the differences: σ̂²_d = ssd/n, σ = √(σ̂²_d/2). */
const mleSigma = (ssd, n) => Math.sqrt(ssd / n / 2);
