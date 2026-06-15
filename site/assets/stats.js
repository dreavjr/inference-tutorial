"use strict";
/* stats.js — distributions and statistical inference built on math.js:
   the Student-t family (native, no library), the normal CDF, the two-sample
   t-test, sampling pools, and the Normal-inverse-gamma conjugate model shared
   by the prior/decision pages. */

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

/* ---- normal distribution ----------------------------------------------- */

const normPdf = (x) => Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);
const normCdf = (x) => 0.5 * (1 + erf(x / Math.SQRT2));

/* ---- two-sample t-test, pooled variance (= scipy.stats.ttest_ind) ------- */

function tTest(xa, xb) {
  const n = xa.length, m = xb.length;
  const sp2 = ((n - 1) * svar(xa) + (m - 1) * svar(xb)) / (n + m - 2);
  const t = (mean(xb) - mean(xa)) / Math.sqrt(sp2 * (1 / n + 1 / m));
  const df = n + m - 2;
  return { t, df, p: 2 * (1 - tCdf(Math.abs(t), df)) };
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
