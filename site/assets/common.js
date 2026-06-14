"use strict";
/* Shared machinery for the three pages: palette, seeded randomness, the
   Student-t distribution (implemented natively — no library needed), and
   the slider builder. Pages keep their own geometry and rendering. */

const PAL = {
  blue: "#3D6FB0", orange: "#E0883A",
  ink: "#333333", muted: "#888888", rule: "#DDDDDD", pale: "#BFBFBF",
  tail: "#C96A5E", green: "#2E7D5B",
};

const POOL = 100; // fixed pool of random draws — dots don't jitter while you drag

/* ---- seeded randomness ------------------------------------------------ */

function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

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

/* ---- small numerics ---------------------------------------------------- */

function linspace(a, b, n) {
  const out = new Array(n);
  for (let i = 0; i < n; i++) out[i] = a + (b - a) * i / (n - 1);
  return out;
}

const mean = (a) => a.reduce((s, v) => s + v, 0) / a.length;

function svar(a) { // sample variance, ddof = 1
  const m = mean(a);
  return a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1);
}

/* ---- Student's t distribution ------------------------------------------ */

function logGamma(x) {
  const g = [76.18009172947146, -86.50532032941677, 24.01409824083091,
             -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let y = x, tmp = x + 5.5;
  tmp -= (x + 0.5) * Math.log(tmp);
  let ser = 1.000000000190015;
  for (let j = 0; j < 6; j++) ser += g[j] / ++y;
  return -tmp + Math.log(2.5066282746310005 * ser / x);
}

function betacf(a, b, x) {
  const FPMIN = 1e-300, EPS = 3e-12;
  const qab = a + b, qap = a + 1, qam = a - 1;
  let c = 1, d = 1 - qab * x / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= 200; m++) {
    const m2 = 2 * m;
    let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d; h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < EPS) break;
  }
  return h;
}

function ibeta(a, b, x) { // regularized incomplete beta I_x(a, b)
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const bt = Math.exp(logGamma(a + b) - logGamma(a) - logGamma(b)
                      + a * Math.log(x) + b * Math.log(1 - x));
  return x < (a + 1) / (a + b + 2)
    ? bt * betacf(a, b, x) / a
    : 1 - bt * betacf(b, a, 1 - x) / b;
}

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

const normPdf = (x) => Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);

/* standard normal CDF via erf (Abramowitz & Stegun 7.1.26, |err| < 1.5e-7) —
   enough to display a teaching-grade power/β computed from the noncentral t. */
function erf(x) {
  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x);
  const t = 1 / (1 + 0.3275911 * x);
  const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t
                  - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return sign * y;
}
const normCdf = (x) => 0.5 * (1 + erf(x / Math.SQRT2));

/* two-sample t-test, pooled variance (same as scipy.stats.ttest_ind) */
function tTest(xa, xb) {
  const n = xa.length, m = xb.length;
  const sp2 = ((n - 1) * svar(xa) + (m - 1) * svar(xb)) / (n + m - 2);
  const t = (mean(xb) - mean(xa)) / Math.sqrt(sp2 * (1 / n + 1 / m));
  const df = n + m - 2;
  return { t, df, p: 2 * (1 - tCdf(Math.abs(t), df)) };
}

/* ---- UI: labelled slider with live readout ------------------------------ */

function makeSlider(parent, { label, min, max, step, value, fmt, onInput }) {
  const row = document.createElement("label");
  row.className = "ctl";
  row.innerHTML =
    `<span class="ctl-label">${label}</span>` +
    `<input type="range" min="${min}" max="${max}" step="${step}" value="${value}">` +
    `<span class="ctl-readout"></span>`;
  const input = row.querySelector("input");
  const out = row.querySelector(".ctl-readout");
  const upd = () => { out.textContent = fmt(+input.value); };
  input.addEventListener("input", () => { upd(); onInput(); });
  upd();
  parent.appendChild(row);
  return {
    get value() { return +input.value; },
    set value(v) { input.value = v; upd(); },
  };
}
