# Refactoring pass — math/display separation + dedup

Goal: complete the **statistical-model layer** so each page keeps only display +
wiring; lift duplicated **display helpers** into `plots.js`; isolate page-5's
"reward sweep" easter egg from the core decision drawing; align naming.
Rendering must stay visually identical (verify in browser at the end).

Stack (unchanged): `math.js` → `stats.js` → `plots.js` → `ui.js`. Pages load all four.

## Plan & status

Legend: [ ] todo · [~] in progress · [x] done · [!] needs verify

### 1. stats.js — new composed-model layer  [x] DONE
- [x] `drawSample(pool, {delta, sigma, n})` → `{xa, xb, d, dbar, ssd}`
- [x] extend `tTest` to also return `{dObs, se, sp}` (t/df/p unchanged)
- [x] `tCrit(alpha, df)` = `tInv(1−alpha/2, df)`
- [x] `tTestPower({delta, sigma, n, alpha, df})` → `{power, beta}`
- [x] `locScaleTpdf(loc, scale, df)`
- [x] `nigMarginalDelta(m, k, a, b)` → `{loc, scale, df}`
- [x] `likelihoodProfile(dbar, ssd, n)` → peak-1 Gaussian profile
- [x] `nigSigmaMode(a, b)` and `mleSigma(ssd, n)`

### 2. plots.js — shared display helpers  [x] DONE
- [x] `fadeTopMask(id, {top, bottom, soft})` → `{defs, mask}`
- [x] `pill(x, y, text, bg)` → rounded-rect + centred white label
- [x] `splayMarker(xv, color, label, side, {lineTopY, textY})`
- [x] marginalSvg now uses fadeTopMask + splayMarker
- [ ] (page 5 uses the shared `svgPoint`, drop its local copy)

### 3. Pages — use the model + helpers, align naming (main fn = `render`)
- [x] 01-the-t-test.html  (refresh→render; drawSample/tTest/tCrit/tTestPower; fadeTopMask/pill) — verified
- [x] 03-bayesian-t-test.html (drawSample/tTest/tCrit; fadeTopMask/pill/splayMarker) — verified
- [x] 04-the-prior.html (drawSample diffs; nigMarginalDelta/locScaleTpdf/likelihoodProfile; nigSigmaMode/mleSigma) — verified
- [x] 05-the-decision.html (computeState via model layer; dropped local svgPoint; fadeTopMask/shared pill;
      egg→sweep, posterior belief region extracted to posteriorBelief() out of core rewardSvg) — verified
      NOTE: fadeTopMask mask now always covers down to BASE (gradient end ≠ mask end) — needed for reward panels.

### 4. Verify  [x] DONE — all via Claude Preview (python http.server on :8765), console clean
- [x] 01: pill + p-value + β render; 03: splayed Δ/Δ̂/0 markers + risk pill
- [x] 04: prior/lik/posterior contours, modes, marginal triptych
- [x] 05: both reward panels (lens, CI bar, dotted links, pills); 40-draw simulation
      (both panels + summary stats); double-click sweep egg (maroon curve, integral 0→accumulating);
      slider drag re-renders + dismisses sim; no console errors
- [x] 02 (untouched) still loads clean against the changed plots.js

## OUTCOME: complete. All four todo sections done and browser-verified.

## Notes / gotchas
- 01 & 03 sampling and 04 & 05 sampling are the SAME draw (verified): 04's
  `Δ + σ(z₂−z₁)` == 05's `xb − xa`. So one `drawSample` serves all.
- `tTest` balanced-n identity: `sp = √((svar_a+svar_b)/2)`, `se = sp√(2/n)` equals
  the general pooled `√(sp²(1/n+1/m))`. Keep the general form in tTest.
- `nigSigmaMode` uses `a+1.5` (= a+3/2): at Δ=m the σ²_d slice is Inv-Gamma(a+½,b)
  (extra ½ from the Normal's σ⁻¹), mode b/(a+3/2); σ = √(mode/2).
- Pill metrics differ trivially between pages (h 22 vs 23, font 12 vs 12.5).
  Unifying to one helper is intentional; visual diff is sub-pixel-ish — confirm in browser.
- The "egg" is a display-only easter egg (sweeps the reward centre on dblclick on
  figRewPost). Keep it page-5-local but stop threading it through the core fn.
