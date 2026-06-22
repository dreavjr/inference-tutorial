"use strict";
/* plots.js — shared data-visualisation: the colour palette, the common 1-D
   plotting geometry used by steps 1, 3, 4 & 5, the curve/SVG helpers and reusable
   chrome (the top-fade mask, the verdict pill, the splayed Δ markers), the two
   reusable figures (the populations panel and the prior·likelihood·posterior
   marginal), and the repeated-sampling interval simulation + hover wiring shared
   by steps 1 & 3. Pages keep only the figures that are unique to them. */

const PAL = {
    blue: "#3D6FB0", orange: "#E0883A",
    ink: "#333333", muted: "#888888", rule: "#DDDDDD", pale: "#BFBFBF",
    tail: "#C96A5E", green: "#2E7D5B",
};

/* the three voices of Bayesian updating — prior, likelihood, posterior */
const PRIOR_C = "#8E6FAE", LIK_C = "#C96A5E", POST_C = "#333333";

/* ---- shared 1-D geometry ----------------------------------------------- */
/* One scale across the pages: SC px per x-unit (width), K px per density-unit
   (height). Box is wider than tall so curves fill the frame. */
const S = 460, SH = 300, HALF = 5.0, YMAX = 0.9;
const SC = S / (2 * HALF), K = SH / YMAX;
const ML = 16, MR = 16, BH = 30;
const FW = ML + S + MR, X0 = ML + S / 2;
const TOP = 14, BASE = TOP + SH, FH = BASE + BH;

const Xp = (x) => X0 + x * SC;
const Yp = (base, dens, cap = SH) => base - Math.min(dens * K, cap);
const clampx = (px, pad = 28) => Math.min(Math.max(px, ML + pad), ML + S - pad);

/* ---- SVG helpers ------------------------------------------------------- */

// "x,y x,y …" point string for a density curve, optionally clamped at `cap` px
function curvePath(xs, dens, base, cap = S) {
    return xs.map((x, i) => `${Xp(x).toFixed(1)},${Yp(base, dens[i], cap).toFixed(1)}`).join(" ");
}

// filled area under a curve + its outline, as one string (steps 1 & 3)
function filledCurve(xs, dens, base, stroke, fill, op, sw = 2) {
    const pts = xs.map((x, i) => [Xp(x), Yp(base, dens[i])]);
    const line = pts.map(p => `${p[0].toFixed(1)},${p[1].toFixed(1)}`).join(" ");
    const poly = `${pts[0][0].toFixed(1)},${base} ${line} ${pts[pts.length - 1][0].toFixed(1)},${base}`;
    return `<polygon points='${poly}' fill='${fill}' opacity='${op}'/>` +
        `<polyline points='${line}' fill='none' stroke='${stroke}' stroke-width='${sw}'/>`;
}

// client coordinates → SVG user space
function svgPoint(svg, e) {
    const pt = svg.createSVGPoint();
    pt.x = e.clientX; pt.y = e.clientY;
    return pt.matrixTransform(svg.getScreenCTM().inverse());
}

/* A luminance mask that dissolves a curve as it climbs out the top of the box, so
   a band tighter/taller than the frame fades instead of clipping flat. The gradient
   runs transparent at `top` → opaque by `soft` of the way down to `bottom`; the mask
   itself covers the whole box (down to BASE) so nothing below the fade is clipped.
   Returns the <defs> to drop in and the `mask=` reference for the group it clips;
   `id` must be unique among masks live in the same SVG document. */
function fadeTopMask(id, { top = TOP, bottom = BASE, soft = 0.22 } = {}) {
    const g = `fadeG-${id}`, m = `fadeM-${id}`;
    return {
        defs: `<defs><linearGradient id='${g}' gradientUnits='userSpaceOnUse' x1='0' y1='${top}' x2='0' y2='${bottom}'>` +
            `<stop offset='0' stop-color='#fff' stop-opacity='0'/>` +
            `<stop offset='${soft}' stop-color='#fff' stop-opacity='1'/>` +
            `<stop offset='1' stop-color='#fff' stop-opacity='1'/></linearGradient>` +
            `<mask id='${m}'><rect x='0' y='${top}' width='${FW}' height='${BASE - top}' fill='url(#${g})'/></mask></defs>`,
        mask: `url(#${m})`,
    };
}

/* a verdict/score badge: a rounded pill sized to its label, white text centred */
function pill(x, y, text, bg) {
    const w = text.length * 6.6 + 23;
    return `<rect x='${x}' y='${y}' width='${w.toFixed(0)}' height='22' rx='11' fill='${bg}'/>` +
        `<text x='${(x + w / 2).toFixed(1)}' y='${y + 15}' text-anchor='middle' ` +
        `font-size='12' font-weight='700' fill='#fff'>${text}</text>`;
}

/* a reference marker on a 1-D panel: a dashed vertical line from the baseline up to
   the label, plus the label itself. `side` splays the text left (−1) / right (+1) /
   centred (0) so neighbouring markers don't overprint. */
function splayMarker(xv, color, label, side, { lineTopY, textY }) {
    const x = Xp(xv);
    const anchor = side < 0 ? "end" : side > 0 ? "start" : "middle";
    const tx = side ? Math.min(Math.max(x + side * 5, ML + 12), ML + S - 12) : clampx(x);
    return `<line x1='${x.toFixed(1)}' x2='${x.toFixed(1)}' y1='${BASE}' y2='${lineTopY}' ` +
        `stroke='${color}' stroke-width='1.4' stroke-dasharray='5 4'/>` +
        `<text x='${tx.toFixed(1)}' y='${textY}' text-anchor='${anchor}' font-size='13.5' ` +
        `font-weight='600' fill='${color}'>${label}</text>`;
}

/* ---- the two populations and the sample in hand ------------------------
   Optionally lays down the hidden implied-mean guides that steps 1 & 3 reveal
   when their right-hand panel is hovered. */

function populationSvg(delta, varr, n, pool, guides = false) {
    const sigma = Math.sqrt(varr);
    const xs = linspace(-HALF, HALF, 201);
    const p = [];
    for (const [mu, col] of [[-delta / 2, PAL.blue], [delta / 2, PAL.orange]]) {
        const pdf = xs.map(x => Math.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * Math.sqrt(2 * Math.PI)));
        const line = curvePath(xs, pdf, BASE);
        p.push(`<polygon points='${Xp(-HALF).toFixed(1)},${BASE} ${line} ${Xp(HALF).toFixed(1)},${BASE}' fill='${col}' opacity='0.10'/>`);
        p.push(`<polyline points='${line}' fill='none' stroke='${col}' stroke-width='2'/>`);
    }
    p.push(`<line x1='${ML}' x2='${ML + S}' y1='${BASE}' y2='${BASE}' stroke='${PAL.rule}' stroke-width='1'/>`);
    for (const [z, j, mu, col, side] of [[pool.z1, pool.j1, -delta / 2, PAL.blue, -1],
    [pool.z2, pool.j2, delta / 2, PAL.orange, 1]]) {
        for (let i = 0; i < n; i++) {
            const x = Math.max(-HALF, Math.min(HALF, mu + sigma * z[i]));
            const y = BASE + side * (11 + j[i] * 8);
            p.push(`<circle cx='${Xp(x).toFixed(1)}' cy='${y.toFixed(1)}' r='3' fill='${col}' opacity='0.75'/>`);
        }
    }
    if (guides) {
        p.push(`<line id='impL' x1='${X0}' x2='${X0}' y1='${BASE}' y2='${TOP + 10}' ` +
            `stroke='${PAL.ink}' stroke-width='1.2' stroke-dasharray='4 3' opacity='0'/>`);
        p.push(`<line id='impR' x1='${X0}' x2='${X0}' y1='${BASE}' y2='${TOP + 10}' ` +
            `stroke='${PAL.ink}' stroke-width='1.2' stroke-dasharray='4 3' opacity='0'/>`);
        p.push(`<text id='impLbl' x='${X0}' y='${TOP + 16}' text-anchor='middle' font-size='11.5' ` +
            `font-style='italic' fill='${PAL.muted}' opacity='0'>implied means</text>`);
    }
    return `<svg viewBox='0 0 ${FW} ${FH}' xmlns='http://www.w3.org/2000/svg'>${p.join("")}</svg>`;
}

/* ---- the marginal over Δ: prior · likelihood · posterior ----------------
   prior & posterior are genuine densities (area = 1) → filled. The likelihood is
   a relative-support weighting (peak = 1): pinned to REF of the panel, drawn
   unfilled and dashed — no area to read. Curves fade out over the top of the box. */

function marginalSvg(prior, lik, post, delta) {
    const p = [];
    const fade = fadeTopMask("marg");
    p.push(fade.defs);

    const xs = linspace(-HALF, HALF, 321);
    const dPrior = xs.map(prior), dPost = xs.map(post);
    const REF = 0.72, refDens = REF * YMAX, refY = Yp(BASE, refDens);
    const dLik = xs.map(x => lik(x) * refDens);

    const filled = (dens, col, op) => {
        const line = curvePath(xs, dens, BASE);
        return `<polygon points='${Xp(-HALF).toFixed(1)},${BASE} ${line} ${Xp(HALF).toFixed(1)},${BASE}' fill='${col}' opacity='${op}'/>`;
    };
    p.push(`<g mask='${fade.mask}'>` +
        filled(dPrior, PRIOR_C, 0.10) +
        filled(dPost, POST_C, 0.12) +
        `<polyline points='${curvePath(xs, dPrior, BASE)}' fill='none' stroke='${PRIOR_C}' stroke-width='2'/>` +
        `<polyline points='${curvePath(xs, dPost, BASE)}' fill='none' stroke='${POST_C}' stroke-width='2.2'/>` +
        `<polyline points='${curvePath(xs, dLik, BASE)}' fill='none' stroke='${LIK_C}' stroke-width='1.8' stroke-dasharray='6 3'/>` +
        `</g>`);
    p.push(`<line x1='${ML}' x2='${ML + S}' y1='${BASE}' y2='${BASE}' stroke='${PAL.rule}' stroke-width='1'/>`);

    p.push(`<line x1='${ML}' x2='${ML + S}' y1='${refY.toFixed(1)}' y2='${refY.toFixed(1)}' ` +
        `stroke='${LIK_C}' stroke-width='1' stroke-dasharray='2 4' opacity='0.65'/>`);
    p.push(`<text x='${ML + 6}' y='${(refY - 5).toFixed(1)}' font-size='10.5' font-style='italic' ` +
        `fill='${LIK_C}'>likelihood: support = 1 at peak</text>`);

    const mark = (xv, color, lbl, side) =>
        splayMarker(xv, color, lbl, side, { lineTopY: TOP + 21, textY: TOP + 16 });
    p.push(mark(0, PAL.muted, "0", delta > 0.05 ? -1 : 0));
    p.push(mark(delta, PAL.green, "Δ", 1));

    const lineX2 = ML + S, lineX1 = lineX2 - 24, labelX = lineX1 - 9;
    const item = (y, col, txt, sw = 2, dash = "") =>
        `<line x1='${lineX1}' x2='${lineX2}' y1='${y}' y2='${y}' stroke='${col}' stroke-width='${sw}'` +
        `${dash ? ` stroke-dasharray='${dash}'` : ""}/>` +
        `<text x='${labelX}' y='${y + 4}' text-anchor='end' font-size='12' fill='${PAL.ink}'>${txt}</text>`;
    p.push(item(TOP + 34, PRIOR_C, "prior"));
    p.push(item(TOP + 53, LIK_C, "likelihood (relative)", 1.8, "6 3"));
    p.push(item(TOP + 72, POST_C, "posterior", 2.2));

    return `<svg viewBox='0 0 ${FW} ${FH}' xmlns='http://www.w3.org/2000/svg'>${p.join("")}</svg>`;
}

/* ---- implied-mean hover (steps 1 & 3) ----------------------------------
   The right-hand panel (null density / posterior) reads as Δ along x; hovering it
   drops a guide to the curve, prints the value, and projects the implied ±Δ/2
   population means onto the left panel. Pages supply the curve-height and label
   formulas for their own statistic. */

function wireImpliedHover(svgId, { heightAt, label, labelX, labelMinY }) {
    const svg = document.getElementById(svgId);
    const hov = document.getElementById("hoverLine"), hl = document.getElementById("hoverLbl");
    const iL = document.getElementById("impL"), iR = document.getElementById("impR");
    const iLbl = document.getElementById("impLbl");
    const hide = () => [hov, hl, iL, iR, iLbl].forEach(el => { el.style.opacity = 0; });

    svg.addEventListener("pointermove", (e) => {
        const u = svgPoint(svg, e);
        if (u.x < ML || u.x > ML + S || u.y < TOP || u.y > BASE + 20) { hide(); return; }
        const x = u.x, d = (x - X0) / SC, cy = heightAt(d);
        hov.setAttribute("x1", x); hov.setAttribute("x2", x);
        hov.setAttribute("y1", BASE); hov.setAttribute("y2", cy);
        hov.style.opacity = 1;
        const ly = Math.min(Math.max(cy - 8, labelMinY), BASE - 6);
        hl.textContent = label(d);
        hl.setAttribute("x", labelX(x)); hl.setAttribute("y", ly);
        hl.style.opacity = 1;
        const xl = X0 - (d / 2) * SC, xr = X0 + (d / 2) * SC;
        iL.setAttribute("x1", xl); iL.setAttribute("x2", xl); iL.style.opacity = 1;
        iR.setAttribute("x1", xr); iR.setAttribute("x2", xr); iR.style.opacity = 1;
        iLbl.style.opacity = 1;
    });
    svg.addEventListener("pointerleave", hide);
}

/* ---- repeated-sampling interval simulation (steps 1 & 3) ----------------
   Draws N_SIM repeated samples as thick vertical interval bars (a confidence
   interval in step 1, a credible HDI in step 3 — numerically the same under flat
   priors). Colour encodes significance, dash encodes a "false" verdict; the run
   animates the panels above through each draw, then makes the bars hoverable so
   pointing at one restores that exact sample. Only the caption and the
   "false" predicate differ between the pages.

   When `tTest` is supplied (step 1) a second panel is stacked above the interval
   strip, sharing the same per-sample data and synchronised hover: it plots each
   sample's |t| as a bar rising from t = 0, against one dashed critical-value line
   (df = 2n−2 and α are fixed across draws, so the threshold is constant). A bar
   crossing the threshold is significant — exactly when the interval clears 0. */

function makeIntervalSim({ sliders, simArea, simBtn, resampleBtn, getParams, adopt, caption, isFalse, tTest = null }) {
    const SIMW = 960, SIMH = 165, SML = 46, SMR = 24, SMT = 16, SMB = 16;
    const N_SIM = 40, SIM_WAIT = 150;
    const MAROON = "#8C2F39";
    const simMidY = SIMH / 2, simSpacing = (SIMW - SML - SMR) / N_SIM;
    const NS = "http://www.w3.org/2000/svg";
    const sleep = (ms) => new Promise(r => setTimeout(r, ms));

    // stacked t-statistic panel (tTest mode only): t = 0 baseline near the bottom,
    // bars climb to |t|; same width/x-grid as the interval strip so columns align.
    const TIMH = 150, TMT = 16, TMB = 26;
    const tBase = TIMH - TMB, tPlotTop = TMT, tPlotH = tBase - tPlotTop;

    let simRunning = false;
    let simEntries = [];

    function setControlsDisabled(d) {
        sliders.querySelectorAll("input").forEach(i => { i.disabled = d; });
        resampleBtn.disabled = d;
        simBtn.disabled = d;
    }

    function dismiss() {
        if (simRunning) return;
        simEntries = [];
        simArea.hidden = true;
        simArea.innerHTML = "";
    }

    async function run() {
        if (simRunning) return;
        dismiss();
        simRunning = true;
        setControlsDisabled(true);

        const { delta, varr, n, alpha } = getParams();
        const sigma = Math.sqrt(varr);
        const df = 2 * n - 2;
        const tc = tInv(1 - alpha / 2, df);   // two-sided critical value, constant across draws

        // draw all samples up front so the vertical scale can be fixed sensibly
        const draws = [];
        for (let i = 0; i < N_SIM; i++) {
            const p = makePool();
            const xa = p.z1.slice(0, n).map(z => -delta / 2 + sigma * z);
            const xb = p.z2.slice(0, n).map(z => delta / 2 + sigma * z);
            const dObs = mean(xb) - mean(xa);
            const sp = Math.sqrt((svar(xa) + svar(xb)) / 2);
            const se = sp * Math.sqrt(2 / n);
            const half = tc * se;                                // (1−α) interval half-width
            draws.push({ pool: p, dObs, lo: dObs - half, hi: dObs + half, t: dObs / se });
        }

        // robust scale: 92nd-percentile spread around Δ, always wide enough to keep
        // the 0-line in frame
        const devs = draws.flatMap(d => [Math.abs(d.lo - delta), Math.abs(d.hi - delta)]).sort((a, b) => a - b);
        const spread = devs[Math.floor(devs.length * 0.92)] * 1.06;
        const ymax = Math.max(spread, Math.abs(delta) + 0.5, 0.5);
        const yScale = (simMidY - SMT) / ymax;
        const clampY = (d) => simMidY - Math.max(-(simMidY - SMT), Math.min(simMidY - SMT, (d - delta) * yScale));

        const trueY = clampY(delta).toFixed(1);
        const zeroY = clampY(0).toFixed(1);
        const x0 = SML, x1 = SIMW - SMR;
        const shell = [];

        // intervals go down first → guides, labels and legend paint ON TOP of them
        shell.push(`<g id='simLines'></g>`);
        shell.push(`<line id='simMarker' x1='0' x2='0' y1='${SMT}' y2='${SIMH - SMB}' stroke='${PAL.ink}' ` +
            `stroke-width='1' stroke-dasharray='3 3' opacity='0'/>`);
        shell.push(`<line x1='${x0}' x2='${x1}' y1='${trueY}' y2='${trueY}' stroke='${PAL.green}' ` +
            `stroke-width='1.2' stroke-dasharray='4 4' opacity='0.7'/>`);
        shell.push(`<line x1='${x0}' x2='${x1}' y1='${zeroY}' y2='${zeroY}' stroke='${PAL.ink}' ` +
            `stroke-width='1.2' stroke-dasharray='4 4' opacity='0.7'/>`);
        shell.push(`<text x='${x1}' y='${(+trueY - 5).toFixed(1)}' text-anchor='end' font-size='11' ` +
            `fill='${PAL.green}'>true Δ</text>`);
        shell.push(`<text x='${x1}' y='${(+zeroY + 12).toFixed(1)}' text-anchor='end' font-size='11' ` +
            `fill='${PAL.ink}'>0</text>`);

        // legend, bottommost row (also on top): colour = significance, style = correctness
        const lgY = SIMH - 6;
        let lx = x0;
        for (const [col, dash, label] of [[MAROON, false, "significant"], [PAL.pale, false, "not significant"],
        [PAL.ink, false, "significance matches Δ being inside interval"], [PAL.ink, true, "significance mismatches Δ coverage"]]) {
            shell.push(`<line x1='${lx}' x2='${lx + 16}' y1='${lgY - 3}' y2='${lgY - 3}' stroke='${col}' ` +
                `stroke-width='3' stroke-linecap='round'${dash ? " stroke-dasharray='0 6'" : ""}/>`);
            shell.push(`<text x='${lx + 22}' y='${lgY}' font-size='10.5' fill='${PAL.muted}'>${label}</text>`);
            lx += 22 + label.length * 5.8 + 26;
        }

        // optional top panel: |t| bars against the constant critical-value line
        let topShell = null, yT = null;
        if (tTest) {
            const absT = draws.map(d => Math.abs(d.t)).sort((a, b) => a - b);
            const tSpread = absT[Math.floor(absT.length * 0.92)] * 1.08;
            const tmax = Math.max(tSpread, tc * 1.3, 1);
            yT = (v) => tBase - Math.min(tPlotH, Math.max(0, v) * (tPlotH / tmax));
            const yc = yT(tc).toFixed(1);
            topShell = [];
            topShell.push(`<g id='simTLines'></g>`);
            topShell.push(`<line id='simTMarker' x1='0' x2='0' y1='${tPlotTop}' y2='${tBase}' stroke='${PAL.ink}' ` +
                `stroke-width='1' stroke-dasharray='3 3' opacity='0'/>`);
            topShell.push(`<line x1='${x0}' x2='${x1}' y1='${tBase}' y2='${tBase}' stroke='${PAL.rule}' stroke-width='1'/>`);
            topShell.push(`<line x1='${x0}' x2='${x1}' y1='${yc}' y2='${yc}' stroke='${MAROON}' ` +
                `stroke-width='1.2' stroke-dasharray='4 4' opacity='0.75'/>`);
            topShell.push(`<text x='${x1}' y='${(+yc - 5).toFixed(1)}' text-anchor='end' font-size='11' ` +
                `fill='${MAROON}'>critical value (p-value = α)</text>`);
            topShell.push(`<text x='${x1}' y='${(tBase - 5).toFixed(1)}' text-anchor='end' font-size='11' ` +
                `fill='${PAL.ink}'>t = 0</text>`);
            topShell.push(`<text x='${(SML - 8).toFixed(1)}' y='${(tPlotTop + 9).toFixed(1)}' text-anchor='end' ` +
                `font-size='9.5' fill='${PAL.muted}'>|t|</text>`);
        }

        simArea.hidden = false;
        if (tTest) {
            simArea.innerHTML =
                `<div class="figcap">3 · repeated-sampling simulation</div>` +
                `<div class="sim-stack">` +
                `<div class="sim-subcap">${tTest.caption}</div>` +
                `<svg id='simTSvg' viewBox='0 0 ${SIMW} ${TIMH}' xmlns='${NS}'>${topShell.join("")}</svg>` +
                `<div class="sim-subcap">${caption(alpha)}</div>` +
                `<svg id='simSvg' viewBox='0 0 ${SIMW} ${SIMH}' xmlns='${NS}'>${shell.join("")}</svg>` +
                `</div>`;
        } else {
            simArea.innerHTML =
                `<div class="figcap">3 · ${caption(alpha)}</div>` +
                `<svg id='simSvg' viewBox='0 0 ${SIMW} ${SIMH}' xmlns='${NS}'>${shell.join("")}</svg>`;
        }
        const svg = simArea.querySelector("#simSvg");
        const gLines = svg.querySelector("#simLines");
        const tSvg = tTest ? simArea.querySelector("#simTSvg") : null;
        const gTopLines = tSvg ? tSvg.querySelector("#simTLines") : null;

        simEntries = [];
        for (let i = 0; i < N_SIM; i++) {
            const d = draws[i];

            // (1) adopt this sample and refresh the panels above
            adopt(d.pool);

            // (2) lay down its interval as a thick vertical line
            const x = SML + (i + 0.5) * simSpacing;
            const significant = !(d.lo <= 0 && d.hi >= 0);     // interval clears 0 → reject H₀
            const coversTrue = d.lo <= delta && d.hi >= delta; // does the interval capture true Δ?
            const line = document.createElementNS(NS, "line");
            line.setAttribute("x1", x.toFixed(1)); line.setAttribute("x2", x.toFixed(1));
            line.setAttribute("y1", clampY(d.lo).toFixed(1)); line.setAttribute("y2", clampY(d.hi).toFixed(1));
            line.setAttribute("stroke", significant ? MAROON : PAL.pale);
            line.setAttribute("stroke-width", "3"); line.setAttribute("stroke-linecap", "round");
            if (isFalse(significant, coversTrue)) line.setAttribute("stroke-dasharray", "0 6");
            line.setAttribute("data-idx", i);
            gLines.appendChild(line);
            const dot = document.createElementNS(NS, "circle");
            dot.setAttribute("cx", x.toFixed(1)); dot.setAttribute("cy", clampY(d.dObs).toFixed(1));
            dot.setAttribute("r", "2"); dot.setAttribute("fill", PAL.ink);
            gLines.appendChild(dot);

            // (2b) and its |t| bar on the stacked panel, same colour encoding
            if (gTopLines) {
                const tline = document.createElementNS(NS, "line");
                tline.setAttribute("x1", x.toFixed(1)); tline.setAttribute("x2", x.toFixed(1));
                tline.setAttribute("y1", tBase.toFixed(1)); tline.setAttribute("y2", yT(Math.abs(d.t)).toFixed(1));
                tline.setAttribute("stroke", significant ? MAROON : PAL.pale);
                tline.setAttribute("stroke-width", "3"); tline.setAttribute("stroke-linecap", "round");
                tline.setAttribute("data-idx", i);
                gTopLines.appendChild(tline);
            }
            simEntries.push({ pool: d.pool, x });

            // (3) hold for a beat
            await sleep(SIM_WAIT);
        }

        simRunning = false;
        setControlsDisabled(false);
        wireSimHover(svg, gLines, tSvg, gTopLines);
    }

    // once the run is complete, hovering either panel restores its sample above and
    // highlights the matching column in both
    function wireSimHover(svg, gLines, tSvg, gTopLines) {
        const marker = svg.querySelector("#simMarker");
        const tMarker = tSvg ? tSvg.querySelector("#simTMarker") : null;
        let last = -1;
        const setW = (g, i, w) => { if (!g) return; const l = g.querySelector(`[data-idx='${i}']`); if (l) l.setAttribute("stroke-width", w); };
        const restore = (i) => { setW(gLines, i, "3"); setW(gTopLines, i, "3"); };
        const onMove = (src, e) => {
            if (simRunning) return;
            const idx = Math.floor((svgPoint(src, e).x - SML) / simSpacing);
            if (idx < 0 || idx >= simEntries.length || idx === last) return;
            if (last >= 0) restore(last);
            last = idx;
            setW(gLines, idx, "5"); setW(gTopLines, idx, "5");
            const x = simEntries[idx].x;
            marker.setAttribute("x1", x); marker.setAttribute("x2", x); marker.style.opacity = 1;
            if (tMarker) { tMarker.setAttribute("x1", x); tMarker.setAttribute("x2", x); tMarker.style.opacity = 1; }
            adopt(simEntries[idx].pool);
        };
        const onLeave = () => {
            if (last >= 0) restore(last);
            marker.style.opacity = 0;
            if (tMarker) tMarker.style.opacity = 0;
            last = -1;
        };
        for (const s of [svg, tSvg]) {
            if (!s) continue;
            s.addEventListener("pointermove", (e) => onMove(s, e));
            s.addEventListener("pointerleave", onLeave);
            s.style.cursor = "pointer";
        }
    }

    simBtn.addEventListener("click", run);
    return { run, dismiss, get running() { return simRunning; } };
}
