"use strict";
/* ui.js — reusable page controls and the small translators from slider values to
   model inputs. The labelled slider with a live readout is the one interactive
   widget every demo page shares. */

/* the σ²₀ slider → the variance prior's Inv-Gamma (a₀, b₀), shared by steps 3 & 4.
   The shape is DECOUPLED from κ₀ (which then weights only the mean) and held fixed at
   a₀ = ν₀/2, with ν₀ a gentle pseudo-count kept > 2 so the mean stays finite; b₀ pins
   the EXPECTED variance E[σ²] = b₀/(a₀−1)/2 to the group variance s0 (the ×2 maps
   σ² → σ²_d) — no a₀ = 1 singularity. */
const NIG_NU0 = 4;
function variancePrior(s0) {
  const a0 = NIG_NU0 / 2;
  return { a0, b0: (a0 - 1) * 2 * s0 };
}

function makeSlider(parent, { label, min, max, step, value, fmt, onInput }) {
  const row = document.createElement("label");
  row.className = "ctl";
  row.innerHTML =
    `<span class="ctl-label">${label}</span>` +
    `<input type="range" min="${min}" max="${max}" step="${step}" value="${value}">` +
    `<span class="ctl-readout"></span>`;
  const input = row.querySelector("input");
  input.dataset.init = value;            // remembered so a reset button can restore it
  const out = row.querySelector(".ctl-readout");
  const upd = () => { out.textContent = fmt(+input.value); };
  input.addEventListener("input", () => { upd(); onInput(); });
  upd();
  parent.appendChild(row);
  return {
    get value() { return +input.value; },
    set value(v) { input.value = v; upd(); },
    reset() { input.value = input.dataset.init; input.dispatchEvent(new Event("input")); },
  };
}

/* Wire a reset button: restore every slider inside `container` to its initial
   value and let each fire its own onInput so the views redraw. No-op while a
   simulation has the controls disabled. */
function wireReset(btn, container) {
  if (!btn) return;
  btn.addEventListener("click", () => {
    if (container.querySelector("input[type=range]:disabled")) return;
    container.querySelectorAll("input[type=range]").forEach(inp => {
      inp.value = inp.dataset.init;
      inp.dispatchEvent(new Event("input"));
    });
  });
}

/* Presenter mode: a distraction-free, full-width stage for talks. Shift+P
   toggles it, Esc leaves it. The styling lives in style.css under `.presenter`;
   here we only flip the body class and flash a brief confirmation toast. */
(function presenterMode() {
  let toast;
  const flash = (msg) => {
    if (!toast) {
      toast = document.createElement("div");
      toast.className = "presenter-toast";
      document.body.appendChild(toast);
    }
    toast.textContent = msg;
    toast.classList.add("show");
    clearTimeout(flash._t);
    flash._t = setTimeout(() => toast.classList.remove("show"), 1100);
  };

  const set = (on) => {
    document.body.classList.toggle("presenter", on);
    flash(on ? "Presenter mode · Esc to exit" : "Presenter mode off");
  };

  document.addEventListener("keydown", (e) => {
    if (e.metaKey || e.ctrlKey || e.altKey) return;       // leave browser shortcuts alone
    if (e.shiftKey && e.code === "KeyP") {
      e.preventDefault();
      set(!document.body.classList.contains("presenter"));
    } else if (e.key === "Escape" && document.body.classList.contains("presenter")) {
      set(false);
    }
  });
})();
