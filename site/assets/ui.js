"use strict";
/* ui.js — reusable page controls. The labelled slider with a live readout is
   the one interactive widget every demo page shares. */

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
