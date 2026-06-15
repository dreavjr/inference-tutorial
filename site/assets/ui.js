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
