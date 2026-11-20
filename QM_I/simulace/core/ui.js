/* Shared UI building blocks. Every engine builds its controls from these so all
   simulations look and behave the same. */
(function () {
  'use strict';
  const QM = (window.QM = window.QM || {});

  function el(tag, attrs, ...children) {
    const e = document.createElement(tag);
    if (attrs) for (const [k, v] of Object.entries(attrs)) {
      if (v === undefined || v === null || v === false) continue;
      if (k === 'class') e.className = v;
      else if (k === 'style' && typeof v === 'object') Object.assign(e.style, v);
      else if (k.startsWith('on') && typeof v === 'function') e.addEventListener(k.slice(2), v);
      else if (k === 'html') e.innerHTML = v;
      else e.setAttribute(k, v === true ? '' : v);
    }
    children.flat().forEach((c) => { if (c === null || c === undefined || c === false) return; e.append(c instanceof Node ? c : document.createTextNode(String(c))); });
    return e;
  }

  /** render $..$ / $$..$$ with KaTeX auto-render if it loaded; otherwise strip delimiters */
  function tex(root) {
    if (window.renderMathInElement) {
      try {
        window.renderMathInElement(root, { delimiters: [{ left: '$$', right: '$$', display: true }, { left: '$', right: '$', display: false }], throwOnError: false });
        return;
      } catch (e) { /* fall through */ }
    }
    // offline: KaTeX not available – the raw $…$ source stays readable as plain text
  }

  /** layout: two columns (controls | views). returns {controls, views} */
  function layout(root) {
    const controls = el('div', { class: 'lab-controls' });
    const views = el('div', { class: 'lab-views' });
    root.append(el('div', { class: 'lab' }, controls, views));
    return { controls, views };
  }
  function panel(title, opts = {}) {
    const body = el('div', { class: 'panel-body' });
    const head = title ? el('div', { class: 'panel-head' }, el('span', { html: title }), opts.right || null) : null;
    const p = el('section', { class: 'panel' + (opts.class ? ' ' + opts.class : '') }, head, body);
    return { el: p, body, head };
  }
  function canvas(height = 260, cls = '') { return el('canvas', { class: 'plot ' + cls, style: { height: height + 'px' } }); }

  /** Parameter controls from a spec list.
      spec: {key, label, type: 'range'|'select'|'check'|'number', min, max, step, options:[{value,label}], fmt(v), hint}
      returns {el, values, set(key, v, silent), refresh()} */
  function params(specs, values, onChange) {
    const wrap = el('div', { class: 'params' });
    const setters = {};
    specs.forEach((s) => {
      if (s.type === 'heading') { wrap.append(el('div', { class: 'param-heading', html: s.label })); return; }
      if (values[s.key] === undefined) values[s.key] = s.value !== undefined ? s.value : s.type === 'check' ? false : s.type === 'select' ? s.options[0].value : s.min;
      const lab = el('label', { class: 'param-label', html: s.label });
      let row;
      if (s.type === 'select') {
        const sel = el('select', { class: 'param-select' }, s.options.map((o) => el('option', { value: o.value }, o.label)));
        sel.value = values[s.key];
        sel.addEventListener('change', () => { const o = s.options.find((o) => String(o.value) === sel.value); values[s.key] = o ? o.value : sel.value; onChange(s.key); });
        setters[s.key] = (v) => { sel.value = v; };
        row = el('div', { class: 'param' }, lab, sel);
      } else if (s.type === 'check') {
        const cb = el('input', { type: 'checkbox' }); cb.checked = !!values[s.key];
        cb.addEventListener('change', () => { values[s.key] = cb.checked; onChange(s.key); });
        setters[s.key] = (v) => { cb.checked = !!v; };
        row = el('div', { class: 'param param-check' }, el('label', null, cb, el('span', { html: ' ' + s.label })));
      } else {
        const fmt = s.fmt || ((v) => QM.fmt(v, 3));
        const range = el('input', { type: 'range', min: s.min, max: s.max, step: s.step || (s.max - s.min) / 200 });
        const num = el('input', { type: 'number', class: 'param-num', step: s.step || 'any' });
        range.value = values[s.key]; num.value = +(+values[s.key]).toFixed(6);
        const out = s.fmt ? el('span', { class: 'param-out' }, fmt(values[s.key])) : null;
        const upd = (v, from) => { values[s.key] = v; if (from !== range) range.value = v; if (from !== num) num.value = +(+v).toFixed(6); if (out) out.textContent = fmt(v); onChange(s.key); };
        range.addEventListener('input', () => upd(+range.value, range));
        num.addEventListener('change', () => { const v = parseFloat(num.value.replace(',', '.')); if (isFinite(v)) upd(v, num); });
        setters[s.key] = (v) => { range.value = v; num.value = +(+v).toFixed(6); if (out) out.textContent = fmt(v); };
        row = el('div', { class: 'param' }, el('div', { class: 'param-top' }, lab, out, num), range);
      }
      if (s.hint) row.append(el('div', { class: 'param-hint', html: s.hint }));
      wrap.append(row);
    });
    tex(wrap);
    return {
      el: wrap, values,
      set(key, v, silent) { values[key] = v; if (setters[key]) setters[key](v); if (!silent) onChange(key); },
    };
  }

  function button(label, onClick, cls = '') { return el('button', { class: 'btn ' + cls, onclick: onClick, html: label }); }
  function buttons(...btns) { return el('div', { class: 'btn-row' }, ...btns); }

  /** Complex amplitude editor (Re/Im inputs, auto-normalised).
      opts: {labels, vector, presets:[{label, vector}], onChange(normalisedVector)} */
  function amplitudes(opts) {
    const n = opts.labels.length;
    const wrap = el('div', { class: 'amps' });
    const grid = el('div', { class: 'amp-grid' });
    const inputs = [];
    grid.append(el('span', { class: 'amp-h' }, ''), el('span', { class: 'amp-h' }, 'Re'), el('span', { class: 'amp-h' }, 'Im'));
    for (let i = 0; i < n; i++) {
      const re = el('input', { type: 'number', step: '0.1' }), im = el('input', { type: 'number', step: '0.1' });
      inputs.push([re, im]);
      grid.append(el('span', { class: 'amp-lab', html: opts.labels[i] }), re, im);
      [re, im].forEach((inp) => inp.addEventListener('change', fire));
    }
    const msg = el('div', { class: 'amp-msg' });
    if (opts.presets && opts.presets.length) {
      const sel = el('select', { class: 'param-select' }, el('option', { value: '' }, 'Preset state…'), opts.presets.map((p, i) => el('option', { value: i }, p.label)));
      sel.addEventListener('change', () => { if (sel.value === '') return; set(opts.presets[+sel.value].vector); fire(); sel.value = ''; });
      wrap.append(sel);
    }
    wrap.append(grid, msg);
    function read() { return QM.V.from(inputs.map(([re, im]) => [parseFloat(re.value) || 0, parseFloat(im.value) || 0])); }
    function set(v) { for (let i = 0; i < n; i++) { inputs[i][0].value = +v.re[i].toFixed(4); inputs[i][1].value = +v.im[i].toFixed(4); } }
    function fire() {
      const v = read(), nn = QM.V.norm2(v);
      if (nn < 1e-12) { msg.textContent = 'Zero vector – enter at least one non-zero amplitude.'; return; }
      msg.textContent = Math.abs(nn - 1) > 1e-6 ? `Normalised (‖ψ‖² was ${QM.fmt(nn, 4)}).` : '';
      opts.onChange(QM.V.normalize(v));
    }
    if (opts.vector) set(opts.vector);
    return { el: wrap, set, read };
  }

  /** Play/pause + time slider + speed. opts: {tMax, t0, speed, label, onTime(t), unit} */
  function transport(opts) {
    let t = opts.t0 || 0, tMax = opts.tMax, playing = false, last = 0, raf = 0, speed = opts.speed || tMax / 10;
    const playBtn = button('▶ Play', toggle, 'primary');
    const range = el('input', { type: 'range', min: 0, max: tMax, step: tMax / 1000, value: t });
    const out = el('span', { class: 'param-out' });
    const spd = el('input', { type: 'range', min: -2, max: 1.5, step: 0.05, value: 0 });
    range.addEventListener('input', () => setT(+range.value));
    const box = el('div', { class: 'transport' },
      el('div', { class: 'param-top' }, el('label', { class: 'param-label', html: opts.label || 'Time $t$' }), out),
      range,
      el('div', { class: 'btn-row' }, playBtn, button('⟲ t = 0', () => setT(0)), el('label', { class: 'speed' }, 'speed', spd)));
    tex(box);
    function fmtT() { out.textContent = (opts.fmt ? opts.fmt(t) : QM.fmt(t, 3)) + (opts.unit || ''); }
    function setT(v, silent) { t = Math.max(0, Math.min(tMax, v)); range.value = t; fmtT(); if (!silent) opts.onTime(t); }
    function frame(now) {
      if (!playing) return;
      const dt = Math.min(0.08, (now - last) / 1000); last = now;
      let nt = t + dt * speed * Math.pow(10, +spd.value);
      if (nt >= tMax) nt = opts.loop === false ? tMax : nt % tMax;
      setT(nt);
      if (nt >= tMax && opts.loop === false) { pause(); return; }
      raf = requestAnimationFrame(frame);
    }
    function play() { if (playing) return; playing = true; playBtn.textContent = '❚❚ Pause'; last = performance.now(); raf = requestAnimationFrame(frame); }
    function pause() { playing = false; playBtn.textContent = '▶ Play'; cancelAnimationFrame(raf); }
    function toggle() { playing ? pause() : play(); }
    fmtT();
    return {
      el: box, get t() { return t; }, set: setT, play, pause, toggle,
      setMax(m) { tMax = m; range.max = m; range.step = m / 1000; speed = opts.speed || m / 10; if (t > m) setT(m); },
      destroy: pause,
    };
  }

  /** Stat tiles: returns {el, set(key, html)} */
  function stats(items) {
    const cells = {};
    const wrap = el('div', { class: 'stats' }, items.map((it) => { const v = el('div', { class: 'stat-v' }, '–'); cells[it.key] = v; return el('div', { class: 'stat' }, el('div', { class: 'stat-k', html: it.label }), v); }));
    tex(wrap);
    return { el: wrap, set(key, html) { if (cells[key]) cells[key].innerHTML = html; } };
  }

  /** redraw on resize + on theme change. returns cleanup fn */
  function onResize(elm, fn) {
    const ro = new ResizeObserver(() => fn());
    ro.observe(elm);
    return () => ro.disconnect();
  }

  /** small tabs widget: tabs([{label, render(container)}]) */
  function tabs(list, onSwitch) {
    const head = el('div', { class: 'tabs' }), body = el('div', { class: 'tab-body' });
    const btns = list.map((t, i) => el('button', { class: 'tab', onclick: () => sel(i), html: t.label }));
    head.append(...btns);
    function sel(i) { btns.forEach((b, j) => b.classList.toggle('active', i === j)); list.forEach((t, j) => t.el && (t.el.style.display = i === j ? '' : 'none')); onSwitch && onSwitch(i); }
    return { head, body, select: sel };
  }

  /** format a complex vector as ket expansion */
  function ket(v, labels, digits = 3) {
    const parts = [];
    for (let i = 0; i < v.n; i++) { const z = QM.V.get(v, i); if (QM.C.abs(z) < 1e-6) continue; parts.push(`(${QM.C.fmt(z, digits)})|${labels[i]}⟩`); }
    return parts.join(' + ') || '0';
  }

  QM.ui = { el, tex, layout, panel, canvas, params, button, buttons, amplitudes, transport, stats, onResize, tabs, ket };
})();
