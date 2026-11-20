/* Shared projective-measurement panel.
   QM.measurePanel({
     observables: () => [{key, label, op (Hermitian Matrix), fmtValue?(v)}],
     getState: () => Vector            (pure state)   – or –  getRho: () => Matrix (density matrix)
     onCollapse: (newStateOrRho, outcome) => void,
     basisLabels: [...]                 (for ket display of post-measurement states)
   }) -> {el, refresh() (observables changed), update() (state changed)} */
(function () {
  'use strict';
  const QM = (window.QM = window.QM || {});
  const { el } = QM.ui;

  function measurePanel(opts) {
    let obsList = [], current = null, groups = [], shots = null;
    const sel = el('select', { class: 'param-select' });
    sel.addEventListener('change', () => { pick(sel.value); });
    const table = el('table', { class: 'meas-table' });
    const result = el('div', { class: 'meas-result' });
    const hist = QM.ui.canvas(150);
    hist.style.display = 'none';
    const btnOnce = QM.ui.button('Measure &amp; collapse', measureOnce, 'primary');
    const btnShots = QM.ui.button('1000 shots (no collapse)', () => sampleShots(1000));
    const wrap = el('div', { class: 'measure' },
      el('div', { class: 'param' }, el('label', { class: 'param-label' }, 'Observable'), sel),
      el('div', { class: 'meas-scroll' }, table), QM.ui.buttons(btnOnce, btnShots), result, hist);

    function fmtVal(v) { return current && current.fmtValue ? current.fmtValue(v) : QM.fmt(v, 4); }
    function probs() {
      if (opts.getRho) { const r = opts.getRho(); return groups.map((g) => Math.max(0, QM.M.trace(QM.M.mul(g.P, r)).re)); }
      return QM.measureProbs(groups, opts.getState());
    }
    function refresh() {
      obsList = opts.observables();
      const prev = sel.value;
      sel.innerHTML = '';
      obsList.forEach((o) => { const opt = el('option', { value: o.key }); opt.innerHTML = o.label; sel.append(opt); });
      pick(obsList.some((o) => o.key === prev) ? prev : (obsList[0] && obsList[0].key));
    }
    function pick(key) {
      current = obsList.find((o) => o.key === key); if (!current) return;
      sel.value = key; groups = QM.spectral(current.op); shots = null; hist.style.display = 'none'; result.innerHTML = '';
      update();
    }
    function update() {
      if (!current) return;
      const p = probs();
      table.innerHTML = '';
      table.append(el('tr', null, el('th', null, 'value'), el('th', null, 'probability'), el('th', null, 'degeneracy')));
      const many = groups.length > 16;
      let hidden = 0;
      groups.forEach((g, i) => {
        if (many && p[i] < 1e-3) { hidden++; return; }
        table.append(el('tr', null,
          el('td', { class: 'mono' }, fmtVal(g.value)),
          el('td', null, el('div', { class: 'pbar' }, el('div', { class: 'pbar-fill', style: { width: (100 * p[i]).toFixed(1) + '%' } }), el('span', null, QM.fmt(p[i], 4)))),
          el('td', { class: 'mono' }, String(g.vectors.length))));
      });
      if (hidden) table.append(el('tr', null, el('td', { colspan: 3, class: 'note' }, `${hidden} further outcomes with p < 0.001 not listed (${groups.length} in total).`)));
      if (shots) drawHist(p);
    }
    function measureOnce() {
      const p = probs(), k = QM.sample(p), g = groups[k];
      if (opts.getRho) {
        const r = opts.getRho(), PrP = QM.M.mul(QM.M.mul(g.P, r), g.P), tr = QM.M.trace(PrP).re;
        opts.onCollapse(QM.M.scale(PrP, 1 / tr), g.value);
        result.innerHTML = `Outcome <b>${fmtVal(g.value)}</b> (p = ${QM.fmt(p[k], 4)}). State collapsed onto the eigenspace.`;
      } else {
        const psi = QM.collapse(g, opts.getState());
        opts.onCollapse(psi, g.value);
        result.innerHTML = `Outcome <b>${fmtVal(g.value)}</b> (p = ${QM.fmt(p[k], 4)}). New state: <span class="mono">${opts.basisLabels ? QM.ui.ket(psi, opts.basisLabels) : ''}</span>`;
      }
      update();
    }
    function sampleShots(N) {
      const p = probs(); const counts = new Array(groups.length).fill(0);
      for (let i = 0; i < N; i++) counts[QM.sample(p)]++;
      shots = { counts, N }; hist.style.display = ''; drawHist(p);
    }
    function drawHist(p) {
      const labels = groups.map((g) => fmtVal(g.value));
      const r = QM.plot.bars(hist, { labels, values: shots.counts.map((c) => c / shots.N), yMax: 1, colors: [QM.plot.theme().series[1] || '#e8590c'], ylabel: `relative frequency (${shots.N} shots) · ▬ Born rule` });
      const ctx = hist.getContext('2d'); ctx.strokeStyle = QM.plot.theme().ink; ctx.lineWidth = 2;
      p.forEach((pi, i) => { ctx.beginPath(); ctx.moveTo(r.X(i + 0.1), r.Y(pi)); ctx.lineTo(r.X(i + 0.9), r.Y(pi)); ctx.stroke(); });
    }
    return { el: wrap, refresh, update };
  }
  QM.measurePanel = measurePanel;
})();
