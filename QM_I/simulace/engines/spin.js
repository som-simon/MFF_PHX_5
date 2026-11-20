/* Engine 'spin' – single spin-½ lab (ħ = 1).
   State ρ = ½(I + P·σ) (pure ⇔ |P| = 1), measurement axis n, Hamiltonian H = ½ Ω·σ
   (Bloch vector precesses as dP/dt = Ω × P), ensemble averaging over flight times,
   Pauli-decomposition tool.

   Preset schema (all keys optional):
   {
     state:    {theta, phi, r?}  – initial Bloch vector P0 = r·(sinθcosφ, sinθsinφ, cosθ), r = 1 pure (default)
             | {bloch: [x, y, z]},
     axis:     {theta, phi}      – measurement axis n (default z)
     field:    {theta, phi, omega, sign?: ±1, fixed?: bool}
               – H = sign·(ω/2) n_B·σ, i.e. Ω = sign·ω·n_B; fixed hides the direction sliders
     dynamics: bool              – show the dynamics section (time t, precession charts)
     tMax:     number            – time-slider range (default 4π)
     t0:       number            – initial time
     ensemble: {on, dT, T?}      – ensemble over t ∈ [T, T+ΔT] (T = transport time; T sets t0)
     tools:    ['pauli', 'eigen']– extra panels: Pauli decomposition, eigenbasis of S_n
     pauli:    {matrix: [[[re,im],[re,im]],[[re,im],[re,im]]]} – initial matrix of the Pauli tool
     showSections: [...]         – override visible sections; keys:
                                   'state','info','measure','eigen','dynamics','ensemble','pauli'
   }
   Problems with sandbox:true show everything; otherwise a "Show all controls" checkbox reveals
   the hidden sections. */
(function () {
  'use strict';
  const QM = window.QM;
  const { el } = QM.ui;
  const { C, V, M } = QM;
  const PI = Math.PI, S2 = Math.SQRT1_2;

  // ---------- pure helpers (also used by the self-tests) ----------
  const dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross3 = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const norm3 = (a) => Math.hypot(a[0], a[1], a[2]);
  const dirOf = (th, ph) => [Math.sin(th) * Math.cos(ph), Math.sin(th) * Math.sin(ph), Math.cos(th)];
  function anglesOf(P) {
    const r = norm3(P); if (r < 1e-12) return { theta: 0, phi: 0, r: 0 };
    let ph = Math.atan2(P[1], P[0]); if (ph < 0) ph += 2 * PI;
    return { theta: Math.acos(Math.max(-1, Math.min(1, P[2] / r))), phi: Math.abs(P[0]) + Math.abs(P[1]) < 1e-12 ? 0 : ph, r: Math.min(1, r) };
  }
  // P = (P·k)k + c (P − (P·k)k) + s (k × P)
  function combine(P, k, c, s) { const pk = dot3(P, k), kx = cross3(k, P); return [0, 1, 2].map((i) => pk * k[i] + c * (P[i] - pk * k[i]) + s * kx[i]); }
  /** Bloch vector at time t under H = ½ Ω·σ (dP/dt = Ω × P) */
  function precess(P, Om, t) { const m = norm3(Om); if (m < 1e-14) return P.slice(); return combine(P, Om.map((x) => x / m), Math.cos(m * t), Math.sin(m * t)); }
  /** exact average of P(t) over t ∈ [T, T+dT] */
  function ensembleAvg(P, Om, T, dT) {
    const m = norm3(Om); if (m < 1e-14) return P.slice();
    if (dT * m < 1e-9) return precess(P, Om, T);
    const a = m * T, b = m * (T + dT);
    return combine(P, Om.map((x) => x / m), (Math.sin(b) - Math.sin(a)) / (b - a), (Math.cos(a) - Math.cos(b)) / (b - a));
  }
  function rhoOf(P) { return M.from([[(1 + P[2]) / 2, [P[0] / 2, -P[1] / 2]], [[P[0] / 2, P[1] / 2], (1 - P[2]) / 2]]); }
  const sig = [M.eye(2), QM.ops.sx(), QM.ops.sy(), QM.ops.sz()];
  const decompose = (Mx) => sig.map((s) => C.scale(M.trace(M.mul(s, Mx)), 0.5));
  const compose = (a) => M.add(...sig.map((s, i) => M.scale(s, a[i])));
  const Sn = (th, ph) => M.scale(QM.ops.sn(th, ph), 0.5);
  const Hof = (Om) => M.add(M.scale(sig[1], Om[0] / 2), M.scale(sig[2], Om[1] / 2), M.scale(sig[3], Om[2] / 2));
  const nPlus = (th, ph) => V.from([Math.cos(th / 2), C.polar(Math.sin(th / 2), ph)]);
  const nMinus = (th, ph) => V.from([C.polar(-Math.sin(th / 2), -ph), Math.cos(th / 2)]);
  /** eigenbasis of S_n; for n = x, y, z the textbook Pauli phase convention */
  function eigBasis(th, ph) {
    const n = dirOf(th, ph), near = (a) => Math.hypot(n[0] - a[0], n[1] - a[1], n[2] - a[2]) < 1e-6;
    if (near([1, 0, 0])) return { name: 'x', plus: V.from([S2, S2]), minus: V.from([S2, -S2]) };
    if (near([0, 1, 0])) return { name: 'y', plus: V.from([S2, [0, S2]]), minus: V.from([S2, [0, -S2]]) };
    if (near([0, 0, 1])) return { name: 'z', plus: V.from([1, 0]), minus: V.from([0, 1]) };
    return { name: 'n', plus: nPlus(th, ph), minus: nMinus(th, ph) };
  }
  /** state vector (cos θ/2, e^{iφ} sin θ/2) with the Bloch direction of P */
  function ketOf(P) { const a = anglesOf(P); return nPlus(a.theta, a.phi); }
  const csqrt = (z) => C.polar(Math.sqrt(C.abs(z)), C.arg(z) / 2);
  const vErr = (u, w) => { let e = 0; for (let i = 0; i < u.n; i++) e = Math.max(e, Math.hypot(u.re[i] - w.re[i], u.im[i] - w.im[i])); return e; };

  // ---------- small renderers ----------
  const bracket = { borderLeft: '1.5px solid currentColor', borderRight: '1.5px solid currentColor', display: 'inline-table', verticalAlign: 'middle' };
  function matEl(m, d = 3) {
    const rows = []; for (let i = 0; i < m.n; i++) { const tds = []; for (let j = 0; j < m.n; j++) tds.push(el('td', null, C.fmt(M.get(m, i, j), d))); rows.push(el('tr', null, tds)); }
    return el('table', { class: 'matrix', style: bracket }, rows);
  }
  const colEl = (v, d = 3) => el('table', { class: 'matrix', style: bracket }, [0, 1].map((i) => el('tr', null, el('td', null, C.fmt(V.get(v, i), d)))));
  const mrow = (...c) => el('div', { style: { display: 'flex', alignItems: 'center', gap: '6px', flexWrap: 'wrap', margin: '4px 0' } }, ...c);
  const f3 = (x) => QM.fmt(x, 3), fPi = (x) => QM.fmt(x / PI, 3) + 'π';
  const vec3 = (P) => '(' + P.map(f3).join(', ') + ')';
  const half = (x) => (x > 0 ? '+½' : '−½');
  const col = (i, d) => QM.plot.theme().series[i] || d;

  QM.registerEngine('spin', {
    title: 'Spin-½ lab',
    mount(host, preset, problem) {
      preset = preset || {};
      const sandbox = !!(problem && problem.sandbox);
      const tools = preset.tools || [];
      const s0 = preset.state || {};
      const a0 = s0.bloch ? anglesOf(s0.bloch) : { theta: s0.theta || 0, phi: s0.phi || 0, r: s0.r === undefined ? 1 : s0.r };
      const ax = preset.axis || {}, fd = preset.field || {}, en = preset.ensemble || {};
      const v = {
        theta: a0.theta, phi: a0.phi, r: a0.r, thn: ax.theta || 0, phn: ax.phi || 0,
        thB: fd.theta === undefined ? 0 : fd.theta, phB: fd.phi || 0, om: fd.omega === undefined ? 1 : fd.omega, sign: fd.sign === -1 ? -1 : 1,
        ens: !!en.on, dT: en.dT === undefined ? 20 : en.dT, showAll: sandbox,
      };
      const shown = new Set(preset.showSections || ['state', 'info', 'measure'].concat(
        tools.includes('eigen') ? ['eigen'] : [], preset.dynamics ? ['dynamics'] : [], preset.ensemble ? ['ensemble'] : [], tools.includes('pauli') ? ['pauli'] : []));
      const vis = (k) => (k === 'fieldDir' ? v.showAll || !fd.fixed : v.showAll || shown.has(k));
      const tMax = preset.tMax || 4 * PI;
      const sections = [];
      const sec = (key, p) => { sections.push([key, p.el]); return p; };
      const { controls, views } = QM.ui.layout(host);
      let cur = null, pauliA = null;

      // ----- show-all toggle -----
      if (!sandbox) {
        const p = QM.ui.panel(null);
        p.body.append(QM.ui.params([{ key: 'showAll', type: 'check', label: 'Show all controls (sandbox)' }], v, () => { applyVis(); meas.refresh(); update(); }).el);
        controls.append(p.el);
      }

      // ----- initial state -----
      const pState = sec('state', QM.ui.panel('Initial state $\\rho_0=\\tfrac12(I+\\vec P_0\\cdot\\vec\\sigma)$'));
      const stateParams = QM.ui.params([
        { key: 'theta', label: '$\\theta$ (polar angle)', type: 'range', min: 0, max: PI, step: PI / 180, fmt: fPi },
        { key: 'phi', label: '$\\varphi$ (azimuth)', type: 'range', min: 0, max: 2 * PI, step: PI / 180, fmt: fPi },
        { key: 'r', label: '$|\\vec P_0|$ (1 pure, &lt;1 mixed)', type: 'range', min: 0, max: 1, step: 0.01 },
      ], v, () => update());
      const quick = [['|z+⟩', 0, 0, 1], ['|z−⟩', PI, 0, 1], ['|x+⟩', PI / 2, 0, 1], ['|x−⟩', PI / 2, PI, 1], ['|y+⟩', PI / 2, PI / 2, 1], ['|y−⟩', PI / 2, 1.5 * PI, 1], ['I/2', 0, 0, 0]];
      const ketOut = el('div', { class: 'note mono' });
      pState.body.append(stateParams.el, QM.ui.buttons(...quick.map(([l, th, ph, r]) => QM.ui.button(l, () => setInitial(th, ph, r)))), ketOut);
      controls.append(pState.el);
      function setInitial(th, ph, r) { stateParams.set('theta', th, true); stateParams.set('phi', ph, true); stateParams.set('r', r, true); update(); }
      /** new initial state (after collapse / dephasing): P0 = P, t = 0, ensemble off */
      function setStateBloch(P) {
        const a = anglesOf(P); transport.set(0, true); ensParams.set('ens', false, true);
        setInitial(a.theta, a.phi, a.r);
      }

      // ----- measurement -----
      const pMeas = sec('measure', QM.ui.panel('Measurement axis $\\vec n$ &amp; Born rule'));
      const axisParams = QM.ui.params([
        { key: 'thn', label: '$\\theta_n$', type: 'range', min: 0, max: PI, step: PI / 180, fmt: fPi },
        { key: 'phn', label: '$\\varphi_n$', type: 'range', min: 0, max: 2 * PI, step: PI / 180, fmt: fPi },
      ], v, () => { meas.refresh(); update(); });
      const setAxis = (th, ph) => { axisParams.set('thn', th, true); axisParams.set('phn', ph, true); meas.refresh(); update(); };
      const probOut = el('div', { class: 'note mono' });
      const meas = QM.measurePanel({
        observables: () => [
          { key: 'sz', label: 'S_z', op: M.scale(sig[3], 0.5), fmtValue: half },
          { key: 'sx', label: 'S_x', op: M.scale(sig[1], 0.5), fmtValue: half },
          { key: 'sy', label: 'S_y', op: M.scale(sig[2], 0.5), fmtValue: half },
          { key: 'sn', label: 'S_n (axis n)', op: Sn(v.thn, v.phn), fmtValue: half },
        ],
        getRho: () => rhoOf(cur.P),
        onCollapse: (r) => setStateBloch(QM.rho.bloch(r)),
      });
      const dephase = (n) => { const d = dot3(cur.P, n); setStateBloch(n.map((x) => x * d)); };
      pMeas.body.append(axisParams.el,
        QM.ui.buttons(QM.ui.button('n = x', () => setAxis(PI / 2, 0)), QM.ui.button('n = y', () => setAxis(PI / 2, PI / 2)), QM.ui.button('n = z', () => setAxis(0, 0))),
        probOut, meas.el,
        el('div', { class: 'param-hint', html: 'Non-selective measurement (outcome not recorded): $\\rho\\to P_+\\rho P_+ + P_-\\rho P_-$' }),
        QM.ui.buttons(QM.ui.button('Non-selective $S_z$', () => dephase([0, 0, 1])), QM.ui.button('Non-selective $S_{\\vec n}$', () => dephase(cur.n))));
      controls.append(pMeas.el);

      // ----- dynamics -----
      const pDyn = sec('dynamics', QM.ui.panel('Dynamics $\\hat H=\\tfrac12\\vec\\Omega\\cdot\\vec\\sigma$'));
      const fieldDir = QM.ui.params([
        { key: 'thB', label: 'field $\\theta_B$', type: 'range', min: 0, max: PI, step: PI / 180, fmt: fPi },
        { key: 'phB', label: 'field $\\varphi_B$', type: 'range', min: 0, max: 2 * PI, step: PI / 180, fmt: fPi },
      ], v, () => update());
      sections.push(['fieldDir', fieldDir.el]);
      const fieldParams = QM.ui.params([
        { key: 'sign', label: 'Hamiltonian', type: 'select', options: [{ value: 1, label: 'H = +(ω/2) n_B·σ' }, { value: -1, label: 'H = −(ω/2) n_B·σ' }] },
        { key: 'om', label: '$\\omega$ (field strength)', type: 'range', min: 0, max: 4, step: 0.01 },
      ], v, () => update());
      const transport = QM.ui.transport({ tMax, t0: Math.min(tMax, preset.t0 || en.T || 0), onTime: () => update(), label: 'Time $t$ (ensemble: start $T$)' });
      const fieldOut = el('div', { class: 'note mono' });
      pDyn.body.append(fieldDir.el, fieldParams.el, fieldOut, transport.el);
      controls.append(pDyn.el);

      // ----- ensemble -----
      const pEns = sec('ensemble', QM.ui.panel('Ensemble / dephasing'));
      const ensParams = QM.ui.params([
        { key: 'ens', type: 'check', label: 'Average over flight times $t\\in[T,T+\\Delta T]$' },
        { key: 'dT', label: '$\\Delta T$', type: 'range', min: 0, max: 60, step: 0.1 },
      ], v, () => update());
      const ensOut = el('div', { class: 'note mono' });
      pEns.body.append(ensParams.el, ensOut);
      controls.append(pEns.el);

      // ----- views: Bloch + state info -----
      const pBloch = QM.ui.panel('Bloch sphere');
      const cvB = QM.ui.canvas(340);
      const bloch = QM.bloch(cvB);
      pBloch.body.append(cvB, el('div', { class: 'note' }, 'Drag to rotate. Blue: state P; green dashed: axis n; orange: field direction Ω̂; purple: Pauli tool â; dots: ensemble members.'));
      const pInfo = sec('info', QM.ui.panel('State at time $t$'));
      const stats = QM.ui.stats([
        { key: 'x', label: '$P_x=\\langle\\sigma_x\\rangle$' }, { key: 'y', label: '$P_y=\\langle\\sigma_y\\rangle$' }, { key: 'z', label: '$P_z=\\langle\\sigma_z\\rangle$' },
        { key: 'len', label: '$|\\vec P|$' }, { key: 'pur', label: '$\\mathrm{Tr}\\,\\rho^2$' },
        { key: 'pz', label: '$P(S_z=+\\frac12)$' }, { key: 'pn', label: '$P(S_{\\vec n}=+\\frac12)$' },
      ]);
      const rhoBox = el('div'), stateOut = el('div', { class: 'note mono' });
      pInfo.body.append(stats.el, rhoBox, stateOut);
      views.append(el('div', { class: 'row2' }, pBloch.el, pInfo.el));

      // ----- eigenbasis -----
      const pEig = sec('eigen', QM.ui.panel('Eigenbasis of $\\hat S_{\\vec n}=\\tfrac12\\vec n\\cdot\\vec\\sigma$'));
      const eigBox = el('div');
      pEig.body.append(eigBox, el('div', { class: 'note', html: 'General axis: $|n+\\rangle=(\\cos\\frac\\theta2,\\,e^{i\\varphi}\\sin\\frac\\theta2)$, $|n-\\rangle=(-e^{-i\\varphi}\\sin\\frac\\theta2,\\,\\cos\\frac\\theta2)$. For $\\vec n=x,y$ the Pauli convention $|x-\\rangle=(1,-1)/\\sqrt2$, $|y-\\rangle=(1,-i)/\\sqrt2$ is used (differs by a global phase).' }));
      views.append(pEig.el);

      // ----- charts -----
      const pCh = sec('dynamics', QM.ui.panel('Time evolution'));
      const cv1 = QM.ui.canvas(210), cv2 = QM.ui.canvas(190);
      pCh.body.append(cv1, cv2);
      views.append(pCh.el);

      // ----- Pauli tool -----
      const pPau = sec('pauli', QM.ui.panel('Pauli decomposition $M=a_0I+\\vec a\\cdot\\vec\\sigma$, $a_\\mu=\\tfrac12\\mathrm{Tr}(\\sigma_\\mu M)$'));
      const pm = (preset.pauli && preset.pauli.matrix) || [[[2, 0], [1, -1]], [[1, 1], [0, 0]]];
      const inputs = [];
      const grid = el('div', { class: 'amp-grid', style: { maxWidth: '320px' } }, el('span', { class: 'amp-h' }, ''), el('span', { class: 'amp-h' }, 'Re'), el('span', { class: 'amp-h' }, 'Im'));
      for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) {
        const re = el('input', { type: 'number', step: '0.1' }), im = el('input', { type: 'number', step: '0.1' });
        re.value = pm[i][j][0]; im.value = pm[i][j][1] || 0;
        [re, im].forEach((x) => x.addEventListener('input', pauliUpdate));
        inputs.push([re, im]); grid.append(el('span', { class: 'amp-lab' }, `M${i + 1}${j + 1}`), re, im);
      }
      const pPresets = [
        ['σx', () => sig[1]], ['σy', () => sig[2]], ['σz', () => sig[3]], ['I', () => sig[0]],
        ['Example [[2, 1−i], [1+i, 0]]', () => M.from([[2, [1, -1]], [[1, 1], 0]])],
        ['σ₊ = [[0, 1], [0, 0]] (non-Hermitian)', () => M.from([[0, 1], [0, 0]])],
        ['current S_n', () => Sn(v.thn, v.phn)], ['current ρ(t)', () => rhoOf(cur.P)],
      ];
      const pSel = el('select', { class: 'param-select', style: { maxWidth: '320px', marginBottom: '8px' } }, el('option', { value: '' }, 'Load matrix…'), pPresets.map((p, i) => el('option', { value: i }, p[0])));
      pSel.addEventListener('change', () => {
        if (pSel.value === '') return; const m = pPresets[+pSel.value][1]();
        for (let k = 0; k < 4; k++) { inputs[k][0].value = +m.re[k].toFixed(4); inputs[k][1].value = +m.im[k].toFixed(4); }
        pSel.value = ''; pauliUpdate();
      });
      const pauliOut = el('div', { class: 'mono', style: { margin: '8px 0' } });
      const names = ['I', 'σx', 'σy', 'σz'];
      const hs = el('table', { class: 'matrix' }, el('tr', null, el('td', null, '(A,B)'), names.map((n) => el('td', null, el('b', null, n)))),
        sig.map((A, i) => el('tr', null, el('td', null, el('b', null, names[i])), sig.map((B) => el('td', null, C.fmt(M.trace(M.mul(M.dagger(A), B))))))));
      pPau.body.append(el('div', { class: 'row2' },
        el('div', null, pSel, grid, pauliOut),
        el('div', null, el('div', { class: 'param-heading', html: 'Hilbert–Schmidt products $(A,B)=\\mathrm{Tr}(A^\\dagger B)$' }), hs,
          el('div', { class: 'note' }, 'Orthogonal, each with norm² = 2 ⇒ {I, σx, σy, σz} is an orthogonal basis of all complex 2×2 matrices.'))));
      views.append(pPau.el);

      function pauliUpdate() {
        const Mx = M.from([0, 1].map((i) => [0, 1].map((j) => { const [re, im] = inputs[i * 2 + j]; return [parseFloat(re.value) || 0, parseFloat(im.value) || 0]; })));
        const a = decompose(Mx), herm = M.isHermitian(Mx, 1e-9), err = M.maxAbs(M.sub(Mx, compose(a)));
        const aa = C.add(C.add(C.mul(a[1], a[1]), C.mul(a[2], a[2])), C.mul(a[3], a[3])), sq = csqrt(aa);
        const lines = a.map((z, i) => `a${i} = ${C.fmt(z, 4)}`);
        lines.push(`reconstruction error ‖M − (a0 I + a·σ)‖ = ${err.toExponential(1)}`);
        lines.push(`eigenvalues a0 ± √(a·a): ${C.fmt(C.add(a[0], sq), 4)}, ${C.fmt(C.sub(a[0], sq), 4)}`);
        pauliA = null;
        if (herm) {
          const av = [a[1].re, a[2].re, a[3].re], m = norm3(av);
          lines.push(`Hermitian: yes (all a_μ real). |a| = ${QM.fmt(m, 4)}, λ = a0 ± |a| = ${QM.fmt(a[0].re + m, 4)}, ${QM.fmt(a[0].re - m, 4)}`);
          lines.push(`numeric eigh: ${QM.eigh(Mx).values.map((x) => QM.fmt(x, 4)).join(', ')}`);
          if (m > 1e-12) { pauliA = av.map((x) => x / m); lines.push(`â = a/|a| = ${vec3(pauliA)} (eigenvector |â+⟩ ↔ λ = a0 + |a|)`); }
        } else lines.push('Hermitian: no (some a_μ not real) – no Bloch direction.');
        pauliOut.innerHTML = ''; lines.forEach((l) => pauliOut.append(el('div', null, l)));
        drawBloch();
      }

      // ----- logic -----
      function applyVis() { sections.forEach(([k, e]) => { e.style.display = vis(k) ? '' : 'none'; }); }
      function compute() {
        const P0 = dirOf(v.theta, v.phi).map((x) => x * v.r);
        const dyn = vis('dynamics'), Om = dirOf(v.thB, v.phB).map((x) => x * v.sign * v.om);
        const t = dyn ? transport.t : 0, ens = dyn && vis('ensemble') && v.ens;
        const P = ens ? ensembleAvg(P0, Om, t, v.dT) : dyn ? precess(P0, Om, t) : P0.slice();
        return { P0, Om, t, dyn, ens, P, n: dirOf(v.thn, v.phn) };
      }
      function drawBloch() {
        if (!cur) return;
        const { P, P0, Om, t } = cur, mag = norm3(Om), len = norm3(P);
        const vecs = [{ v: P, color: col(0, '#175ec5'), label: cur.ens ? 'ρ̄' : len > 0.999 ? 'ψ' : 'ρ' }];
        if (cur.dyn && mag > 1e-9) vecs.push({ v: Om.map((x) => x / mag), color: col(1, '#e8590c'), label: 'Ω' });
        if (vis('pauli') && pauliA) vecs.push({ v: pauliA, color: col(3, '#ae3ec9'), label: 'â' });
        let trail = null, points = null;
        if (cur.ens) { points = []; for (let k = 0; k < 80; k++) points.push(precess(P0, Om, t + (v.dT * (k + 0.5)) / 80)); }
        else if (cur.dyn && mag > 1e-9 && t > 0) { trail = []; const s0 = Math.max(0, t - (2 * PI) / mag); for (let k = 0; k <= 120; k++) trail.push(precess(P0, Om, s0 + ((t - s0) * k) / 120)); }
        bloch.draw({ vectors: vecs, axis: vis('measure') || vis('eigen') ? { v: cur.n, label: 'n' } : null, trail, points, caption: `|P| = ${f3(len)}` + (cur.dyn ? `   t = ${f3(t)}` : '') });
      }
      function drawCharts() {
        const { P0, Om, n, t } = cur, tEnd = cur.ens ? Math.max(tMax, t + v.dT) : tMax, N = 400;
        const xs = [], pz = [], pn = [], b = [[], [], []];
        for (let i = 0; i <= N; i++) { const s = (tEnd * i) / N, P = precess(P0, Om, s); xs.push(s); pz.push(0.5 * (1 + P[2])); pn.push(0.5 * (1 + dot3(n, P))); P.forEach((c, k) => b[k].push(c)); }
        const vlines = [{ x: t }], shade = cur.ens ? [{ from: t, to: t + v.dT }] : [];
        QM.plot.line(cv1, { xMin: 0, xMax: tEnd, yMin: 0, yMax: 1, xlabel: 't', ylabel: cur.ens ? 'probability (— — ensemble average)' : 'probability', vlines, shade,
          hlines: cur.ens ? [{ y: 0.5 * (1 + cur.P[2]) }, { y: 0.5 * (1 + dot3(n, cur.P)) }] : [],
          series: [{ x: xs, y: pz, label: 'P(S_z = +½)' }, { x: xs, y: pn, label: 'P(S_n = +½)', dash: [5, 4] }] });
        QM.plot.line(cv2, { xMin: 0, xMax: tEnd, yMin: -1.05, yMax: 1.05, xlabel: 't', ylabel: 'Bloch vector P(t)', vlines, shade,
          series: [{ x: xs, y: b[0], label: 'P_x' }, { x: xs, y: b[1], label: 'P_y' }, { x: xs, y: b[2], label: 'P_z' }] });
      }
      function fillEigen() {
        const eb = eigBasis(v.thn, v.phn), S = Sn(v.thn, v.phn), nm = eb.name, len = norm3(cur.P);
        const res = Math.max(vErr(M.apply(S, eb.plus), V.scale(eb.plus, C.of(0.5))), vErr(M.apply(S, eb.minus), V.scale(eb.minus, C.of(-0.5))));
        const Pp = M.outer(eb.plus, eb.plus), Pm = M.outer(eb.minus, eb.minus), sumErr = M.maxAbs(M.sub(M.add(Pp, Pm), M.eye(2)));
        let amp;
        if (len > 0.9999) {
          const psi = ketOf(cur.P), ap = V.dot(eb.plus, psi), am = V.dot(eb.minus, psi);
          amp = `⟨${nm}+|ψ⟩ = ${C.fmt(ap, 4)}  (|·|² = ${QM.fmt(C.abs2(ap), 4)}),   ⟨${nm}−|ψ⟩ = ${C.fmt(am, 4)}  (|·|² = ${QM.fmt(C.abs2(am), 4)})   [ψ = (cos θ/2, e^{iφ} sin θ/2)]`;
        } else {
          const r = rhoOf(cur.P);
          amp = `mixed state: ⟨${nm}+|ρ|${nm}+⟩ = ${QM.fmt(M.trace(M.mul(Pp, r)).re, 4)},  ⟨${nm}−|ρ|${nm}−⟩ = ${QM.fmt(M.trace(M.mul(Pm, r)).re, 4)}`;
        }
        eigBox.innerHTML = '';
        eigBox.append(
          mrow(`n = ${vec3(cur.n)};  S_${nm} =`, matEl(S), `eigenvalues (numeric): ${QM.eigh(S).values.map((x) => QM.fmt(x, 4)).join(', ')}`),
          mrow(`|${nm}+⟩ =`, colEl(eb.plus), `   |${nm}−⟩ =`, colEl(eb.minus), el('span', { class: 'note' }, `‖S|±⟩ ∓ ½|±⟩‖ = ${res.toExponential(1)}`)),
          mrow('P₊ = ½(I + n·σ) =', matEl(Pp), '  P₋ =', matEl(Pm), el('span', { class: 'note' }, `‖P₊ + P₋ − I‖ = ${sumErr.toExponential(1)}`)),
          el('div', { class: 'mono', style: { fontSize: '12.5px' } }, amp));
      }
      function update() {
        cur = compute();
        const { P, P0, Om, n } = cur, len = norm3(P), np = dot3(n, P);
        drawBloch();
        stats.set('x', f3(P[0])); stats.set('y', f3(P[1])); stats.set('z', f3(P[2]));
        stats.set('len', f3(len)); stats.set('pur', f3(0.5 * (1 + len * len)));
        stats.set('pz', QM.fmt(0.5 * (1 + P[2]), 4)); stats.set('pn', QM.fmt(0.5 * (1 + np), 4));
        rhoBox.innerHTML = ''; rhoBox.append(mrow('ρ =', matEl(rhoOf(P), 4)));
        stateOut.textContent = len > 0.9999 ? `pure: |ψ⟩ ≅ ${QM.ui.ket(ketOf(P), ['↑', '↓'])}  (up to a global phase)` : `mixed: |P| = ${f3(len)} < 1, no state vector`;
        ketOut.textContent = `P0 = ${vec3(P0)}` + (v.r > 0.9999 ? `,  |ψ0⟩ = ${QM.ui.ket(nPlus(v.theta, v.phi), ['↑', '↓'])}` : ' (mixed)');
        probOut.textContent = `n = ${vec3(n)},  n·P = ${f3(np)};  P(S_n = ±½) = ½(1 ± n·P) = ${QM.fmt(0.5 * (1 + np), 4)} / ${QM.fmt(0.5 * (1 - np), 4)}`;
        const mag = norm3(Om);
        fieldOut.textContent = `Ω = ${vec3(Om)},  |Ω| = ${f3(mag)}` + (mag > 1e-9 ? `,  period 2π/|Ω| = ${f3((2 * PI) / mag)}` : '') + '.  dP/dt = Ω × P';
        if (mag > 1e-9) { const k = Om.map((x) => x / mag), lim = k.map((x) => x * dot3(P0, k)); ensOut.textContent = `|Ω|ΔT = ${f3(mag * v.dT)};  average |P̄| = ${f3(len)}, Tr ρ̄² = ${f3(0.5 * (1 + len * len))}` + (cur.ens ? '' : ' (when on)') + `.  Limit |Ω|ΔT → ∞: P̄ → (P0·Ω̂)Ω̂ = ${vec3(lim)}`; }
        else ensOut.textContent = 'No field (ω = 0): nothing to average.';
        if (vis('eigen')) fillEigen();
        meas.update();
        if (vis('dynamics')) drawCharts();
      }
      function redraw() { if (!cur) return; drawBloch(); if (vis('dynamics')) drawCharts(); }

      applyVis();
      cur = compute();
      meas.refresh();
      pauliUpdate();
      update();
      const offResize = QM.ui.onResize(views, redraw);
      return { destroy() { offResize(); transport.destroy(); }, redraw };
    },
  });

  // ---------- self-tests ----------
  const ok = (pass, detail) => ({ pass, detail });
  QM.addTest('Spin: Hilbert–Schmidt products Tr(σμ†σν) = 2δμν', () => {
    let e = 0; sig.forEach((A, i) => sig.forEach((B, j) => { const z = M.trace(M.mul(M.dagger(A), B)); e = Math.max(e, Math.hypot(z.re - (i === j ? 2 : 0), z.im)); }));
    return ok(e < 1e-12, 'max error ' + e.toExponential(1));
  });
  QM.addTest('Spin: Pauli decomposition a_μ = ½Tr(σμ M) reconstructs complex M; example a = (1,1,1,1)', () => {
    const Mx = M.from([[[0.3, -1.2], [2.1, 0.4]], [[-0.7, 0.9], [1.5, -0.2]]]);
    const e1 = M.maxAbs(M.sub(Mx, compose(decompose(Mx))));
    const a = decompose(M.from([[2, [1, -1]], [[1, 1], 0]]));
    const e2 = Math.max(...a.map((z) => Math.hypot(z.re - 1, z.im)));
    return ok(e1 < 1e-12 && e2 < 1e-12, `recon ${e1.toExponential(1)}, example ${e2.toExponential(1)}`);
  });
  QM.addTest('Spin: σx,σy,σz have eigenvalues ±1, eigenvectors |x±⟩,|y±⟩,|z±⟩, P₊+P₋ = I', () => {
    let e = 0; const axes = [[PI / 2, 0], [PI / 2, PI / 2], [0, 0]];
    axes.forEach(([th, ph], i) => {
      const s = sig[i + 1], eb = eigBasis(th, ph), vals = QM.eigh(s).values;
      e = Math.max(e, Math.abs(vals[0] + 1), Math.abs(vals[1] - 1), vErr(M.apply(s, eb.plus), eb.plus), vErr(M.apply(s, eb.minus), V.scale(eb.minus, C.of(-1))));
      const Pp = M.scale(M.add(M.eye(2), s), 0.5), Pm = M.scale(M.sub(M.eye(2), s), 0.5);
      e = Math.max(e, M.maxAbs(M.sub(M.add(Pp, Pm), M.eye(2))), M.maxAbs(M.sub(Pp, M.outer(eb.plus, eb.plus))), M.maxAbs(M.sub(Pm, M.outer(eb.minus, eb.minus))));
    });
    return ok(e < 1e-9, 'max error ' + e.toExponential(1));
  });
  QM.addTest('Spin: ⟨x+|y+⟩ = (1+i)/2, ⟨x−|y+⟩ = (1−i)/2', () => {
    const x = eigBasis(PI / 2, 0), y = eigBasis(PI / 2, PI / 2), a = V.dot(x.plus, y.plus), b = V.dot(x.minus, y.plus);
    const e = Math.max(Math.hypot(a.re - 0.5, a.im - 0.5), Math.hypot(b.re - 0.5, b.im + 0.5));
    return ok(e < 1e-12, `⟨x+|y+⟩ = ${C.fmt(a)}, ⟨x−|y+⟩ = ${C.fmt(b)}`);
  });
  QM.addTest('Spin: S_n has eigenvalues ±½ with analytic |n±⟩ (random axes)', () => {
    let e = 0;
    [[0.3, 1.1], [1.2, 4.0], [2.7, 2.2], [PI / 3, PI / 4]].forEach(([th, ph]) => {
      const S = Sn(th, ph), vals = QM.eigh(S).values, p = nPlus(th, ph), m = nMinus(th, ph);
      e = Math.max(e, Math.abs(vals[0] + 0.5), Math.abs(vals[1] - 0.5), vErr(M.apply(S, p), V.scale(p, C.of(0.5))), vErr(M.apply(S, m), V.scale(m, C.of(-0.5))), C.abs(V.dot(p, m)));
    });
    return ok(e < 1e-9, 'max error ' + e.toExponential(1));
  });
  QM.addTest('Spin: Born rule P(+) = ½(1+n·m); |z+⟩ → cos²(θ/2), |x+⟩ → ½(1+sinθ cosφ)', () => {
    let e = 0; const th = 1.1, ph = 0.8, n = dirOf(th, ph), p = nPlus(th, ph);
    e = Math.max(e, Math.abs(C.abs2(V.dot(p, V.from([1, 0]))) - Math.cos(th / 2) ** 2));
    e = Math.max(e, Math.abs(C.abs2(V.dot(p, V.from([S2, S2]))) - 0.5 * (1 + Math.sin(th) * Math.cos(ph))));
    const mth = 2.3, mph = 5.1, psi = nPlus(mth, mph);
    e = Math.max(e, Math.abs(C.abs2(V.dot(p, psi)) - 0.5 * (1 + dot3(n, dirOf(mth, mph)))));
    const g = QM.spectral(Sn(th, ph)), pr = QM.measureProbs(g, psi);
    e = Math.max(e, Math.abs(pr[1] - 0.5 * (1 + dot3(n, dirOf(mth, mph)))));
    return ok(e < 1e-9, 'max error ' + e.toExponential(1));
  });
  QM.addTest('Spin: precession dP/dt = Ω×P matches exact evolution e^{−iHt}, H = ½Ω·σ', () => {
    const Om = [0.7, -1.3, 0.4], psi0 = nPlus(1.1, 0.7), P0 = QM.blochOf(psi0), ev = QM.evolver(QM.eigh(Hof(Om)), psi0);
    let e = 0; [0.3, 1.7, 5.2, 11].forEach((t) => { const a = QM.blochOf(ev(t)), b = precess(P0, Om, t); e = Math.max(e, ...a.map((x, i) => Math.abs(x - b[i]))); });
    return ok(e < 1e-9, 'max error ' + e.toExponential(1));
  });
  function precPreset() {
    const p = QM.registry.problems.find((q) => q.id === 'spin-precession');
    const pr = p.preset, th = pr.state.theta, w = pr.field.omega;
    return { th, w, Om: dirOf(pr.field.theta, pr.field.phi).map((x) => x * (pr.field.sign || 1) * w) };
  }
  QM.addTest('Spin precession (DÚ 1): H = −(ω/2)σx ⇒ P(T) = (sinθ, cosθ sinωT, cosθ cosωT), P(S_z=+½) = ½(1+cosθ cosωT)', () => {
    const { th, w, Om } = precPreset();
    const H = M.scale(sig[1], -w / 2), eH = M.maxAbs(M.sub(H, Hof(Om)));
    const psi0 = nPlus(th, 0), ev = QM.evolver(QM.eigh(H), psi0);
    let e = 0;
    [0.4, 1.3, 2.9, 7.7].forEach((T) => {
      const ana = [Math.sin(th), Math.cos(th) * Math.sin(w * T), Math.cos(th) * Math.cos(w * T)];
      const P = precess(QM.blochOf(psi0), Om, T), psi = ev(T);
      const U = M.add(M.scale(M.eye(2), Math.cos(w * T / 2)), M.scale(sig[1], C.of(0, Math.sin(w * T / 2))));
      e = Math.max(e, ...P.map((x, i) => Math.abs(x - ana[i])), ...QM.blochOf(psi).map((x, i) => Math.abs(x - ana[i])), vErr(psi, M.apply(U, psi0)),
        Math.abs(C.abs2(V.get(psi, 0)) - 0.5 * (1 + Math.cos(th) * Math.cos(w * T))));
    });
    return ok(e < 1e-9 && eH < 1e-12, `max error ${e.toExponential(1)} (θ = ${QM.fmt(th, 4)}, ω = ${w})`);
  });
  QM.addTest('Spin precession: ensemble over [T, T+ΔT], ωΔT ≫ 1 ⇒ ρ = ½(I + sinθ σx); analytic = sampled average', () => {
    const { th, w, Om } = precPreset(), P0 = dirOf(th, 0);
    const A = ensembleAvg(P0, Om, 0.37, (2 * PI * 50) / w), e1 = Math.max(Math.abs(A[0] - Math.sin(th)), Math.abs(A[1]), Math.abs(A[2]));
    const B = ensembleAvg(P0, Om, 1.1, 137.1), e2 = Math.max(Math.abs(B[0] - Math.sin(th)), Math.abs(B[1]), Math.abs(B[2]));
    const N = 20000, s = [0, 0, 0], T = 0.8, dT = 3.3;
    for (let k = 0; k < N; k++) precess(P0, Om, T + (dT * (k + 0.5)) / N).forEach((x, i) => (s[i] += x / N));
    const C3 = ensembleAvg(P0, Om, T, dT), e3 = Math.max(...s.map((x, i) => Math.abs(x - C3[i])));
    const pur = QM.rho.purity(rhoOf(A)), e4 = Math.abs(pur - 0.5 * (1 + Math.sin(th) ** 2));
    return ok(e1 < 1e-9 && e2 < 2 / (w * 137.1) && e3 < 1e-6 && e4 < 1e-9 && pur < 1, `full periods ${e1.toExponential(1)}, ωΔT=137: ${e2.toExponential(1)}, sampled ${e3.toExponential(1)}, Tr ρ² = ${QM.fmt(pur, 4)}`);
  });
  QM.addTest('Spin: non-selective S_z kills off-diagonals; Tr(P±ρ) = ½(1 ± n·P) for mixed ρ', () => {
    const P = [0.3, -0.4, 0.5], r = rhoOf(P);
    const g = QM.spectral(M.scale(sig[3], 0.5)), rr = g.reduce((acc, x) => M.add(acc, M.mul(M.mul(x.P, r), x.P)), M.zeros(2));
    let e = Math.max(M.maxAbs(M.sub(rr, rhoOf([0, 0, P[2]]))), Math.hypot(rr.re[1], rr.im[1]));
    const n = dirOf(1.9, 2.6), gn = QM.spectral(Sn(1.9, 2.6));
    e = Math.max(e, Math.abs(M.trace(M.mul(gn[1].P, r)).re - 0.5 * (1 + dot3(n, P))), Math.abs(M.trace(M.mul(gn[0].P, r)).re - 0.5 * (1 - dot3(n, P))));
    const b = QM.rho.bloch(r); e = Math.max(e, ...b.map((x, i) => Math.abs(x - P[i])));
    return ok(e < 1e-9, 'max error ' + e.toExponential(1));
  });
})();
