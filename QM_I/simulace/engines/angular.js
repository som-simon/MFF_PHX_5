/* Engine 'angular' – decomposition of an angular function f(θ,φ) into spherical harmonics Y_lm (l ≤ 6),
   measurement of L² and L_z (Born rule + collapse).  ħ = 1.

   f is a polynomial in the direction cosines x = n_x = sinθcosφ, y = n_y = sinθsinφ, z = n_z = cosθ,
   typed as an expression: numbers, x y z (or nx ny nz), i, pi, the parameters c and k, + − * / ^ ( ),
   implicit multiplication (e.g. "2xy", "(x+i*y)^2 z"). Exponents must evaluate to integers 0…12.

   Preset schema (all optional):
   { expr: '(x + c*y)^k',   angular function (string as above)
     c: 1, k: 2,            values of the parameters c, k used in expr
     projL: 0,              l of the projection shown in the right-hand map
     radial: 'HTML'         description of the radial part R(r) (only displayed; L², L_z do not act on it) }

   c_lm = ∫ Y_lm* f dΩ by Gauss–Legendre (32 nodes in cosθ) × uniform φ (64) quadrature (exact for these
   polynomials); Y_lm with Condon–Shortley phase, Y_{l,−m} = (−1)^m Y_lm*. */
(function () {
  'use strict';
  const QM = window.QM;
  const { el } = QM.ui;
  const LMAX = 6, NU = 32, NPHI = 64, HX = 128, HY = 64, NLM = (LMAX + 1) * (LMAX + 1);
  const lmIdx = (l, m) => l * l + l + m; // 0 … NLM−1
  const pIdx = (l, m) => (l * (l + 1)) / 2 + m; // m ≥ 0

  // ---------- polynomials in (x, y, z) with complex coefficients: {"a,b,c": [re, im]} ----------
  const Poly = {
    c(re, im = 0) { return { '0,0,0': [re, im] }; },
    v(a, b, c) { return { [a + ',' + b + ',' + c]: [1, 0] }; },
    add(p, q, s = 1) { const r = {}; for (const k in p) r[k] = p[k].slice(); for (const k in q) { const t = r[k] || [0, 0]; r[k] = [t[0] + s * q[k][0], t[1] + s * q[k][1]]; } return r; },
    mul(p, q) {
      const r = {};
      for (const a in p) for (const b in q) {
        const ea = a.split(',').map(Number), eb = b.split(',').map(Number), e = [ea[0] + eb[0], ea[1] + eb[1], ea[2] + eb[2]];
        if (e[0] + e[1] + e[2] > 12) throw new Error('degree > 12 not supported');
        const k = e.join(','), t = r[k] || [0, 0], x = p[a], y = q[b];
        r[k] = [t[0] + x[0] * y[0] - x[1] * y[1], t[1] + x[0] * y[1] + x[1] * y[0]];
      }
      return r;
    },
    constant(p) { let v = [0, 0]; for (const k in p) { if (k !== '0,0,0' && Math.hypot(p[k][0], p[k][1]) > 1e-14) return null; if (k === '0,0,0') v = p[k]; } return v; },
    terms(p) { const t = []; for (const k in p) { const [a, b, c] = k.split(',').map(Number); if (Math.hypot(p[k][0], p[k][1]) > 1e-14) t.push({ a, b, c, re: p[k][0], im: p[k][1] }); } return t; },
  };

  /** parse an expression into polynomial terms; vars = {c, k} */
  function parse(src, vars) {
    const s = String(src).replace(/\s+/g, '').toLowerCase();
    const toks = s.match(/n[xyz]|pi|\d*\.?\d+(?:e[-+]?\d+)?|[a-z]|[-+*/^()]/g) || [];
    if (toks.join('') !== s) throw new Error('unexpected character');
    if (!toks.length) throw new Error('empty expression');
    let i = 0;
    const peek = () => toks[i], eat = () => toks[i++];
    function expr() {
      let p = term();
      while (peek() === '+' || peek() === '-') { const op = eat(); p = Poly.add(p, term(), op === '-' ? -1 : 1); }
      return p;
    }
    function term() {
      let p = unary();
      for (;;) {
        const t = peek();
        if (t === '*') { eat(); p = Poly.mul(p, unary()); }
        else if (t === '/') {
          eat(); const d = Poly.constant(unary());
          if (!d || Math.hypot(d[0], d[1]) < 1e-14) throw new Error('can only divide by a non-zero number');
          const n2 = d[0] * d[0] + d[1] * d[1]; p = Poly.mul(p, Poly.c(d[0] / n2, -d[1] / n2));
        } else if (t !== undefined && t !== ')' && t !== '+' && t !== '-' && t !== '^') p = Poly.mul(p, power()); // implicit ×
        else return p;
      }
    }
    function unary() { if (peek() === '-') { eat(); return Poly.mul(Poly.c(-1), unary()); } if (peek() === '+') { eat(); return unary(); } return power(); }
    function power() {
      const b = atom();
      if (peek() !== '^') return b;
      eat(); const e = Poly.constant(unary());
      if (!e || Math.abs(e[1]) > 1e-12 || Math.abs(e[0] - Math.round(e[0])) > 1e-9 || e[0] < 0 || e[0] > 12) throw new Error('exponent must be an integer 0…12');
      let r = Poly.c(1); for (let n = 0; n < Math.round(e[0]); n++) r = Poly.mul(r, b);
      return r;
    }
    function atom() {
      const t = eat();
      if (t === undefined) throw new Error('unexpected end');
      if (t === '(') { const p = expr(); if (eat() !== ')') throw new Error('missing )'); return p; }
      if (/^[\d.]/.test(t)) return Poly.c(parseFloat(t));
      if (t === 'x' || t === 'nx') return Poly.v(1, 0, 0);
      if (t === 'y' || t === 'ny') return Poly.v(0, 1, 0);
      if (t === 'z' || t === 'nz') return Poly.v(0, 0, 1);
      if (t === 'i') return Poly.c(0, 1);
      if (t === 'pi') return Poly.c(Math.PI);
      if (t === 'c' || t === 'k') return Poly.c(+vars[t] || 0);
      throw new Error('unknown symbol "' + t + '"');
    }
    const p = expr();
    if (i < toks.length) throw new Error('unexpected "' + toks[i] + '"');
    return Poly.terms(p);
  }

  function evalPoly(terms, nx, ny, nz) {
    let re = 0, im = 0;
    for (const t of terms) { const v = Math.pow(nx, t.a) * Math.pow(ny, t.b) * Math.pow(nz, t.c); re += t.re * v; im += t.im * v; }
    return [re, im];
  }

  // ---------- spherical harmonics ----------
  const fact = [1]; for (let n = 1; n <= 2 * LMAX + 2; n++) fact[n] = fact[n - 1] * n;
  const K = []; for (let l = 0; l <= LMAX; l++) for (let m = 0; m <= l; m++) K[pIdx(l, m)] = Math.sqrt(((2 * l + 1) / (4 * Math.PI)) * (fact[l - m] / fact[l + m]));
  /** normalised K_lm P_l^m(u) (Condon–Shortley phase) for 0 ≤ m ≤ l ≤ LMAX */
  function legendre(u) {
    const s = Math.sqrt(Math.max(0, 1 - u * u)), out = new Float64Array(pIdx(LMAX, LMAX) + 1);
    let pmm = 1;
    for (let m = 0; m <= LMAX; m++) {
      if (m > 0) pmm *= -(2 * m - 1) * s;
      out[pIdx(m, m)] = K[pIdx(m, m)] * pmm;
      if (m === LMAX) break;
      let p2 = pmm, p1 = u * (2 * m + 1) * pmm;
      out[pIdx(m + 1, m)] = K[pIdx(m + 1, m)] * p1;
      for (let l = m + 2; l <= LMAX; l++) { const p = ((2 * l - 1) * u * p1 - (l + m - 1) * p2) / (l - m); out[pIdx(l, m)] = K[pIdx(l, m)] * p; p2 = p1; p1 = p; }
    }
    return out;
  }
  /** Y_lm(θ,φ) = sgn · K P_l^{|m|}(cosθ) e^{imφ},  sgn = (−1)^m for m < 0 */
  const ysgn = (m) => (m < 0 && m % 2 ? -1 : 1);
  function Ylm(l, m, u, phi) { const v = ysgn(m) * legendre(u)[pIdx(l, Math.abs(m))]; return [v * Math.cos(m * phi), v * Math.sin(m * phi)]; }

  function gaussLegendre(n) {
    const x = new Float64Array(n), w = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let z = Math.cos((Math.PI * (i + 0.75)) / (n + 0.5)), pp = 1;
      for (let it = 0; it < 100; it++) {
        let p1 = 1, p2 = 0;
        for (let j = 1; j <= n; j++) { const p3 = p2; p2 = p1; p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j; }
        pp = (n * (z * p1 - p2)) / (z * z - 1);
        const z1 = z; z = z1 - p1 / pp;
        if (Math.abs(z - z1) < 1e-15) break;
      }
      x[i] = z; w[i] = 2 / ((1 - z * z) * pp * pp);
    }
    return { x, w };
  }
  const GL = gaussLegendre(NU), GLP = Array.from(GL.x, legendre);

  /** generic quadrature of a complex field F(u, φ) against all Y_lm*: returns {re, im, norm (∫|F|²)} */
  function project(F) {
    const cre = new Float64Array(NLM), cim = new Float64Array(NLM), dphi = (2 * Math.PI) / NPHI;
    let norm = 0;
    for (let iu = 0; iu < NU; iu++) {
      const u = GL.x[iu], s = Math.sqrt(1 - u * u), w = GL.w[iu];
      const fr = new Float64Array(NPHI), fi = new Float64Array(NPHI);
      for (let ip = 0; ip < NPHI; ip++) { const ph = (ip + 0.5) * dphi, v = F(u, ph, s); fr[ip] = v[0]; fi[ip] = v[1]; norm += w * dphi * (v[0] * v[0] + v[1] * v[1]); }
      for (let m = -LMAX; m <= LMAX; m++) { // A_m = ∫ f e^{−imφ} dφ
        let ar = 0, ai = 0;
        for (let ip = 0; ip < NPHI; ip++) { const ph = (ip + 0.5) * dphi, c = Math.cos(m * ph), sn = Math.sin(m * ph); ar += fr[ip] * c + fi[ip] * sn; ai += fi[ip] * c - fr[ip] * sn; }
        ar *= dphi; ai *= dphi;
        for (let l = Math.abs(m); l <= LMAX; l++) { const y = w * ysgn(m) * GLP[iu][pIdx(l, Math.abs(m))]; cre[lmIdx(l, m)] += y * ar; cim[lmIdx(l, m)] += y * ai; }
      }
    }
    return { re: cre, im: cim, norm };
  }
  function decompose(terms) { return project((u, ph, s) => evalPoly(terms, s * Math.cos(ph), s * Math.sin(ph), u)); }

  function probsL(c) { const P = new Float64Array(LMAX + 1); for (let l = 0; l <= LMAX; l++) for (let m = -l; m <= l; m++) { const j = lmIdx(l, m); P[l] += c.re[j] * c.re[j] + c.im[j] * c.im[j]; } return P; }

  // heat-map grid: u = cosθ rows (bottom → top), φ columns
  const HU = Array.from({ length: HY }, (_, j) => -1 + ((j + 0.5) * 2) / HY), HLP = HU.map(legendre);
  function fieldMap(c, onlyL) { // |Σ c_lm Y_lm|² on the map grid
    const out = new Float32Array(HX * HY);
    for (let j = 0; j < HY; j++) {
      const A = []; // A_m(u) = Σ_l c_lm Y_lm(u) without e^{imφ}
      for (let m = -LMAX; m <= LMAX; m++) {
        let ar = 0, ai = 0;
        for (let l = Math.abs(m); l <= LMAX; l++) { if (onlyL !== undefined && l !== onlyL) continue; const y = ysgn(m) * HLP[j][pIdx(l, Math.abs(m))], k = lmIdx(l, m); ar += y * c.re[k]; ai += y * c.im[k]; }
        A.push([m, ar, ai]);
      }
      for (let i = 0; i < HX; i++) {
        const ph = ((i + 0.5) * 2 * Math.PI) / HX; let re = 0, im = 0;
        for (const [m, ar, ai] of A) { if (!ar && !ai) continue; const cc = Math.cos(m * ph), ss = Math.sin(m * ph); re += ar * cc - ai * ss; im += ar * ss + ai * cc; }
        out[j * HX + i] = re * re + im * im;
      }
    }
    return out;
  }

  // ---------- engine ----------
  const TEMPLATES = [
    ['(x + c*y)^k', '(n_x + c·n_y)^k'], ['(x + i*c*y)^k', '(n_x + i c·n_y)^k'], ['z^k', 'n_z^k'], ['x*y', 'n_x n_y'], ['x*y*z', 'n_x n_y n_z'],
    ['x^2 - y^2', 'n_x² − n_y²'], ['3z^2 - 1', '3n_z² − 1  (∝ Y20)'], ['1 + c*z', '1 + c·n_z'], ['(x + i*y)^2 z', '(n_x + i n_y)² n_z  (∝ Y32)'], ['x^4 + y^4 + z^4', 'n_x⁴ + n_y⁴ + n_z⁴ (cubic)'],
  ];

  QM.registerEngine('angular', {
    title: 'Spherical-harmonic decomposition',
    mount(host, preset) {
      const P = Object.assign({ expr: '(x + c*y)^k', c: 1, k: 2, projL: 0, radial: '' }, preset || {});
      const { controls, views } = QM.ui.layout(host);
      let terms = [], raw = null, state = null, rest = 0, collapsedMsg = '';

      // --- controls ---
      const fPanel = QM.ui.panel('Angular function $f(\\theta,\\varphi)$');
      const tsel = el('select', { class: 'param-select' }, el('option', { value: '' }, 'Template…'), TEMPLATES.map(([e, l]) => el('option', { value: e }, l)));
      const inp = el('input', { type: 'text', class: 'param-select mono', value: P.expr, spellcheck: 'false' });
      const err = el('div', { class: 'amp-msg' });
      tsel.addEventListener('change', () => { if (!tsel.value) return; inp.value = P.expr = tsel.value; tsel.value = ''; rebuild(); });
      inp.addEventListener('change', () => { P.expr = inp.value; rebuild(); });
      const prm = QM.ui.params([
        { key: 'c', label: 'parameter $c$', type: 'range', min: -3, max: 3, step: 0.05 },
        { key: 'k', label: 'parameter $k$ (integer)', type: 'range', min: 0, max: 6, step: 1 },
      ], P, () => rebuild());
      fPanel.body.append(tsel, el('div', { class: 'param', style: { marginTop: '8px' } }, el('label', { class: 'param-label', html: 'expression in $x=n_x,\\ y=n_y,\\ z=n_z$' }), inp), err, prm.el,
        el('div', { class: 'note', html: 'On the unit sphere $x^2+y^2+z^2=1$. Example: <span class="mono">(x+y)^2</span>, <span class="mono">2xy + i z</span>.' }));

      const mPanel = QM.ui.panel('Measurement');
      const mOut = el('div', { class: 'meas-result' });
      const lSel = QM.ui.params([{ key: 'projL', label: 'show projection onto $l$', type: 'select', options: Array.from({ length: LMAX + 1 }, (_, l) => ({ value: l, label: 'l = ' + l })) }], P, () => redraw());
      mPanel.body.append(lSel.el,
        QM.ui.buttons(QM.ui.button('Measure $L^2$', () => measure('L2'), 'primary'), QM.ui.button('Measure $L_z$', () => measure('Lz')), QM.ui.button('⟲ Back to f', () => rebuild())), mOut);
      const rPanel = QM.ui.panel('Full wavefunction');
      const rOut = el('div', { class: 'statement' });
      rPanel.body.append(rOut);
      controls.append(fPanel.el, mPanel.el, rPanel.el);

      // --- views ---
      const cA = QM.ui.canvas(250), cB = QM.ui.canvas(250), cL = QM.ui.canvas(190);
      const pA = QM.ui.panel('Current state $|\\psi(\\theta,\\varphi)|^2$'), pB = QM.ui.panel('Projection onto $l$');
      pA.body.append(cA); pB.body.append(cB);
      const st = QM.ui.stats([{ key: 'L2', label: '$\\langle L^2\\rangle$' }, { key: 'Lz', label: '$\\langle L_z\\rangle$' }, { key: 'dLz', label: '$\\Delta L_z$' }, { key: 'nrm', label: '$\\int|f|^2d\\Omega$' }, { key: 'rest', label: 'weight $l>6$' }]);
      const pL = QM.ui.panel('$P(l)$ – outcomes of $L^2 = l(l+1)$');
      pL.body.append(cL);
      const pT = QM.ui.panel('$P(l,m)=|c_{lm}|^2$ and amplitudes');
      const tbl = el('div', { class: 'matrix-wrap' }), amps = el('div', { class: 'note mono' });
      pT.body.append(tbl, amps);
      const sP = QM.ui.panel('Summary'); sP.body.append(st.el);
      views.append(el('div', { class: 'row2' }, pA.el, pB.el), sP.el, el('div', { class: 'row2' }, pL.el, pT.el));
      views.append(el('div', { class: 'note', html: 'Maps: equal-area projection – horizontal $\\varphi\\in[0,2\\pi)$, vertical $\\cos\\theta$ (top = north pole). Equal areas carry equal solid angle.' }));

      function rebuild() {
        collapsedMsg = ''; mOut.innerHTML = '';
        try {
          terms = parse(P.expr, P); raw = decompose(terms);
          if (raw.norm < 1e-12) throw new Error('f vanishes on the sphere');
          err.textContent = '';
        } catch (e) { err.textContent = 'Error: ' + e.message; return; }
        let cap = 0; for (let j = 0; j < NLM; j++) cap += raw.re[j] * raw.re[j] + raw.im[j] * raw.im[j];
        rest = Math.max(0, 1 - cap / raw.norm);
        const s = Math.sqrt(cap);
        state = { re: raw.re.map((v) => v / s), im: raw.im.map((v) => v / s) };
        redraw();
      }

      function measure(kind) {
        if (!state) return;
        const pl = probsL(state);
        let probs, label;
        if (kind === 'L2') probs = Array.from(pl);
        else { probs = new Array(2 * LMAX + 1).fill(0); for (let l = 0; l <= LMAX; l++) for (let m = -l; m <= l; m++) { const j = lmIdx(l, m); probs[m + LMAX] += state.re[j] ** 2 + state.im[j] ** 2; } }
        const k = QM.sample(probs), pk = probs[k];
        for (let l = 0; l <= LMAX; l++) for (let m = -l; m <= l; m++) {
          const keep = kind === 'L2' ? l === k : m === k - LMAX, j = lmIdx(l, m);
          if (keep) { state.re[j] /= Math.sqrt(pk); state.im[j] /= Math.sqrt(pk); } else { state.re[j] = 0; state.im[j] = 0; }
        }
        if (kind === 'L2') { label = `l = ${k}, L² = ${k * (k + 1)}ħ²`; P.projL = k; lSel.set('projL', k, true); }
        else label = `m = ${k - LMAX}, L<sub>z</sub> = ${k - LMAX}ħ`;
        collapsedMsg = `Outcome <b>${label}</b> (p = ${QM.fmt(pk, 4)}). The angular part collapsed onto the eigenspace${rest > 1e-9 ? ' (components with l &gt; 6 ignored)' : ''}.`;
        redraw();
      }

      function ampText() {
        const parts = [];
        for (let l = 0; l <= LMAX; l++) for (let m = -l; m <= l; m++) { const j = lmIdx(l, m), z = { re: state.re[j], im: state.im[j] }; if (QM.C.abs2(z) > 1e-10) parts.push(`(${QM.C.fmt(z, 4)}) Y${l},${m}`); }
        return parts.join(' + ') || '0';
      }

      function redraw() {
        if (!state) return;
        const th = QM.plot.theme(), pl = probsL(state);
        const hm = { nx: HX, ny: HY, xMin: 0, xMax: 2 * Math.PI, yMin: -1, yMax: 1, piTicks: true, xlabel: 'φ', ylabel: 'cos θ' };
        QM.plot.heatmap(cA, Object.assign({ data: fieldMap(state) }, hm));
        const L = +P.projL;
        pB.head.firstChild.innerHTML = `Projection onto $l=${L}$ (normalised), $P(l{=}${L})=${QM.fmt(pl[L], 4)}$`; QM.ui.tex(pB.head);
        if (pl[L] > 1e-12) QM.plot.heatmap(cB, Object.assign({ data: fieldMap(state, L) }, hm));
        else { const { ctx, w, h } = QM.plot.setup(cB); ctx.fillStyle = th.muted; ctx.textAlign = 'center'; ctx.fillText(`no l = ${L} component`, w / 2, h / 2); }
        QM.plot.bars(cL, { labels: Array.from(pl, (_, l) => 'l=' + l), values: Array.from(pl), yMax: 1, valueLabels: true });
        // table
        let L2 = 0, Lz = 0, Lz2 = 0;
        // show only l up to the highest one carrying weight (keeps the table narrow)
        let LS = 1;
        for (let l = 0; l <= LMAX; l++) for (let m = -l; m <= l; m++) { const j = lmIdx(l, m); if (state.re[j] ** 2 + state.im[j] ** 2 > 1e-10) LS = Math.max(LS, l); }
        for (let l = LS + 1; l <= LMAX; l++) for (let m = -l; m <= l; m++) { const j = lmIdx(l, m), p = state.re[j] ** 2 + state.im[j] ** 2; L2 += l * (l + 1) * p; Lz += m * p; Lz2 += m * m * p; }
        let h = '<table class="matrix"><tr><td></td>' + Array.from({ length: 2 * LS + 1 }, (_, i) => `<td><b>m=${i - LS}</b></td>`).join('') + '</tr>';
        for (let l = 0; l <= LS; l++) {
          h += `<tr><td><b>l=${l}</b></td>`;
          for (let m = -LS; m <= LS; m++) {
            if (Math.abs(m) > l) { h += '<td></td>'; continue; }
            const j = lmIdx(l, m), p = state.re[j] ** 2 + state.im[j] ** 2;
            L2 += l * (l + 1) * p; Lz += m * p; Lz2 += m * m * p;
            h += `<td style="background:${QM.plot.hexA(th.accent, Math.min(0.85, 0.08 + p))}">${p > 1e-10 ? QM.fmt(p, 4) : '·'}</td>`;
          }
          h += '</tr>';
        }
        tbl.innerHTML = h + '</table>';
        amps.textContent = 'ψ_ang = ' + ampText();
        st.set('L2', QM.fmt(L2, 4) + ' ħ²'); st.set('Lz', QM.fmt(Lz, 4) + ' ħ'); st.set('dLz', QM.fmt(Math.sqrt(Math.max(0, Lz2 - Lz * Lz)), 4) + ' ħ');
        st.set('nrm', QM.fmt(raw.norm, 5)); st.set('rest', QM.fmt(rest, 6));
        mOut.innerHTML = collapsedMsg;
        rOut.innerHTML = (P.radial ? `$\\psi(r,\\theta,\\varphi)=R(r)\\,f(\\theta,\\varphi)$ with ${P.radial}.<br>` : 'Only the angular part is modelled; any radial factor $R(r)$ multiplies it unchanged.<br>') +
          'Measuring $L^2$ or $L_z$ acts on the angular part only – the radial factor is unaffected. Current angular part: <span class="mono">' + ampText() + '</span>';
        QM.ui.tex(rOut);
      }

      rebuild();
      const off = QM.ui.onResize(views, redraw);
      return { destroy() { off(); }, redraw };
    },
  });

  // ---------- self-tests ----------
  QM.addTest('angular: (x+y)² → P(l=0) = 5/9, P(l=2) = 4/9; m=0: 1/9, m=±2: 1/6, m=±1: 0', () => {
    const r = decompose(parse('(x+y)^2', {})), p = (l, m) => (r.re[lmIdx(l, m)] ** 2 + r.im[lmIdx(l, m)] ** 2) / r.norm;
    const p0 = p(0, 0); let p2 = 0; for (let m = -2; m <= 2; m++) p2 += p(2, m);
    const ok = Math.abs(p0 - 5 / 9) < 1e-10 && Math.abs(p2 - 4 / 9) < 1e-10 && Math.abs(p(2, 0) - 1 / 9) < 1e-10 && Math.abs(p(2, 2) - 1 / 6) < 1e-10 && Math.abs(p(2, -2) - 1 / 6) < 1e-10 && p(2, 1) < 1e-12 && p(2, -1) < 1e-12;
    return { pass: ok, detail: `P(0) = ${p0.toFixed(8)}, P(2) = ${p2.toFixed(8)}, P(2,0) = ${p(2, 0).toFixed(6)}, P(2,±2) = ${p(2, 2).toFixed(6)}, ${p(2, -2).toFixed(6)}` };
  });
  QM.addTest('angular: Y_lm orthonormal for l ≤ 3 (quadrature)', () => {
    let worst = 0;
    for (let l = 0; l <= 3; l++) for (let m = -l; m <= l; m++) {
      const c = project((u, ph) => Ylm(l, m, u, ph));
      for (let l2 = 0; l2 <= 3; l2++) for (let m2 = -l2; m2 <= l2; m2++) {
        const j = lmIdx(l2, m2), want = l === l2 && m === m2 ? 1 : 0;
        worst = Math.max(worst, Math.hypot(c.re[j] - want, c.im[j]));
      }
    }
    return { pass: worst < 1e-12, detail: `max |⟨Y_l'm'|Y_lm⟩ − δ| = ${worst.toExponential(2)}` };
  });
  QM.addTest('angular: Condon–Shortley phase, Y_11 = −√(3/8π) sinθ e^{iφ}', () => {
    const u = 0.3, ph = 0.7, y = Ylm(1, 1, u, ph), want = -Math.sqrt(3 / (8 * Math.PI)) * Math.sqrt(1 - u * u);
    const ok = Math.abs(y[0] - want * Math.cos(ph)) < 1e-14 && Math.abs(y[1] - want * Math.sin(ph)) < 1e-14;
    return { pass: ok, detail: `Y11 = ${y[0].toFixed(6)} + ${y[1].toFixed(6)}i` };
  });
})();
