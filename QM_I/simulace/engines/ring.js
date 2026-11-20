/* Engine 'ring' – particle on a ring / quantum rotor sandbox (ħ = 1).
   Truncated angular-momentum basis |n⟩, n = -nMax..nMax, vector index k = n + nMax, ψ(φ) = (1/√2π) Σ a_n e^{inφ}.
   Hamiltonians (S.ham):
     'free' : H = L²/2I                                  E_n = n²/2I, analytic phases, revival T = 4πI
     'flux' : H = (L − Φ)²/2I  (Aharonov–Bohm flux Φ in flux quanta)   E_n = (n−Φ)²/2I, still diagonal
     'pot'  : H = L²/2I + V1 cos(φ − φ1) + V2 cos 2φ    ⟨n±1|V|n⟩ = V1/2 e^{∓iφ1}, ⟨n±2|V|n⟩ = V2/2
              → diagonalised with QM.eigh, evolution in the eigenbasis.
   Operators: L = diag(n), P = diag((-1)^n), Â† = d e^{iφ}: Â†|n⟩ = d|n+1⟩, Â = (Â†)†, x̂ = (Â†+Â)/2, ŷ = (Â†−Â)/2i, H = current Ĥ.

   PRESET SCHEMA (all keys optional; old presets without `ham` are the free rotor):
   { nMax: 6            truncation |n| ≤ nMax (1..20)
     I: 0.5             moment of inertia (I = ½ → E_n = n², free revival T = 4πI)
     d: 1               ring radius (x̂ = d cos φ̂, ŷ = d sin φ̂, Â† = d e^{iφ})
     ham: 'free' | 'pot' | 'flux'
       pot : V1: 4, phi1: 0, V2: 0        V(φ) = V1 cos(φ−φ1) + V2 cos 2φ   (V1 < 0: minimum at φ1 – pendulum)
       flux: flux: 0                     Φ in flux quanta (spectrum periodic in Φ with period 1)
     tAuto: true        time window (slider + carpet) = free revival period 4πI;  false → tWin: 2π (user window)
     builder: 'gauss' | 'single' | 'pm' | 'cos' | 'sin' | 'xiy' | 'comb' | 'eigen' | 'super' | 'cat' | 'box' | 'random' | 'custom'
       gauss : phi0 (centre), sigma (angular width), n0 (kick = mean L, may be non-integer)
               a_n ∝ exp(−σ²(n−n0)²/2) e^{−inφ0}  (wrapped Gaussian)
       single: nSingle          → |n⟩
       pm    : nPm, alpha       → |n⟩ + e^{iα}|−n⟩
       cos / sin                → ψ ∝ cos φ / sin φ
       xiy   : cRe, cIm         → (x̂ + c·iŷ)|0⟩ ∝ (1+c)/2 |1⟩ + (1−c)/2 |−1⟩
       comb  : K                → all |n| ≤ K equal (ψ ∝ Dirichlet kernel D_K)
       eigen : kEig             → k-th eigenstate of the current Ĥ (k = 0 ground state, sorted by energy)
       super : eA, eB, eGamma   → (|E_A⟩ + e^{iγ}|E_B⟩)/√2
       cat   : phi0, sigma, n0, catSep, catPhase → Gaussian at φ0 + e^{iχ}·Gaussian at φ0 + sep
       box   : phi0, boxW, n0   → top-hat ψ ∝ e^{in0φ} on the arc |φ−φ0| < w/2 (analytic sinc coefficients)
       random: seed, rK         → seeded random complex a_n for |n| ≤ rK
       custom: amps: [[n, re, im], ...]
     view: 'reim' | 'dens' | 'phase'   ψ(φ,t) plot mode
     tab: 'matrix' | 'measure'         lower panel tab;   matOp: 'x'|'y'|'L'|'Ap'|'Am'|'P'|'H'
     showDirichlet: false   true → Dirichlet-kernel panel on top (else at the bottom);  dirichletN: 10
     t0: 0, play: false, probePhi: π/2 }                                                          */
(function () {
  'use strict';
  const { V, M, C } = QM;
  const TAU = 2 * Math.PI, INV = 1 / Math.sqrt(TAU);
  const fmtN = (n) => (n < 0 ? '−' + -n : String(n));
  const nsOf = (nMax) => { const a = []; for (let n = -nMax; n <= nMax; n++) a.push(n); return a; };

  // ---------- physics (pure helpers, also used by self-tests) ----------
  /** analytic evolution for diagonal H: E_n = (n − Φ)²/2I */
  function evolveC(a, nMax, I, t, Phi = 0) {
    const c = V.zeros(a.n);
    for (let k = 0; k < a.n; k++) {
      const n = k - nMax, th = (-(n - Phi) * (n - Phi) / (2 * I)) * t, co = Math.cos(th), si = Math.sin(th);
      c.re[k] = a.re[k] * co - a.im[k] * si; c.im[k] = a.re[k] * si + a.im[k] * co;
    }
    return c;
  }
  function psiAt(c, nMax, phi) {
    let re = 0, im = 0;
    for (let k = 0; k < c.n; k++) {
      const cr = c.re[k], ci = c.im[k]; if (cr === 0 && ci === 0) continue;
      const th = (k - nMax) * phi, co = Math.cos(th), si = Math.sin(th);
      re += cr * co - ci * si; im += cr * si + ci * co;
    }
    return { re: re * INV, im: im * INV };
  }
  /** moments of the (normalised) coefficient vector c (H = free-rotor value L²/2I) */
  function moments(c, nMax, I, d) {
    let L = 0, L2 = 0, pp = 0, pm = 0, p0 = 0, er = 0, ei = 0;
    for (let k = 0; k < c.n; k++) {
      const n = k - nMax, p = c.re[k] * c.re[k] + c.im[k] * c.im[k];
      L += n * p; L2 += n * n * p; if (n > 0) pp += p; else if (n < 0) pm += p; else p0 += p;
      if (k + 1 < c.n) { // <e^{iφ}> = Σ conj(c_{n+1}) c_n
        er += c.re[k + 1] * c.re[k] + c.im[k + 1] * c.im[k];
        ei += c.re[k + 1] * c.im[k] - c.im[k + 1] * c.re[k];
      }
    }
    return { L, dL: Math.sqrt(Math.max(0, L2 - L * L)), H: L2 / (2 * I), pp, pm, p0, er, ei, x: d * er, y: d * ei, mag: Math.hypot(er, ei) };
  }
  function makeOps(nMax, d, I) {
    const N = 2 * nMax + 1, L = M.zeros(N), H = M.zeros(N), P = M.zeros(N), Ap = M.zeros(N);
    for (let k = 0; k < N; k++) {
      const n = k - nMax;
      L.re[k * N + k] = n; H.re[k * N + k] = n * n / (2 * I); P.re[k * N + k] = n % 2 ? -1 : 1;
      if (k + 1 < N) Ap.re[(k + 1) * N + k] = d;
    }
    const Am = M.dagger(Ap);
    return { L, H, P, Ap, Am, x: M.scale(M.add(Ap, Am), 0.5), y: M.scale(M.sub(Ap, Am), { re: 0, im: -0.5 }) };
  }
  const potV = (S, phi) => (S.ham === 'pot' ? S.V1 * Math.cos(phi - S.phi1) + S.V2 * Math.cos(2 * phi) : 0);
  /** Hamiltonian + eigensystem. numeric = true forces QM.eigh even for diagonal H (tests).
      returns {H, diag, Phi, eig: {values ascending, vectors}} */
  function makeSystem(S, nMax, numeric) {
    const N = 2 * nMax + 1, I = S.I, ham = S.ham || 'free', Phi = ham === 'flux' ? +S.flux || 0 : 0, H = M.zeros(N);
    for (let k = 0; k < N; k++) H.re[k * N + k] = (k - nMax - Phi) ** 2 / (2 * I);
    if (ham === 'pot') {
      const h1 = S.V1 / 2, h2 = S.V2 / 2, c = Math.cos(S.phi1), s = Math.sin(S.phi1);
      for (let k = 0; k + 1 < N; k++) {
        H.re[(k + 1) * N + k] = h1 * c; H.im[(k + 1) * N + k] = -h1 * s; // ⟨n+1|V|n⟩ = V1/2 e^{−iφ1}
        H.re[k * N + k + 1] = h1 * c; H.im[k * N + k + 1] = h1 * s;
        if (k + 2 < N) { H.re[(k + 2) * N + k] = h2; H.re[k * N + k + 2] = h2; }
      }
    }
    const diag = ham !== 'pot' && !numeric;
    let eig;
    if (diag) {
      const idx = [...Array(N).keys()].sort((p, q) => H.re[p * N + p] - H.re[q * N + q] || p - q);
      eig = { values: idx.map((k) => H.re[k * N + k]), vectors: idx.map((k) => V.basis(N, k)), idx };
    } else eig = QM.eigh(H);
    return { H, diag, Phi, eig, N };
  }
  function evolveSys(sys, a, nMax, I, t) { return sys.diag ? evolveC(a, nMax, I, t, sys.Phi) : QM.evolve(sys.eig, a, t); }
  function rng(seed) { // mulberry32
    let s = (Math.round(seed) * 2654435761) >>> 0 || 1;
    return () => { s = (s + 0x6d2b79f5) >>> 0; let z = s; z = Math.imul(z ^ (z >>> 15), z | 1); z ^= z + Math.imul(z ^ (z >>> 7), z | 61); return ((z ^ (z >>> 14)) >>> 0) / 4294967296; };
  }
  const sinc = (x) => (Math.abs(x) < 1e-9 ? 1 : Math.sin(x) / x);
  function buildState(S, nMax, sys) {
    const v = V.zeros(2 * nMax + 1);
    const put = (n, re, im = 0) => { if (Math.abs(n) <= nMax) { v.re[n + nMax] += re; v.im[n + nMax] += im; } };
    const gauss = (n) => Math.exp(-0.5 * S.sigma * S.sigma * (n - S.n0) ** 2);
    const eigVec = (k) => { sys = sys || makeSystem(S, nMax); return sys.eig.vectors[Math.max(0, Math.min(sys.N - 1, Math.round(k)))]; };
    switch (S.builder) {
      case 'gauss':
        for (let n = -nMax; n <= nMax; n++) { const A = gauss(n); put(n, A * Math.cos(n * S.phi0), -A * Math.sin(n * S.phi0)); }
        break;
      case 'single': put(Math.round(S.nSingle), 1); break;
      case 'pm': { const n = Math.abs(Math.round(S.nPm)); if (n === 0) put(0, 1); else { put(n, 1); put(-n, Math.cos(S.alpha), Math.sin(S.alpha)); } break; }
      case 'cos': put(1, 1); put(-1, 1); break;
      case 'sin': put(1, 0, -1); put(-1, 0, 1); break;
      case 'xiy': put(1, (1 + S.cRe) / 2, S.cIm / 2); put(-1, (1 - S.cRe) / 2, -S.cIm / 2); break;
      case 'comb': { const K = Math.round(S.K); for (let n = -K; n <= K; n++) put(n, 1); break; }
      case 'eigen': { const e = eigVec(S.kEig); V.addTo(v, e, { re: 1, im: 0 }); break; }
      case 'super': { V.addTo(v, eigVec(S.eA), { re: 1, im: 0 }); V.addTo(v, eigVec(S.eB), C.expi(S.eGamma)); break; }
      case 'cat':
        for (let n = -nMax; n <= nMax; n++) { // G(n) [e^{−inφ0} + e^{iχ} e^{−in(φ0+sep)}]
          const A = gauss(n), t1 = -n * S.phi0, t2 = S.catPhase - n * (S.phi0 + S.catSep);
          put(n, A * (Math.cos(t1) + Math.cos(t2)), A * (Math.sin(t1) + Math.sin(t2)));
        }
        break;
      case 'box': // a_n ∝ ∫_arc e^{i(n0−n)φ} dφ = w e^{i(n0−n)φ0} sinc((n0−n)w/2)
        for (let n = -nMax; n <= nMax; n++) { const q = S.n0 - n, A = sinc((q * S.boxW) / 2); put(n, A * Math.cos(q * S.phi0), A * Math.sin(q * S.phi0)); }
        break;
      case 'random': {
        const r = rng(S.seed), K = Math.min(nMax, Math.round(S.rK));
        for (let n = -K; n <= K; n++) { const rr = Math.sqrt(-2 * Math.log(1 - r())), th = TAU * r(); put(n, rr * Math.cos(th), rr * Math.sin(th)); }
        break;
      }
      default: (S.amps || []).forEach(([n, re, im]) => put(n, re, im || 0));
    }
    if (V.norm2(v) < 1e-14) put(0, 1);
    return V.normalize(v);
  }
  function dirichlet(x, N) {
    const a = x - TAU * Math.round(x / TAU);
    if (Math.abs(a) < 1e-9) return (2 * N + 1) / TAU;
    return Math.sin((N + 0.5) * a) / (TAU * Math.sin(a / 2));
  }
  QM.ring = { evolveC, psiAt, moments, makeOps, buildState, dirichlet, makeSystem, evolveSys, potV };

  const DEFAULTS = { nMax: 6, I: 0.5, d: 1, ham: 'free', V1: 4, phi1: 0, V2: 0, flux: 0, tAuto: true, tWin: TAU,
    builder: 'gauss', phi0: Math.PI / 2, sigma: 0.6, n0: 2, nSingle: 1, nPm: 2, alpha: 0, cRe: -3, cIm: 0, K: 3,
    kEig: 0, eA: 0, eB: 1, eGamma: 0, catSep: Math.PI, catPhase: 0, boxW: Math.PI / 2, seed: 1, rK: 4,
    view: 'reim', tab: 'matrix', matOp: 'x', showDirichlet: false, dirichletN: 10, t0: 0, play: false, probePhi: Math.PI / 2 };

  QM.registerEngine('ring', {
    title: 'Particle-on-a-ring lab',
    mount(host, preset) {
      const { el, panel, params, button, buttons, canvas, tex } = QM.ui;
      const S = Object.assign({}, DEFAULTS, preset || {});
      let nMax = Math.max(1, Math.min(20, Math.round(S.nMax))); S.nMax = nMax;
      let Icur = S.I, sys = makeSystem(S, nMax), a = buildState(S, nMax, sys), cur = a, ops = null;
      let carpet = null, probePhi = S.probePhi, mode = S.view, lowTab = S.tab === 'measure' ? 1 : 0;
      let psiGeom = null, carGeom = null, ringGeom = null, animTimer = null;
      const basisLabels = nsOf(nMax).map(fmtN);
      const tRev = () => 4 * Math.PI * S.I;
      const win = () => (S.tAuto ? tRev() : Math.max(0.05, S.tWin));
      const rebuildOps = () => { ops = makeOps(nMax, S.d, S.I); ops.H = sys.H; };
      const back = (w) => evolveSys(sys, w, nMax, S.I, -tr.t); // state at t → coefficients at t = 0
      rebuildOps();
      const { controls, views } = QM.ui.layout(host);

      // ---------- controls: Hamiltonian ----------
      const pH = panel('Hamiltonian');
      const hSel = params([{ key: 'ham', label: 'Hamiltonian $\\hat H$', type: 'select', options: [
        { value: 'free', label: 'Free rotor  L²/2I' }, { value: 'pot', label: 'Rotor in a potential V(φ) (pendulum / hindered)' },
        { value: 'flux', label: 'Flux Φ through the ring (Aharonov–Bohm)' }] }], S, () => { showHam(); onHam(); });
      const hp = {
        pot: params([
          { key: 'V1', label: '$V_1$ in $V_1\\cos(\\varphi-\\varphi_1)$', type: 'range', min: -30, max: 30, step: 0.1, hint: '$V_1<0$: minimum at $\\varphi_1$ (pendulum, $\\varphi_1 = 3\\pi/2$ hangs down)' },
          { key: 'phi1', label: '$\\varphi_1$', type: 'range', min: 0, max: TAU, step: 0.01 },
          { key: 'V2', label: '$V_2$ in $V_2\\cos 2\\varphi$', type: 'range', min: -30, max: 30, step: 0.1, hint: 'double well (hindered rotor); couples $n\\leftrightarrow n\\pm2$' }], S, onHam),
        flux: params([{ key: 'flux', label: 'Flux $\\Phi$ (flux quanta)', type: 'range', min: -1, max: 1, step: 0.01, hint: '$\\hat H = (\\hat L-\\Phi)^2/2I$ – minimal coupling $\\hat p\\to\\hat p - qA$ with $\\oint A = \\Phi\\,\\Phi_0$' }], S, onHam),
      };
      const hNote = el('div', { class: 'note' });
      pH.body.append(hSel.el, hp.pot.el, hp.flux.el, hNote);
      function showHam() {
        hp.pot.el.style.display = S.ham === 'pot' ? '' : 'none'; hp.flux.el.style.display = S.ham === 'flux' ? '' : 'none';
        hNote.innerHTML = { free: '$E_n = n^2/2I$, analytic phases $e^{-iE_nt}$.', flux: '$E_n = (n-\\Phi)^2/2I$ – still diagonal in $|n\\rangle$; the parabolas shift with $\\Phi$.',
          pot: 'Matrix in the truncated $|n\\rangle$ basis ($\\langle n\\pm1|V|n\\rangle = \\tfrac{V_1}{2}e^{\\mp i\\varphi_1}$, $\\langle n\\pm2|V|n\\rangle = \\tfrac{V_2}{2}$) diagonalised numerically; $\\psi(t) = \\sum_k e^{-iE_kt}|E_k\\rangle\\langle E_k|\\psi(0)\\rangle$.' }[S.ham] || '';
        tex(hNote);
        revBtns.style.display = S.ham === 'pot' ? 'none' : '';
      }
      function onHam() {
        sys = makeSystem(S, nMax); rebuildOps();
        if (S.builder === 'eigen' || S.builder === 'super') { a = buildState(S, nMax, sys); editor.set(a); }
        renderMatrix(); meas.refresh(); carNote(); onState();
      }
      controls.append(pH.el);

      // ---------- controls: rotor parameters + time ----------
      const pRot = panel('Rotor &amp; time');
      const par = params([
        { key: 'nMax', label: 'Truncation $n_{\\max}$ ($|n|\\le n_{\\max}$)', type: 'range', min: 1, max: 20, step: 1 },
        { key: 'I', label: 'Moment of inertia $I$', type: 'range', min: 0.1, max: 3, step: 0.05, hint: '$E_n = n^2/2I$; free revival at $T = 4\\pi I$' },
        { key: 'd', label: 'Radius $d$', type: 'range', min: 0.2, max: 3, step: 0.05 },
      ], S, onParam);
      const wPar = params([
        { key: 'tAuto', label: 'Time window = free revival period $4\\pi I$', type: 'check' },
        { key: 'tWin', label: 'Time window $T_w$ (slider &amp; carpet)', type: 'range', min: 0.2, max: 100, step: 0.1 }], S, (k) => {
        if (k === 'tWin') { S.tAuto = false; wPar.set('tAuto', false, true); } else if (S.tAuto) wPar.set('tWin', +tRev().toFixed(3), true);
        tr.setMax(win()); computeCarpet(); redraw();
      });
      if (S.tAuto) wPar.set('tWin', +tRev().toFixed(3), true);
      const tr = QM.ui.transport({
        tMax: win(), t0: S.t0 || 0, label: 'Time $t$', onTime: redraw,
        fmt: (t) => QM.fmt(t, 3) + '  (τ = ' + QM.fmt(t / (2 * S.I) / Math.PI, 3) + 'π)',
      });
      const revBtns = buttons(button('t = T/4', () => tr.set(tRev() / 4)), button('T/2 (shift π)', () => tr.set(tRev() / 2)), button('T (revival)', () => tr.set(tRev())));
      pRot.body.append(par.el, wPar.el, el('div', { style: { height: '8px' } }), tr.el, revBtns);
      controls.append(pRot.el);

      // ---------- controls: state builder ----------
      const pB = panel('Initial state');
      const bSel = params([{ key: 'builder', label: 'Builder', type: 'select', options: [
        { value: 'gauss', label: 'Gaussian packet (φ₀, σ, kick n₀)' }, { value: 'single', label: 'Single |n⟩' },
        { value: 'pm', label: '|n⟩ + e^(iα)|−n⟩' }, { value: 'cos', label: 'ψ ∝ cos φ' }, { value: 'sin', label: 'ψ ∝ sin φ' },
        { value: 'xiy', label: '(x̂ + c·iŷ)|0⟩' }, { value: 'comb', label: 'Comb |n| ≤ K (ψ ∝ D_K)' },
        { value: 'eigen', label: 'Eigenstate |E_k⟩ of current Ĥ' }, { value: 'super', label: '(|E_A⟩ + e^(iγ)|E_B⟩)/√2' },
        { value: 'cat', label: 'Cat: two packets (separation, phase)' }, { value: 'box', label: 'Top-hat on an arc' },
        { value: 'random', label: 'Random state (seeded)' }, { value: 'custom', label: 'Custom (edit coefficients)' }] }],
      S, () => { showGroup(); applyBuilder(); });
      const R = (key, label, min, max, step, hint) => ({ key, label, type: 'range', min, max, step, hint });
      const gp = {
        gauss: params([R('phi0', 'Centre $\\varphi_0$', 0, TAU, 0.01), R('sigma', 'Width $\\sigma_\\varphi$', 0.15, 3, 0.01), R('n0', 'Kick $n_0$ (≈ $\\langle L\\rangle$)', -10, 10, 0.1)], S, applyBuilder),
        single: params([R('nSingle', '$n$', -20, 20, 1)], S, applyBuilder),
        pm: params([R('nPm', '$n$', 0, 20, 1), R('alpha', 'Relative phase $\\alpha$', 0, TAU, 0.01)], S, applyBuilder),
        xiy: params([R('cRe', '$\\mathrm{Re}\\,c$', -5, 5, 0.1),
          R('cIm', '$\\mathrm{Im}\\,c$', -5, 5, 0.1, '$\\hat x + c\\,i\\hat y = d\\big[\\tfrac{1+c}{2}e^{i\\varphi} + \\tfrac{1-c}{2}e^{-i\\varphi}\\big]$, so $(\\hat x + c\\,i\\hat y)|0\\rangle \\propto \\tfrac{1+c}{2}|1\\rangle + \\tfrac{1-c}{2}|{-1}\\rangle$')], S, applyBuilder),
        comb: params([R('K', '$K$', 0, 20, 1)], S, applyBuilder),
        eigen: params([R('kEig', 'Level $k$ ($k=0$: ground state)', 0, 40, 1, 'Levels sorted by energy; degenerate pairs in arbitrary order (free: $-n$ before $+n$).')], S, applyBuilder),
        super: params([R('eA', 'Level $A$', 0, 40, 1), R('eB', 'Level $B$', 0, 40, 1), R('eGamma', 'Relative phase $\\gamma$', 0, TAU, 0.01, 'Beats with frequency $E_B - E_A$.')], S, applyBuilder),
        cat: params([R('phi0', 'Centre of first packet $\\varphi_0$', 0, TAU, 0.01), R('catSep', 'Separation $\\Delta\\varphi$', 0, TAU, 0.01), R('sigma', 'Width $\\sigma_\\varphi$', 0.15, 3, 0.01),
          R('n0', 'Kick $n_0$', -10, 10, 0.1), R('catPhase', 'Relative phase $\\chi$', 0, TAU, 0.01, '$\\psi \\propto g(\\varphi-\\varphi_0) + e^{i\\chi}g(\\varphi-\\varphi_0-\\Delta\\varphi)$; $\\Delta\\varphi=\\pi$, $\\chi = 0/\\pi$ → only even / odd $n$.')], S, applyBuilder),
        box: params([R('phi0', 'Arc centre $\\varphi_0$', 0, TAU, 0.01), R('boxW', 'Arc width $w$', 0.05, TAU, 0.01), R('n0', 'Phase gradient $n_0$ ($\\psi\\propto e^{in_0\\varphi}$)', -10, 10, 0.1, '$a_n \\propto e^{i(n_0-n)\\varphi_0}\\,\\mathrm{sinc}\\big((n-n_0)w/2\\big)$ – Gibbs ringing from the truncation.')], S, applyBuilder),
        random: params([R('seed', 'Seed', 1, 999, 1), R('rK', 'Max $|n|$', 0, 20, 1)], S, applyBuilder),
      };
      const gNotes = {
        cos: '$\\psi \\propto \\cos\\varphi$: $a_{\\pm1} = 1/\\sqrt2$ (this is $\\hat x|0\\rangle$).',
        sin: '$\\psi \\propto \\sin\\varphi$: $a_{\\pm1} = \\mp i/\\sqrt2$ (this is $\\hat y|0\\rangle$).',
        custom: 'Edit the coefficients below (they are normalised automatically).',
      };
      const groupEls = {};
      Object.keys(gp).forEach((k) => { groupEls[k] = gp[k].el; });
      Object.keys(gNotes).forEach((k) => { groupEls[k] = el('div', { class: 'note', html: gNotes[k] }); });
      pB.body.append(bSel.el, ...Object.values(groupEls));
      function showGroup() { Object.keys(groupEls).forEach((k) => { groupEls[k].style.display = k === S.builder ? '' : 'none'; }); }
      function applyBuilder(key) {
        if (S.builder === 'custom') return;
        if (['phi0', 'sigma', 'n0'].includes(key)) Object.values(gp).forEach((g) => g.set(key, S[key], true)); // shared by gauss/cat/box
        a = buildState(S, nMax, sys); editor.set(a); onState();
      }
      function markCustom() { S.builder = 'custom'; bSel.set('builder', 'custom', true); showGroup(); }
      controls.append(pB.el);

      // ---------- controls: coefficient editor ----------
      const pE = panel('Coefficients $a_n$ at $t = 0$');
      let editor = null;
      function rebuildEditor() {
        pE.body.innerHTML = '';
        editor = QM.ui.amplitudes({ labels: nsOf(nMax).map((n) => 'n = ' + fmtN(n)), vector: a,
          onChange(v) { a = v; markCustom(); onState(); } });
        pE.body.append(editor.el, el('div', { class: 'note', html: 'Free / flux rotor: each $a_n$ picks up $e^{-iE_nt}$. With a potential the state is evolved in the eigenbasis of $\\hat H$. Ladder operators and measurements act on the state at the current $t$; the result is propagated back to $t=0$ here.' }));
        tex(pE.body);
      }
      rebuildEditor();
      controls.append(pE.el);

      // ---------- controls: ladder + parity ----------
      const pL = panel('Ladder operators &amp; parity');
      const lOut = el('div', { class: 'meas-result' });
      pL.body.append(buttons(button('$\\hat A^\\dagger$ (raise)', () => applyOp('up'), 'primary'), button('$\\hat A$ (lower)', () => applyOp('down')), button('$\\hat P$ (parity)', () => applyOp('P'))),
        lOut, el('div', { class: 'note', html: '$\\hat A^\\dagger = \\hat x + i\\hat y = d\\,e^{i\\varphi}$: $a_n \\to d\\,a_{n-1}$ (whole spectrum shifts up by one). $\\hat A = d\\,e^{-i\\varphi}$ shifts down. Because $\\hat A^\\dagger\\hat A = \\hat A\\hat A^\\dagger = d^2$, the norm is always multiplied by exactly $d$ – unlike the harmonic oscillator there is <b>no lowest state</b>; the ladder is unbounded in both directions (here only cut by the truncation). $\\hat P$: $\\psi(\\varphi) \\to \\psi(\\varphi+\\pi)$, i.e. $a_n \\to (-1)^n a_n$. With a potential or flux these are no longer ladders of $\\hat H$ – they still shift $\\langle L\\rangle$ by $\\pm1$.' }));
      controls.append(pL.el);
      function applyOp(kind) {
        const c = cur, before = moments(c, nMax, S.I, S.d);
        const w = M.apply(kind === 'P' ? ops.P : kind === 'up' ? ops.Ap : ops.Am, c);
        const norm = Math.sqrt(V.norm2(w));
        if (norm < 1e-12) { lOut.textContent = 'Result is the zero vector (state sits entirely at the truncation edge) – unchanged.'; return; }
        const lost = kind === 'P' ? 0 : C.abs2(V.get(c, kind === 'up' ? c.n - 1 : 0));
        V.normalize(w);
        const after = moments(w, nMax, S.I, S.d), name = kind === 'P' ? 'P̂' : kind === 'up' ? 'Â†' : 'Â';
        a = V.normalize(back(w)); markCustom(); editor.set(a);
        lOut.innerHTML = `${name}: ⟨L⟩ ${QM.fmt(before.L, 3)} → <b>${QM.fmt(after.L, 3)}</b>` +
          (kind === 'P' ? '' : `; ‖${name}ψ‖ = ${QM.fmt(norm, 4)} = d · ${QM.fmt(norm / S.d, 4)}`) +
          (lost > 1e-9 ? `<br><span style="color:var(--s5)">Truncation: weight ${QM.fmt(lost, 4)} at n = ${kind === 'up' ? '' : '−'}${nMax} fell off the basis (increase n<sub>max</sub>). The state was renormalised.</span>` : '');
        onState();
      }

      // ---------- views ----------
      const dirP = panel('Dirichlet kernel $D_N(\\varphi) \\to$ periodic $\\delta$');
      if (S.showDirichlet) views.append(dirP.el);

      const modeKeys = ['reim', 'dens', 'phase'];
      const modeTabs = QM.ui.tabs([{ label: 'Re / Im' }, { label: '|ψ|²' }, { label: '|ψ|² + phase' }], (i) => { mode = modeKeys[i]; drawPsi(); drawRing(); });
      const pPsi = panel('Wave function $\\psi(\\varphi, t)$', { right: modeTabs.head });
      const cvPsi = canvas(250), probeNote = el('div', { class: 'note' });
      pPsi.body.append(cvPsi, probeNote);
      views.append(pPsi.el);

      const pRing = panel('Ring: $|\\psi|^2$ band, arrow $\\langle e^{i\\varphi}\\rangle$'), cvRing = canvas(260);
      pRing.body.append(cvRing, el('div', { class: 'note' }, 'Click the ring to move the probe angle.'));
      const pBars = panel('Spectrum $|a_n|^2$ (colour = phase of $a_n(t)$)'), cvBars = canvas(260);
      pBars.body.append(cvBars);
      views.append(el('div', { class: 'row2' }, pRing.el, pBars.el));

      const pSpec = panel('Energy levels &amp; weights $|\\langle E|\\psi\\rangle|^2$'), cvLev = canvas(280), cvW = canvas(280), specGrid = el('div', null, cvLev, cvW), specNote = el('div', { class: 'note' });
      pSpec.body.append(specGrid, specNote);
      views.append(pSpec.el);

      const pCar = panel('Quantum carpet $|\\psi(\\varphi,t)|^2$ – horizontal $\\varphi$, vertical $t \\in [0, T_w]$'), cvCar = canvas(280), carNoteEl = el('div', { class: 'note' });
      pCar.body.append(cvCar, carNoteEl);
      views.append(pCar.el);
      function carNote() {
        carNoteEl.innerHTML = (S.ham === 'free' ? 'Exact revival at $T = 4\\pi I$ ($\\tau = 2\\pi$, all phases $e^{-2\\pi i n^2} = 1$). At $T/2$ ($\\tau=\\pi$): $e^{-i\\pi n^2} = (-1)^n$, so $\\psi(\\varphi, T/2) = \\psi(\\varphi+\\pi, 0)$ – the packet reappears on the opposite side (dashed line).'
          : S.ham === 'flux' ? 'With flux, at $T = 4\\pi I$ the phases are $e^{-2\\pi i(n-\\Phi)^2} = e^{-2\\pi i\\Phi^2}e^{4\\pi i n\\Phi}$: $\\psi(\\varphi,T) \\propto \\psi(\\varphi + 4\\pi\\Phi, 0)$ – the revival is <b>rotated</b> by $-4\\pi\\Phi$ (Aharonov–Bohm phase); at $T/2$ the packet sits at $\\varphi_0 + \\pi - 2\\pi\\Phi$.'
            : 'With a potential the levels are incommensurate – no exact revivals in general. Set the window $T_w$ yourself (e.g. several classical periods $2\\pi/\\omega$).') + ' Click to set $\\varphi$ and $t$.';
        tex(carNoteEl);
      }

      const pSt = panel('Expectation values at time $t$');
      const st = QM.ui.stats([
        { key: 'L', label: '$\\langle L\\rangle$' }, { key: 'dL', label: '$\\Delta L$' }, { key: 'H', label: '$\\langle H\\rangle$' }, { key: 'dH', label: '$\\Delta H$' },
        { key: 'x', label: '$\\langle x\\rangle = d\\,\\mathrm{Re}\\langle e^{i\\varphi}\\rangle$' }, { key: 'y', label: '$\\langle y\\rangle = d\\,\\mathrm{Im}\\langle e^{i\\varphi}\\rangle$' },
        { key: 'e', label: '$|\\langle e^{i\\varphi}\\rangle|$' }, { key: 'pp', label: '$P(\\ell>0)$' }, { key: 'pm', label: '$P(\\ell<0)$' },
        { key: 'p0', label: '$P(\\ell=0)$' }, { key: 'tau', label: '$\\tau = t/2I$' }]);
      pSt.body.append(st.el);
      views.append(pSt.el);

      // lower: matrix elements | measurement
      const matEl = el('div'), measEl = el('div', { style: { display: 'none' } });
      const low = QM.ui.tabs([{ label: 'Matrix elements', el: matEl }, { label: 'Measurement', el: measEl }], (i) => { lowTab = i; if (i === 1) meas.update(); });
      const pLow = panel('Operators in the $|n\\rangle$ basis', { right: low.head });
      pLow.body.append(matEl, measEl);
      views.append(pLow.el);
      if (!S.showDirichlet) views.append(dirP.el);

      const matSel = params([{ key: 'matOp', label: 'Operator', type: 'select', options: [
        { value: 'x', label: 'x̂ = d cos φ̂' }, { value: 'y', label: 'ŷ = d sin φ̂' }, { value: 'L', label: 'L̂ = −i d/dφ' },
        { value: 'Ap', label: 'Â† = x̂ + iŷ = d e^(iφ)' }, { value: 'Am', label: 'Â = x̂ − iŷ = d e^(−iφ)' }, { value: 'P', label: 'P̂ (parity)' }, { value: 'H', label: 'Ĥ (current Hamiltonian)' }] }],
      S, renderMatrix);
      const matBox = el('div', { class: 'matrix-wrap' }), matInfo = el('div', { class: 'note' }), matChecks = el('div', { class: 'note' });
      matEl.append(matSel.el, matBox, matInfo, matChecks);
      const MAT_INFO = {
        x: '$\\langle m\'|\\hat x|m\\rangle = \\tfrac d2(\\delta_{m\',m+1}+\\delta_{m\',m-1})$. The true spectrum of $\\hat x$ is continuous, $[-d,d]$, each inner value doubly degenerate ($\\varphi$ and $2\\pi-\\varphi$); the truncated matrix only has the $2n_{\\max}+1$ eigenvalues $d\\cos\\frac{k\\pi}{2n_{\\max}+2}$.',
        y: '$\\langle m\'|\\hat y|m\\rangle = \\tfrac{d}{2i}(\\delta_{m\',m+1}-\\delta_{m\',m-1})$ – Hermitian, purely imaginary. Truncation caveat as for $\\hat x$.',
        L: '$\\langle m\'|\\hat L|m\\rangle = m\\,\\delta_{m\'m}$ – the basis is the eigenbasis of $\\hat L$.',
        Ap: '$\\langle m\'|\\hat A^\\dagger|m\\rangle = d\\,\\delta_{m\',m+1}$, i.e. $\\hat A^\\dagger|m\\rangle = d|m+1\\rangle$. $[\\hat L,\\hat A^\\dagger] = \\hat A^\\dagger$. No state is annihilated – no lowest rung.',
        Am: '$\\langle m\'|\\hat A|m\\rangle = d\\,\\delta_{m\',m-1}$, $[\\hat L,\\hat A] = -\\hat A$.',
        P: '$\\hat P = e^{i\\pi\\hat L}$: $\\langle m\'|\\hat P|m\\rangle = (-1)^m\\delta_{m\'m}$. Commutes with the free $\\hat L,\\hat H$; anticommutes with $\\hat x,\\hat y$.',
        H: () => (S.ham === 'pot' ? '$\\langle m\'|\\hat H|m\\rangle = \\frac{m^2}{2I}\\delta_{m\'m} + \\frac{V_1}{2}\\big(e^{-i\\varphi_1}\\delta_{m\',m+1} + e^{i\\varphi_1}\\delta_{m\',m-1}\\big) + \\frac{V_2}{2}(\\delta_{m\',m+2}+\\delta_{m\',m-2})$ – a band matrix; $\\hat L$ is no longer conserved.'
          : S.ham === 'flux' ? '$\\langle m\'|\\hat H|m\\rangle = \\frac{(m-\\Phi)^2}{2I}\\delta_{m\'m}$ – the $\\pm n$ degeneracy is lifted unless $\\Phi$ is an integer or half-integer (then $n \\leftrightarrow 2\\Phi - n$).'
            : '$\\langle m\'|\\hat H|m\\rangle = \\frac{m^2}{2I}\\delta_{m\'m}$ – levels $n\\ne0$ doubly degenerate ($\\pm n$).'),
      };
      function renderMatrix() {
        const K = Math.min(nMax, 4), op = ops[S.matOp], th = { style: { color: 'var(--muted)' } };
        const tbl = el('table', { class: 'matrix' });
        const hr = el('tr', null, el('td', th, ''));
        for (let m = -K; m <= K; m++) hr.append(el('td', th, '|' + fmtN(m) + '⟩'));
        tbl.append(hr);
        for (let mp = -K; mp <= K; mp++) {
          const row = el('tr', null, el('td', th, '⟨' + fmtN(mp) + '|'));
          for (let m = -K; m <= K; m++) { const z = M.get(op, mp + nMax, m + nMax); row.append(el('td', null, C.abs(z) < 1e-12 ? '·' : C.fmt(z, 3))); }
          tbl.append(row);
        }
        matBox.innerHTML = ''; matBox.append(tbl);
        const info = MAT_INFO[S.matOp];
        matInfo.innerHTML = (K < nMax ? `Showing $m, m' \\in [-${K}, ${K}]$ of the ${2 * nMax + 1}-dimensional truncation. ` : '') + (typeof info === 'function' ? info() : info);
        const e = (A, B) => { const v = M.maxAbs(M.sub(A, B)); return v < 1e-12 ? '0' : v.toExponential(1); };
        const AdA = M.mul(ops.Ap, ops.Am);
        matChecks.innerHTML = 'Numerical checks on the truncated matrices: ' +
          `$\\|[\\hat L,\\hat A^\\dagger]-\\hat A^\\dagger\\| = ${e(M.comm(ops.L, ops.Ap), ops.Ap)}$, ` +
          `$\\|[\\hat L,\\hat x]-i\\hat y\\| = ${e(M.comm(ops.L, ops.x), M.scale(ops.y, { re: 0, im: 1 }))}$, ` +
          `$\\|[\\hat L,\\hat y]+i\\hat x\\| = ${e(M.comm(ops.L, ops.y), M.scale(ops.x, { re: 0, im: -1 }))}$, ` +
          `$\\|\\hat P\\hat x+\\hat x\\hat P\\| = ${e(M.mul(ops.P, ops.x), M.scale(M.mul(ops.x, ops.P), -1))}$, ` +
          `$\\|[\\hat H,\\hat L]\\| = ${e(M.comm(ops.H, ops.L), M.zeros(ops.L.n))}$. ` +
          `$\\hat A^\\dagger\\hat A$ has diagonal $d^2 = ${QM.fmt(S.d * S.d, 3)}$ everywhere except the corner $m=-n_{\\max}$ (value ${QM.fmt(AdA.re[0], 3)}) – a truncation artefact.`;
        tex(matEl);
      }

      const meas = QM.measurePanel({
        observables: () => [
          { key: 'L', label: 'L̂ (angular momentum)', op: ops.L }, { key: 'H', label: 'Ĥ (energy, current Hamiltonian)', op: ops.H },
          { key: 'P', label: 'P̂ (parity)', op: ops.P }, { key: 'x', label: 'x̂ = d cos φ̂ (truncated!)', op: ops.x }, { key: 'y', label: 'ŷ = d sin φ̂ (truncated!)', op: ops.y }],
        getState: () => cur,
        onCollapse: (psi) => { a = V.normalize(back(psi)); markCustom(); editor.set(a); onState(); },
        basisLabels,
      });
      measEl.append(meas.el, el('div', { class: 'note', html: 'Measurement acts on $\\psi(t)$ in the truncated basis. $\\hat L,\\hat P$ are diagonal; $\\hat H$ is the current Hamiltonian (numerically diagonalised with a potential). $\\hat x,\\hat y$ have a continuous spectrum on the full space – their truncated matrices give only a discrete approximation.' }));

      // ---------- Dirichlet kernel ----------
      const dirV = { N: Math.round(S.dirichletN) };
      const dirPar = params([{ key: 'N', label: 'Number of frequencies $N$', type: 'range', min: 0, max: 200, step: 1 }], dirV, () => { dirV.N = Math.round(dirV.N); drawDir(); });
      const animBtn = button('▶ Animate N', toggleAnim);
      const cvDir = canvas(300);
      const dirSt = QM.ui.stats([{ key: 'terms', label: 'terms $2N+1$' }, { key: 'peak', label: 'peak $D_N(0) = \\frac{2N+1}{2\\pi}$' }, { key: 'int', label: '$\\int_{-\\pi}^{\\pi} D_N\\,d\\varphi$' }, { key: 'w', label: 'first zero $\\frac{2\\pi}{2N+1}$' }]);
      dirP.body.append(
        el('div', { class: 'statement', html: '$$D_N(\\varphi) = \\frac{1}{2\\pi}\\sum_{n=-N}^{N} e^{in\\varphi} = \\frac{\\sin\\big((N+\\tfrac12)\\varphi\\big)}{2\\pi\\sin(\\varphi/2)} \\;\\xrightarrow{N\\to\\infty}\\; \\delta_{2\\pi}(\\varphi) = \\sum_{k}\\delta(\\varphi-2\\pi k) = \\frac{1}{2\\pi}\\sum_{n\\in\\mathbb Z} e^{in\\varphi}$$' }),
        dirPar.el, buttons(animBtn, button('−1', () => setDirN(dirV.N - 1)), button('+1', () => setDirN(dirV.N + 1)), button('Load ψ ∝ D_N into the ring', loadComb)),
        dirSt.el, cvDir,
        el('div', { class: 'note', html: 'Dashed: $\\varphi = 2\\pi k$, positions of the periodic deltas. $D_N$ does not converge pointwise; it converges <i>as a distribution</i>: $\\int_{-\\pi}^{\\pi} D_N(\\varphi) f(\\varphi)\\,d\\varphi = \\sum_{|n|\\le N}\\hat f_n \\to f(0)$ for smooth periodic $f$. On the ring this is the completeness relation $\\sum_n \\psi_n(\\varphi)\\psi_n^*(\\varphi\') = \\delta_{2\\pi}(\\varphi-\\varphi\')$; the state with equal $a_n$ for $|n|\\le N$ has $\\psi = \\sqrt{2\\pi/(2N+1)}\\,D_N$.' }));
      function setDirN(n) { dirV.N = Math.max(0, Math.min(200, Math.round(n))); dirPar.set('N', dirV.N, true); drawDir(); }
      function stopAnim() { if (animTimer !== null) { clearInterval(animTimer); animTimer = null; } animBtn.textContent = '▶ Animate N'; }
      function toggleAnim() {
        if (animTimer !== null) { stopAnim(); return; }
        if (dirV.N >= 200) setDirN(0);
        animBtn.textContent = '❚❚ Pause';
        animTimer = setInterval(() => { if (dirV.N >= 200) { stopAnim(); return; } setDirN(dirV.N + 1); }, 140);
      }
      function loadComb() {
        S.K = Math.min(dirV.N, nMax); S.builder = 'comb';
        bSel.set('builder', 'comb', true); gp.comb.set('K', S.K, true); showGroup(); applyBuilder();
      }
      function drawDir() {
        const N = dirV.N, X = 6 * Math.PI, th = QM.plot.theme(), top = (2 * N + 1) / TAU;
        const ns = Math.max(Math.ceil((cvDir.clientWidth || 600) * 2.2), Math.min(10000, 50 * (2 * N + 1)));
        const xs = new Float64Array(ns + 1), ys = new Float64Array(ns + 1);
        for (let j = 0; j <= ns; j++) { xs[j] = -X + (2 * X * j) / ns; ys[j] = dirichlet(xs[j], N); }
        const vl = []; for (let k = -3; k <= 3; k++) vl.push({ x: TAU * k, color: th.series[1], dash: [4, 5] });
        QM.plot.line(cvDir, { xMin: -X, xMax: X, yMin: -Math.max(0.3, top * 0.29), yMax: Math.max(0.4, top * 1.12), piTicks: true, xlabel: 'φ',
          series: [{ x: xs, y: ys, color: th.accent, width: 1.6, label: 'D_N(φ)' }], vlines: vl });
        // Simpson integral over [−π, π]
        const m = 2 * Math.max(1000, 20 * (2 * N + 1)), hh = TAU / m; let s = dirichlet(-Math.PI, N) + dirichlet(Math.PI, N);
        for (let j = 1; j < m; j++) s += (j % 2 ? 4 : 2) * dirichlet(-Math.PI + j * hh, N);
        dirSt.set('terms', String(2 * N + 1)); dirSt.set('peak', QM.fmt(top, 4)); dirSt.set('int', QM.fmt((s * hh) / 3, 6)); dirSt.set('w', QM.fmt(TAU / (2 * N + 1), 4));
      }

      // ---------- computations ----------
      function computeCarpet() {
        const nx = 200, ny = 160, N = 2 * nMax + 1, T = win(), data = new Float32Array(nx * ny);
        const cs = new Float64Array(N * nx), sn = new Float64Array(N * nx);
        for (let k = 0; k < N; k++) for (let i = 0; i < nx; i++) { const th = (k - nMax) * TAU * i / (nx - 1); cs[k * nx + i] = Math.cos(th); sn[k * nx + i] = Math.sin(th); }
        const act = []; for (let k = 0; k < N; k++) if (!sys.diag || a.re[k] * a.re[k] + a.im[k] * a.im[k] > 1e-16) act.push(k);
        const ev = sys.diag ? (t) => evolveC(a, nMax, S.I, t, sys.Phi) : QM.evolver(sys.eig, a);
        let vMax = 0;
        for (let j = 0; j < ny; j++) {
          const c = ev((T * j) / (ny - 1));
          for (let i = 0; i < nx; i++) {
            let re = 0, im = 0;
            for (const k of act) { const q = k * nx + i; re += c.re[k] * cs[q] - c.im[k] * sn[q]; im += c.re[k] * sn[q] + c.im[k] * cs[q]; }
            const p = (re * re + im * im) / TAU; data[j * nx + i] = p; if (p > vMax) vMax = p;
          }
        }
        carpet = { data, nx, ny, vMax: vMax || 1 / TAU, T };
      }
      function onParam(key) {
        if (key === 'nMax') {
          const nn = Math.max(1, Math.min(20, Math.round(S.nMax))); S.nMax = nn;
          if (nn === nMax) return;
          const old = nMax; nMax = nn; sys = makeSystem(S, nMax);
          if (S.builder === 'custom') {
            const b = V.zeros(2 * nn + 1);
            for (let n = -Math.min(nn, old); n <= Math.min(nn, old); n++) { b.re[n + nn] = a.re[n + old]; b.im[n + nn] = a.im[n + old]; }
            if (V.norm2(b) < 1e-14) b.re[nn] = 1;
            a = V.normalize(b);
          } else a = buildState(S, nMax, sys);
          basisLabels.length = 0; nsOf(nMax).forEach((n) => basisLabels.push(fmtN(n)));
          rebuildEditor();
        } else if (key === 'I') {
          const tau = tr.t / (2 * Icur); Icur = S.I;
          sys = makeSystem(S, nMax);
          if (S.tAuto) wPar.set('tWin', +tRev().toFixed(3), true);
          tr.setMax(win()); tr.set(Math.min(win(), tau * 2 * S.I), true);
          if (S.builder === 'eigen' || S.builder === 'super') { a = buildState(S, nMax, sys); editor.set(a); }
        }
        rebuildOps();
        renderMatrix(); meas.refresh(); onState();
      }
      function onState() { computeCarpet(); redraw(); }

      // ---------- drawing ----------
      function redraw() {
        cur = evolveSys(sys, a, nMax, S.I, tr.t);
        drawPsi(); drawRing(); drawBars(); drawSpec(); drawCarpet(); updateStats();
        if (lowTab === 1) meas.update();
      }
      function drawPsi() {
        const Np = 400, th = QM.plot.theme(), xs = new Float64Array(Np + 1), re = new Float64Array(Np + 1), im = new Float64Array(Np + 1), ab = new Float64Array(Np + 1), p = new Float64Array(Np + 1);
        for (let i = 0; i <= Np; i++) { xs[i] = (TAU * i) / Np; const z = psiAt(cur, nMax, xs[i]); re[i] = z.re; im[i] = z.im; p[i] = z.re * z.re + z.im * z.im; ab[i] = Math.sqrt(p[i]); }
        const vMax = carpet.vMax, vl = [{ x: probePhi, color: th.muted }];
        let o;
        if (mode === 'reim') {
          const A = 1.1 * Math.sqrt(vMax);
          o = { yMin: -A, yMax: A, series: [{ x: xs, y: re, color: th.series[0], label: 'Re ψ' }, { x: xs, y: im, color: th.series[1], label: 'Im ψ' }, { x: xs, y: ab, color: th.muted, width: 1, dash: [3, 3], label: '|ψ|' }] };
        } else if (mode === 'dens') o = { yMin: 0, yMax: 1.1 * vMax, series: [{ x: xs, y: p, color: th.accent, fill: true, label: '|ψ|²' }] };
        else o = { yMin: 0, yMax: 1.1 * vMax, series: [{ x: xs, y: p, color: th.ink, width: 1.5, label: '|ψ|², colour = arg ψ' }] };
        if (S.ham === 'pot' && mode !== 'reim') { // scaled potential overlay
          const vs = xs.map((x) => potV(S, x)), lo = Math.min(...vs), hi = Math.max(...vs);
          if (hi - lo > 1e-9) o.series.push({ x: xs, y: vs.map((v) => (0.95 * vMax * (v - lo)) / (hi - lo)), color: th.series[2], width: 1.2, dash: [5, 4], label: 'V(φ) (scaled)' });
        }
        psiGeom = QM.plot.line(cvPsi, Object.assign(o, { xMin: 0, xMax: TAU, piTicks: true, xlabel: 'φ', vlines: vl }));
        if (mode === 'phase') {
          const ctx = cvPsi.getContext('2d'), r = psiGeom; ctx.save(); ctx.globalCompositeOperation = 'destination-over';
          for (let i = 0; i < Np; i++) { const x0 = r.X(xs[i]), y1 = r.Y(p[i]); ctx.fillStyle = QM.plot.phaseColor(Math.atan2(im[i], re[i]), 60); ctx.fillRect(x0, y1, r.X(xs[i + 1]) - x0 + 0.6, r.Y(0) - y1); }
          ctx.restore();
        }
        const z = psiAt(cur, nMax, probePhi);
        probeNote.textContent = `Probe φ = ${QM.fmt(probePhi, 3)} (${QM.fmt(probePhi / Math.PI, 3)}π): ψ = ${C.fmt(z, 3)}, |ψ|² = ${QM.fmt(C.abs2(z), 4)}, arg ψ = ${QM.fmt(C.arg(z), 3)}. Click the plot to move the probe.`;
      }
      function drawRing() {
        const { ctx, w, h } = QM.plot.setup(cvRing), th = QM.plot.theme();
        const cx = w / 2, cy = h / 2, sz = Math.min(w, h), R = sz * 0.27, band = sz * 0.17, K = 240, pMax = carpet.vMax;
        ringGeom = { cx, cy, R };
        ctx.strokeStyle = th.grid; ctx.lineWidth = 1; ctx.beginPath();
        ctx.moveTo(cx - R - band - 8, cy); ctx.lineTo(cx + R + band + 8, cy); ctx.moveTo(cx, cy - R - band - 8); ctx.lineTo(cx, cy + R + band + 8); ctx.stroke();
        ctx.fillStyle = th.muted; ctx.textAlign = 'left'; ctx.textBaseline = 'bottom';
        ctx.fillText('x', cx + R + band + 4, cy - 3); ctx.fillText('y', cx + 4, cy - R - band - 2);
        if (S.ham === 'flux' && Math.abs(S.flux) > 1e-9) { // flux tube through the centre
          ctx.fillStyle = QM.plot.hexA(th.series[3] || th.accent, 0.35); ctx.beginPath(); ctx.arc(cx, cy, R * 0.22, 0, TAU); ctx.fill();
          ctx.fillStyle = th.ink; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText('Φ = ' + QM.fmt(S.flux, 2), cx, cy + R * 0.38);
        }
        for (let k = 0; k < K; k++) {
          const phi = (TAU * (k + 0.5)) / K, z = psiAt(cur, nMax, phi), q = Math.min(1, C.abs2(z) / pMax), thick = 1 + band * q;
          ctx.strokeStyle = mode === 'phase' ? QM.plot.phaseColor(C.arg(z)) : QM.plot.hexA(th.accent, 0.25 + 0.7 * Math.sqrt(q));
          ctx.lineWidth = thick; ctx.beginPath(); ctx.arc(cx, cy, R + thick / 2, -phi - (Math.PI / K) * 1.1, -phi + (Math.PI / K) * 1.1); ctx.stroke();
        }
        ctx.strokeStyle = th.axis; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(cx, cy, R, 0, TAU); ctx.stroke();
        if (S.ham === 'pot') { // potential minimum marker
          let best = 0, bv = Infinity; for (let i = 0; i < 360; i++) { const f = (TAU * i) / 360, v = potV(S, f); if (v < bv) { bv = v; best = f; } }
          ctx.fillStyle = th.series[2]; ctx.beginPath(); ctx.arc(cx + (R - 9) * Math.cos(best), cy - (R - 9) * Math.sin(best), 3.5, 0, TAU); ctx.fill();
        }
        const m = moments(cur, nMax, S.I, S.d);
        if (m.mag > 1e-6) {
          const x = cx + R * m.er, y = cy - R * m.ei, ang = Math.atan2(y - cy, x - cx), Lh = 9;
          ctx.strokeStyle = th.series[1]; ctx.lineWidth = 3; ctx.lineCap = 'round'; ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(x, y);
          ctx.moveTo(x, y); ctx.lineTo(x - Lh * Math.cos(ang - 0.5), y - Lh * Math.sin(ang - 0.5)); ctx.moveTo(x, y); ctx.lineTo(x - Lh * Math.cos(ang + 0.5), y - Lh * Math.sin(ang + 0.5)); ctx.stroke();
        }
        ctx.fillStyle = th.ink; ctx.beginPath(); ctx.arc(cx + R * Math.cos(probePhi), cy - R * Math.sin(probePhi), 4.5, 0, TAU); ctx.fill();
        ctx.beginPath(); ctx.arc(cx, cy, 2.3, 0, TAU); ctx.fill();
      }
      function drawBars() {
        const p = V.probs(cur), ph = [];
        for (let k = 0; k < cur.n; k++) ph.push(Math.atan2(cur.im[k], cur.re[k]));
        QM.plot.bars(cvBars, { labels: basisLabels.slice(), values: Array.from(p), phases: ph, yMax: Math.max(1e-6, ...p) * 1.15, xlabel: 'n', ylabel: '|a_n|²', valueLabels: nMax <= 7 });
      }
      /** energy-level diagram: free → E_n vs n; flux → parabolas E_n(Φ); pot → V(φ) with levels E_k */
      function drawSpec() {
        const th = QM.plot.theme(), N = 2 * nMax + 1, I = S.I, E = sys.eig.values;
        const wE = sys.eig.vectors.map((v) => C.abs2(V.dot(v, a))), eH = M.expect(sys.H, cur).re;
        const pot = S.ham === 'pot';
        specGrid.className = pot ? 'row2' : ''; cvW.style.display = pot ? '' : 'none';
        let Ew = E[Math.min(N - 1, 2)]; wE.forEach((w, k) => { if (w > 1e-3) Ew = Math.max(Ew, E[k]); });
        const dot = (r, x, y, w) => { const ctx = cvLev.getContext('2d'); ctx.fillStyle = QM.plot.hexA(th.accent, 0.35 + 0.65 * Math.min(1, 3 * w)); ctx.strokeStyle = th.ink; ctx.lineWidth = 0.8; ctx.beginPath(); ctx.arc(r.X(x), r.Y(y), 2.5 + 11 * Math.sqrt(w), 0, TAU); ctx.fill(); ctx.stroke(); };
        const hl = [{ y: eH, color: th.series[1], dash: [6, 4] }];
        if (S.ham === 'flux') {
          const xs = new Float64Array(241), lo = -1.5, hi = 1.5, series = [], w = V.probs(a);
          for (let i = 0; i <= 240; i++) xs[i] = lo + ((hi - lo) * i) / 240;
          for (let n = -nMax; n <= nMax; n++) if (Math.abs(n) <= 4 || w[n + nMax] > 1e-4) series.push({ x: xs, y: xs.map((f) => (n - f) ** 2 / (2 * I)), color: w[n + nMax] > 1e-3 ? th.accent : th.axis, width: w[n + nMax] > 1e-3 ? 1.8 : 1 });
          const yMax = 1.15 * Math.max(Ew, 1.2 / I);
          const r = QM.plot.line(cvLev, { xMin: lo, xMax: hi, yMin: 0, yMax, series, legend: false, xlabel: 'Φ (flux quanta)', ylabel: 'E_n(Φ) = (n−Φ)²/2I, dot size ∝ |a_n|²', vlines: [{ x: S.flux, color: th.ink, dash: [3, 3] }], hlines: hl });
          const ctx = cvLev.getContext('2d'); ctx.fillStyle = th.muted; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
          for (let n = -1; n <= 1; n++) ctx.fillText('n=' + fmtN(n), r.X(n), r.Y(0) - 3);
          for (let n = -nMax; n <= nMax; n++) { const En = (n - S.flux) ** 2 / (2 * I); if (En <= yMax) dot(r, S.flux, En, w[n + nMax]); }
        } else if (!pot) {
          const xs = new Float64Array(201), w = V.probs(a);
          for (let i = 0; i <= 200; i++) xs[i] = -nMax - 0.5 + ((2 * nMax + 1) * i) / 200;
          const yMax = 1.15 * Math.max(Ew, 1 / I);
          const r = QM.plot.line(cvLev, { xMin: -nMax - 0.5, xMax: nMax + 0.5, yMin: -0.05 * yMax, yMax, xlabel: 'n', ylabel: 'E_n = n²/2I, dot size ∝ |a_n|²', legend: false, hlines: hl,
            series: [{ x: xs, y: xs.map((x) => (x * x) / (2 * I)), color: th.axis, width: 1, dash: [4, 4] }] });
          for (let n = -nMax; n <= nMax; n++) { const En = (n * n) / (2 * I); if (En <= yMax) dot(r, n, En, w[n + nMax]); }
        } else {
          const xs = new Float64Array(361); for (let i = 0; i <= 360; i++) xs[i] = (TAU * i) / 360;
          const vs = xs.map((x) => potV(S, x)), vLo = Math.min(...vs), vHi = Math.max(...vs);
          const yLo = Math.min(vLo, E[0]), yHi = Math.max(vHi, Ew, E[Math.min(N - 1, 3)]), pad = 0.08 * Math.max(1, yHi - yLo);
          const r = QM.plot.line(cvLev, { xMin: 0, xMax: TAU, yMin: yLo - pad, yMax: yHi + pad, piTicks: true, xlabel: 'φ', ylabel: 'V(φ) and levels E_k (thickness ∝ weight)', legend: false, hlines: hl,
            series: [{ x: xs, y: vs, color: th.series[2], width: 2 }] });
          const ctx = cvLev.getContext('2d'); ctx.save(); ctx.font = '11px system-ui, sans-serif'; ctx.textAlign = 'right'; ctx.textBaseline = 'bottom';
          E.forEach((En, k) => {
            if (En > yHi + pad) return; const y = r.Y(En), w = wE[k];
            ctx.strokeStyle = QM.plot.hexA(th.accent, 0.3 + 0.7 * Math.min(1, 3 * w)); ctx.lineWidth = 1 + 6 * Math.sqrt(w);
            ctx.beginPath(); ctx.moveTo(r.box.x, y); ctx.lineTo(r.box.x + r.box.w, y); ctx.stroke();
            if (w > 0.02) { ctx.fillStyle = th.ink; ctx.fillText(`k=${k}: ${QM.fmt(w, 2)}`, r.box.x + r.box.w - 3, y - 2); }
          });
          ctx.restore();
          const K = Math.min(N, Math.max(8, 1 + wE.reduce((m, w, k) => (w > 1e-3 ? k : m), 0)));
          QM.plot.bars(cvW, { labels: [...Array(K).keys()].map(String), values: wE.slice(0, K), yMax: Math.max(1e-6, ...wE) * 1.15, xlabel: 'level k', ylabel: '|⟨E_k|ψ⟩|²', valueLabels: K <= 10 });
        }
        // harmonic estimate around the potential minimum
        let extra = '';
        if (pot) {
          let fm = 0, vm = Infinity; for (let i = 0; i < 720; i++) { const f = (TAU * i) / 720, v = potV(S, f); if (v < vm) { vm = v; fm = f; } }
          const hh = 1e-3, v2 = (potV(S, fm + hh) - 2 * potV(S, fm) + potV(S, fm - hh)) / (hh * hh);
          if (v2 > 1e-9) { const om = Math.sqrt(v2 / I); extra = ` Harmonic approximation at the minimum φ = ${QM.fmt(fm / Math.PI, 3)}π: ω = √(V''/I) = ${QM.fmt(om, 4)}, E₀ ≈ V_min + ω/2 = ${QM.fmt(vm + om / 2, 4)}, classical period 2π/ω = ${QM.fmt(TAU / om, 4)}.`; }
        }
        const deg = S.ham === 'flux' ? ` Spectrum is periodic in Φ (period 1 flux quantum); at Φ = ½ the levels n and 1−n are degenerate. Ground-state persistent current ∝ ⟨L−Φ⟩/I.` : '';
        specNote.textContent = `E₀ = ${QM.fmt(E[0], 4)}, E₁ − E₀ = ${QM.fmt(E[1] - E[0], 4)}, E₂ − E₁ = ${QM.fmt(E[2] - E[1], 4)}; ⟨H⟩ = ${QM.fmt(eH, 4)} (dashed). Weights are conserved in time.` + deg + extra;
      }
      function drawCarpet() {
        const T = carpet.T;
        carGeom = QM.plot.heatmap(cvCar, { data: carpet.data, nx: carpet.nx, ny: carpet.ny, xMin: 0, xMax: TAU, yMin: 0, yMax: T, vMax: carpet.vMax, piTicks: true, xlabel: 'φ', hline: Math.min(tr.t, T), marker: { x: probePhi, y: Math.min(tr.t, T) } });
        if (S.ham === 'pot' || Math.abs(T - tRev()) > 1e-9) return;
        const ctx = cvCar.getContext('2d'), r = carGeom;
        ctx.save(); ctx.strokeStyle = 'rgba(255,255,255,.55)'; ctx.setLineDash([2, 4]); ctx.lineWidth = 1;
        ctx.beginPath(); ctx.moveTo(r.box.x, r.Y(T / 2)); ctx.lineTo(r.box.x + r.box.w, r.Y(T / 2)); ctx.stroke();
        ctx.setLineDash([]); ctx.fillStyle = 'rgba(255,255,255,.8)'; ctx.textAlign = 'right'; ctx.textBaseline = 'bottom';
        ctx.fillText(S.ham === 'flux' ? 'T/2: ψ(φ+π+2πΦ)' : 'T/2: ψ(φ+π)', r.box.x + r.box.w - 4, r.Y(T / 2) - 2); ctx.restore();
      }
      function updateStats() {
        const m = moments(cur, nMax, S.I, S.d), f = (v) => QM.fmt(v, 4), Hc = M.apply(sys.H, cur), eH = V.dot(cur, Hc).re;
        st.set('L', f(m.L)); st.set('dL', f(m.dL)); st.set('H', f(eH)); st.set('dH', f(Math.sqrt(Math.max(0, V.norm2(Hc) - eH * eH)))); st.set('x', f(m.x)); st.set('y', f(m.y)); st.set('e', f(m.mag));
        st.set('pp', f(m.pp)); st.set('pm', f(m.pm)); st.set('p0', f(m.p0)); st.set('tau', QM.fmt(tr.t / (2 * S.I) / Math.PI, 4) + 'π');
      }

      // ---------- pointer interaction ----------
      const local = (cv, ev) => { const r = cv.getBoundingClientRect(); return { x: ev.clientX - r.left, y: ev.clientY - r.top }; };
      const clamp01 = (v) => Math.max(0, Math.min(1, v));
      cvPsi.addEventListener('pointerdown', (ev) => { if (!psiGeom) return; const { x } = local(cvPsi, ev), b = psiGeom.box; probePhi = TAU * clamp01((x - b.x) / b.w); redraw(); });
      cvRing.addEventListener('pointerdown', (ev) => { if (!ringGeom) return; const { x, y } = local(cvRing, ev); probePhi = (Math.atan2(-(y - ringGeom.cy), x - ringGeom.cx) + TAU) % TAU; redraw(); });
      cvCar.addEventListener('pointerdown', (ev) => {
        if (!carGeom) return; const { x, y } = local(cvCar, ev), b = carGeom.box;
        probePhi = TAU * clamp01((x - b.x) / b.w); tr.set(carpet.T * clamp01((b.y + b.h - y) / b.h));
      });

      // ---------- start ----------
      showGroup(); showHam(); carNote();
      computeCarpet();
      renderMatrix();
      meas.refresh();
      modeTabs.select(Math.max(0, modeKeys.indexOf(mode)));
      low.select(lowTab);
      redraw(); drawDir();
      const stopResize = QM.ui.onResize(views, () => { redraw(); drawDir(); });
      if (S.play) tr.play();
      return {
        destroy() { tr.destroy(); stopAnim(); stopResize(); },
        redraw() { redraw(); drawDir(); },
      };
    },
  });

  // ---------- self-tests ----------
  const close = (x, y, tol = 1e-9) => Math.abs(x - y) < tol;
  const St = (o) => Object.assign({}, DEFAULTS, o);
  QM.addTest('ring: (x̂ − 3iŷ)|0⟩ gives P(ℓ>0) = 1/5, P(ℓ<0) = 4/5', () => {
    const m = moments(buildState(St({ builder: 'xiy', cRe: -3, cIm: 0 }), 6), 6, 0.5, 1);
    return { pass: close(m.pp, 0.2) && close(m.pm, 0.8), detail: `P(ℓ>0) = ${QM.fmt(m.pp, 6)}, P(ℓ<0) = ${QM.fmt(m.pm, 6)}` };
  });
  QM.addTest('ring: full revival at τ = 2π, half revival = shift by π (parity)', () => {
    const S = St({ builder: 'gauss', phi0: 1, sigma: 0.4, n0: 2.3 }), a = buildState(S, 8);
    let err = 0, errH = 0;
    [0.5, 0.73].forEach((I) => {
      const T = 4 * Math.PI * I, full = evolveC(a, 8, I, T), half = evolveC(a, 8, I, T / 2);
      const P = M.apply(makeOps(8, 1, I).P, a);
      for (let k = 0; k < a.n; k++) { err = Math.max(err, Math.hypot(full.re[k] - a.re[k], full.im[k] - a.im[k])); errH = Math.max(errH, Math.hypot(half.re[k] - P.re[k], half.im[k] - P.im[k])); }
    });
    return { pass: err < 1e-9 && errH < 1e-9, detail: `max |a(T) − a(0)| = ${err.toExponential(1)}, max |a(T/2) − P̂a| = ${errH.toExponential(1)}` };
  });
  QM.addTest('ring: Â† on |2⟩ raises ⟨L⟩ by 1, norm d', () => {
    const d = 1.7, ops = makeOps(6, d, 0.5), w = M.apply(ops.Ap, V.basis(13, 2 + 6)), nrm = Math.sqrt(V.norm2(w));
    const m = moments(V.normalize(w), 6, 0.5, d);
    return { pass: close(m.L, 3) && close(nrm, d), detail: `⟨L⟩ = ${QM.fmt(m.L, 6)}, ‖Â†|2⟩‖ = ${QM.fmt(nrm, 6)} (d = ${d})` };
  });
  QM.addTest('ring: truncated [L̂, x̂] = iŷ, [L̂, Â†] = Â†', () => {
    const o = makeOps(5, 1.3, 0.5);
    const e1 = M.maxAbs(M.sub(M.comm(o.L, o.x), M.scale(o.y, { re: 0, im: 1 }))), e2 = M.maxAbs(M.sub(M.comm(o.L, o.Ap), o.Ap));
    return { pass: e1 < 1e-12 && e2 < 1e-12, detail: `errors ${e1.toExponential(1)}, ${e2.toExponential(1)}` };
  });
  QM.addTest('ring: ⟨e^{iφ}⟩ from coefficients matches ∫|ψ|² e^{iφ} dφ', () => {
    const a = buildState(St({ phi0: 2, sigma: 0.7, n0: 1 }), 6), m = moments(a, 6, 0.5, 1);
    let re = 0, im = 0; const Nq = 2000;
    for (let i = 0; i < Nq; i++) { const phi = (TAU * i) / Nq, z = psiAt(a, 6, phi), p = z.re * z.re + z.im * z.im; re += p * Math.cos(phi); im += p * Math.sin(phi); }
    re *= TAU / Nq; im *= TAU / Nq;
    return { pass: close(re, m.er, 1e-8) && close(im, m.ei, 1e-8), detail: `quadrature ${QM.fmt(re, 5)} + ${QM.fmt(im, 5)}i vs ${QM.fmt(m.er, 5)} + ${QM.fmt(m.ei, 5)}i` };
  });
  QM.addTest('ring: Dirichlet kernel D_N(0) = (2N+1)/2π and equals the Fourier sum', () => {
    let ok = true, det = [];
    [0, 3, 10, 57].forEach((N) => {
      const d0 = dirichlet(0, N), d2 = dirichlet(4 * Math.PI, N); let s = 0; for (let n = -N; n <= N; n++) s += Math.cos(n * 0.7);
      ok = ok && close(d0, (2 * N + 1) / TAU) && close(d2, d0) && close(dirichlet(0.7, N), s / TAU, 1e-9);
      det.push(`N=${N}: ${QM.fmt(d0, 4)}`);
    });
    return { pass: ok, detail: det.join(', ') };
  });
  QM.addTest('ring: V = 0 potential (numerical eigh path) evolves exactly like the analytic phases', () => {
    const S = St({ ham: 'pot', V1: 0, V2: 0, I: 0.7, builder: 'gauss', phi0: 1.1, sigma: 0.5, n0: 1.6 }), sys = makeSystem(S, 8), a = buildState(S, 8, sys);
    let err = 0;
    [0.37, 2.9, 11.3].forEach((t) => { const c1 = evolveSys(sys, a, 8, 0.7, t), c2 = evolveC(a, 8, 0.7, t); for (let k = 0; k < a.n; k++) err = Math.max(err, Math.hypot(c1.re[k] - c2.re[k], c1.im[k] - c2.im[k])); });
    const eErr = Math.max(...sys.eig.values.map((E, i) => Math.abs(E - makeSystem(St({ I: 0.7 }), 8).eig.values[i])));
    return { pass: !sys.diag && err < 1e-9 && eErr < 1e-10, detail: `max |a_num(t) − a_ana(t)| = ${err.toExponential(1)}, max level error ${eErr.toExponential(1)}` };
  });
  QM.addTest('ring: flux Φ = ½ → E_n = E_{1−n} (degenerate pairs), numerical spectrum pairs up', () => {
    const nMax = 6, S = St({ ham: 'flux', flux: 0.5, I: 0.8 }), E = (n) => (n - 0.5) ** 2 / 1.6;
    let ok = true; for (let n = -5; n <= 6; n++) ok = ok && close(E(n), E(1 - n));
    const num = makeSystem(S, nMax, true).eig.values; let pe = 0;
    for (let j = 0; j + 1 < 12; j += 2) pe = Math.max(pe, Math.abs(num[j] - num[j + 1]));
    const gap = num[2] - num[1];
    return { pass: ok && pe < 1e-10 && gap > 0.5, detail: `levels ${num.slice(0, 6).map((v) => QM.fmt(v, 4)).join(', ')}…; max pair splitting ${pe.toExponential(1)}` };
  });
  QM.addTest('ring: pendulum V₁cos φ ground state ≈ ω/2 − V₁ (harmonic, ω = √(V₁/I))', () => {
    const I = 1, V1 = 50, sys = makeSystem(St({ ham: 'pot', V1, phi1: 0, V2: 0, I }), 20), om = Math.sqrt(V1 / I), Eh = om / 2 - V1;
    const E0 = sys.eig.values[0], gap = sys.eig.values[1] - E0;
    return { pass: Math.abs(E0 - Eh) < 0.02 * om && Math.abs(gap - om) < 0.05 * om, detail: `E₀ = ${QM.fmt(E0, 5)}, ω/2 − V₁ = ${QM.fmt(Eh, 5)} (quartic correction ≈ −1/32I), E₁−E₀ = ${QM.fmt(gap, 4)} vs ω = ${QM.fmt(om, 4)}` };
  });
  QM.addTest('ring: potential evolution is unitary, U(−t)U(t) = 1, ⟨H⟩ conserved', () => {
    const S = St({ ham: 'pot', V1: -3, phi1: 1, V2: 2, I: 0.6, builder: 'cat', phi0: 0.5, sigma: 0.4, n0: 1, catSep: 2, catPhase: 1 });
    const sys = makeSystem(S, 10), a = buildState(S, 10, sys), c = evolveSys(sys, a, 10, 0.6, 3.7), b = evolveSys(sys, c, 10, 0.6, -3.7);
    let err = 0; for (let k = 0; k < a.n; k++) err = Math.max(err, Math.hypot(b.re[k] - a.re[k], b.im[k] - a.im[k]));
    const h0 = M.expect(sys.H, a).re, h1 = M.expect(sys.H, c).re;
    return { pass: err < 1e-9 && close(V.norm2(c), 1, 1e-10) && close(h0, h1, 1e-9) && M.isHermitian(sys.H), detail: `round-trip error ${err.toExponential(1)}, ⟨H⟩ ${QM.fmt(h0, 6)} → ${QM.fmt(h1, 6)}` };
  });
  QM.addTest('ring: top-hat arc state – analytic sinc coefficients put ≳ 95 % of |ψ|² on the arc', () => {
    const w = Math.PI / 2, phi0 = 2, a = buildState(St({ builder: 'box', phi0, boxW: w, n0: 0 }), 20);
    let pin = 0; const Nq = 4000;
    for (let i = 0; i < Nq; i++) { const f = (TAU * i) / Nq, z = psiAt(a, 20, f); if (Math.abs(((f - phi0 + 3 * Math.PI) % TAU) - Math.PI) < w / 2) pin += C.abs2(z) * TAU / Nq; }
    return { pass: pin > 0.95, detail: `P(on arc) = ${QM.fmt(pin, 4)} (nMax = 20)` };
  });
})();
