/* Engine 'wave1d' – general 1D wavefunction sandbox (x and p representation, split-step evolution).
   ħ = 1.  Grid: N points on [−L/2, L/2) (periodic), momentum grid p = 2πk/L via FFT.
   ψ̃(p) = (2π)^{-1/2} ∫ψ(x)e^{−ipx}dx, so ∫|ψ̃|²dp = ∫|ψ|²dx = 1.
   Propagation: split-step Fourier for smooth V; potentials with jumps / δ / walls use the finite-difference H with
   Crank–Nicolson (Dirichlet at the grid ends or walls); the infinite box uses the exact sine-series evolution.
   Eigenstates of the current V: tridiagonal FD H (Sturm bisection + inverse iteration) on the simulation grid,
   for smooth V refined by Rayleigh–Ritz with the exact (FFT) kinetic energy.

   Preset schema (every key optional, defaults in DEF below):
   { pot: potential V(x)
       'free'
       'harmonic'  ½mω²x²                          omega
       'delta'     −λδ(x)                          lambda   (one grid site of depth λ/Δx)
       'well'      −V0 for |x| < a/2               V0, a
       'box'       infinite well |x| < wBox/2      wBox     (walls of height 1000)
       'barrier'   Vb for |x| < bw/2               Vb, bw   (transmission / reflection shown)
       'step'      Vs for x > 0                    Vs
       'double'    dwA (x² − dwB²)²                dwA, dwB
       'linear'    F x  (floor: wall for x < 0 → quantum bouncer)   F, floor
       'cos'       Vc sin²(πx/per)  (per rounded so that L/per is an integer)   Vc, per
     m: 1
     state: initial state ψ(t=0)
       'gauss'    ψ(p) ∝ exp(−(p−p0)²/(4Δ²))·e^{−ipx0}     p0, Delta, x0
       'delta'    ψ = √κ e^{−κ|x|}, κ = mλ                  lambda
       'coherent' U|0⟩ = e^{α(â−â†)}|0⟩ = |β = −α⟩         aRe, aIm (α), omega; x0 = 1/√(mω)
       'hoN'      HO eigenstate n (Hermite function) displaced by x0       nHO, x0, omega
       'cat'      G(x−x0−sep/2) + e^{iφ} G(x−x0+sep/2), both with momentum p0, width 1/(2Δ)   sep, phi, p0, Delta, x0
       'box'      sin(nπ(x+w/2)/w) inside |x| < wBox/2                   nBox, wBox
       'ground'   ground state of the current V (numerical)
       'eigen'    eigenstate n1 of the current V (numerical)            n1
       'eigen2'   (|n1⟩ + e^{iφ}|n2⟩)/√2 of the current V                n1, n2, phi
       'exp'      e^{−κ|x−x0|} e^{ip0x}                                   kappa, x0, p0
     L: 40, N: 1024 (power of 2), xView: 8, pView: 8     box length / plotted half-ranges
     xStar: 0, pStar: 0                                   tail thresholds: P(x > x*), P(p > p*)
     tMax: 10 (auto: harmonic 2 periods, box revival time, eigen2 2 beat periods)
     showReIm, showPReIm, analytic: true (analytic overlays) } */
(function () {
  'use strict';
  const QM = window.QM;
  const { el } = QM.ui;
  const SQ2 = Math.SQRT2, DT = 0.01, WALL = 1000, NLEV = 16;

  const DEF = { pot: 'free', m: 1, omega: 1, lambda: 1, V0: 5, a: 2, wBox: 6, Vb: 2, bw: 1, Vs: 2, dwA: 0.25, dwB: 2, F: 0.5, floor: false, Vc: 4, per: 2,
    state: 'gauss', p0: 2, Delta: 1, x0: 0, aRe: 1, aIm: 0, nHO: 1, sep: 6, phi: 0, nBox: 1, n1: 0, n2: 1, kappa: 1,
    L: 40, N: 1024, xView: 8, pView: 8, xStar: 0, pStar: 0, tMax: 10, showReIm: false, showPReIm: false, analytic: true };
  const POT_KEYS = ['pot', 'm', 'omega', 'lambda', 'V0', 'a', 'wBox', 'Vb', 'bw', 'Vs', 'dwA', 'dwB', 'F', 'floor', 'Vc', 'per', 'L', 'N'];
  const BOUND = ['delta', 'ground', 'eigen', 'eigen2', 'box', 'hoN'];
  const CONTINUUM = (P) => ['free', 'barrier', 'step'].includes(P.pot) || (P.pot === 'linear' && !P.floor);

  // ---------- math helpers ----------
  function erfc(x) { // Numerical Recipes erfcc, |rel err| < 1.2e-7
    const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
    const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 +
      t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? r : 2 - r;
  }
  const erf = (x) => 1 - erfc(x);

  function makeGrid(N, L) {
    const dx = L / N, dp = (2 * Math.PI) / L, x = new Float64Array(N), k = new Float64Array(N);
    for (let j = 0; j < N; j++) { x[j] = -L / 2 + j * dx; k[j] = (j < N / 2 ? j : j - N) * dp; }
    return { N, L, dx, dp, x, k };
  }

  function potential(g, P) {
    const V = new Float64Array(g.N), per = g.L / Math.max(1, Math.round(g.L / P.per));
    for (let j = 0; j < g.N; j++) {
      const x = g.x[j];
      switch (P.pot) {
        case 'harmonic': V[j] = 0.5 * P.m * P.omega * P.omega * x * x; break;
        case 'well': V[j] = Math.abs(x) < P.a / 2 ? -P.V0 : 0; break;
        case 'box': V[j] = Math.abs(x) < P.wBox / 2 ? 0 : WALL; break;
        case 'barrier': V[j] = Math.abs(x) < P.bw / 2 ? P.Vb : 0; break;
        case 'step': V[j] = x >= -1e-12 ? P.Vs : 0; break;
        case 'double': V[j] = P.dwA * (x * x - P.dwB * P.dwB) ** 2; break;
        case 'linear': V[j] = P.floor && x <= 1e-12 ? WALL : P.F * x; break;
        case 'cos': V[j] = P.Vc * Math.sin((Math.PI * x) / per) ** 2; break;
      }
    }
    if (P.pot === 'delta') V[g.N / 2] = -P.lambda / g.dx; // x = 0 is grid point N/2
    // Sharp potentials (jumps, δ, walls) defeat the split-step splitting (spurious high-k energy), so they use the
    // finite-difference H on the index range inner = {i0, i1} (Dirichlet outside): Crank–Nicolson, or exact sine series in the box.
    if (['box', 'delta', 'well', 'barrier', 'step'].includes(P.pot) || (P.pot === 'linear' && P.floor)) {
      let i0 = -1, i1 = -1; for (let j = 0; j < g.N; j++) if (V[j] < WALL) { if (i0 < 0) i0 = j; i1 = j; }
      if (i1 > i0) V.inner = { i0, i1, box: P.pot === 'box' };
    }
    return V;
  }
  function boxEnergy(g, V, m, n) { const W = (V.inner.i1 - V.inner.i0 + 2) * g.dx; return ((n + 1) * Math.PI / W) ** 2 / (2 * m); }

  /** propagator for the current V: exact sine-series evolution in the box, Crank–Nicolson with walls, split-step otherwise */
  function propagator(g, V, m) { return !V.inner ? stepper(g, V, m) : V.inner.box ? boxProp(g, V, m) : cnProp(g, V, m); }
  function boxProp(g, V, m) { // ψ = Σ c_n φ_n, φ_n = discrete sines (exact eigenvectors of the FD box), E_n = ideal n²π²/2mW²
    const { i0, i1 } = V.inner, M = i1 - i0 + 1, c = Math.sqrt(2 / ((M + 1) * g.dx)), S = new Float64Array(M * M), cr = new Float64Array(M), ci = new Float64Array(M);
    for (let n = 0; n < M; n++) for (let j = 0; j < M; j++) S[n * M + j] = c * Math.sin(((n + 1) * Math.PI * (j + 1)) / (M + 1));
    return function (re, im, dt, k) {
      const t = dt * k;
      for (let n = 0; n < M; n++) {
        let a = 0, b = 0; for (let j = 0; j < M; j++) { a += S[n * M + j] * re[i0 + j]; b += S[n * M + j] * im[i0 + j]; }
        const ph = -boxEnergy(g, V, m, n) * t, cs = Math.cos(ph), sn = Math.sin(ph); a *= g.dx; b *= g.dx;
        cr[n] = a * cs - b * sn; ci[n] = a * sn + b * cs;
      }
      for (let j = 0; j < M; j++) { let a = 0, b = 0; for (let n = 0; n < M; n++) { a += S[n * M + j] * cr[n]; b += S[n * M + j] * ci[n]; } re[i0 + j] = a; im[i0 + j] = b; }
    };
  }
  function cnProp(g, V, m) { // Crank–Nicolson (1 + iτH)ψ' = (1 − iτH)ψ, τ = dt/2, FD H on the inner range, Dirichlet outside
    const { i0, i1 } = V.inner, M = i1 - i0 + 1, h = 1 / (m * g.dx * g.dx), e = -h / 2, d = Float64Array.from(V.subarray(i0, i1 + 1), (v) => v + h);
    const cpr = new Float64Array(M), cpi = new Float64Array(M), ivr = new Float64Array(M), ivi = new Float64Array(M), yr = new Float64Array(M), yi = new Float64Array(M);
    let cdt = NaN, tau = 0;
    function prep(dt) { // forward elimination of A = tridiag(a, 1 + iτd, a), a = iτe
      cdt = dt; tau = dt / 2; const ai = tau * e;
      for (let j = 0; j < M; j++) {
        let dr = 1, di = tau * d[j]; if (j) { dr -= -ai * cpi[j - 1]; di -= ai * cpr[j - 1]; } // α − a c'
        const q = dr * dr + di * di; ivr[j] = dr / q; ivi[j] = -di / q;
        cpr[j] = -ai * ivi[j]; cpi[j] = ai * ivr[j]; // a / (α − a c')
      }
    }
    return function (re, im, dt, k) {
      const sub = 2, ds = dt / sub; if (ds !== cdt) prep(ds); const ai = tau * e;
      for (let s = 0; s < k * sub; s++) {
        for (let j = 0; j < M; j++) {
          const x = i0 + j, l = j ? x - 1 : -1, r = j < M - 1 ? x + 1 : -1;
          const hr = d[j] * re[x] + e * ((l < 0 ? 0 : re[l]) + (r < 0 ? 0 : re[r])), hi = d[j] * im[x] + e * ((l < 0 ? 0 : im[l]) + (r < 0 ? 0 : im[r]));
          let rr = re[x] + tau * hi, ri = im[x] - tau * hr;
          if (j) { rr -= -ai * yi[j - 1]; ri -= ai * yr[j - 1]; }
          yr[j] = rr * ivr[j] - ri * ivi[j]; yi[j] = rr * ivi[j] + ri * ivr[j];
        }
        for (let j = M - 1; j >= 0; j--) {
          const x = i0 + j; let a = yr[j], b = yi[j];
          if (j < M - 1) { a -= cpr[j] * re[x + 1] - cpi[j] * im[x + 1]; b -= cpr[j] * im[x + 1] + cpi[j] * re[x + 1]; }
          re[x] = a; im[x] = b;
        }
      }
    };
  }

  /** Gaussian-type initial states: centre xm, pm and position width s0 (minimum uncertainty, unchirped) */
  function gaussInfo(P) {
    if (P.state === 'gauss') return { xm: P.x0, pm: P.p0, s0: 1 / (2 * P.Delta) };
    if (P.state === 'coherent') { const x0 = 1 / Math.sqrt(P.m * P.omega); return { xm: -SQ2 * x0 * P.aRe, pm: (-SQ2 * P.aIm) / x0, s0: x0 / SQ2 }; }
    return null;
  }

  function normalize(g, re, im) {
    let s = 0; for (let j = 0; j < g.N; j++) s += re[j] * re[j] + im[j] * im[j];
    s = Math.sqrt(s * g.dx); if (s > 0) for (let j = 0; j < g.N; j++) { re[j] /= s; im[j] /= s; }
  }

  function phaseMul(re, im, cr, ci) {
    for (let j = 0; j < re.length; j++) { const a = re[j], b = im[j]; re[j] = a * cr[j] - b * ci[j]; im[j] = a * ci[j] + b * cr[j]; }
  }

  /** split-step Fourier propagator e^{−iVdt/2} e^{−ip²dt/2m} e^{−iVdt/2}; step(re, im, dt, nSteps) */
  function stepper(g, V, m) {
    const N = g.N, vr = new Float64Array(N), vi = new Float64Array(N), kr = new Float64Array(N), ki = new Float64Array(N);
    let cdt = NaN;
    function prep(dt) {
      cdt = dt;
      for (let j = 0; j < N; j++) {
        const a = (-V[j] * dt) / 2, b = (-g.k[j] * g.k[j] * dt) / (2 * m);
        vr[j] = Math.cos(a); vi[j] = Math.sin(a); kr[j] = Math.cos(b); ki[j] = Math.sin(b);
      }
    }
    return function step(re, im, dt, n) {
      if (dt !== cdt) prep(dt);
      for (let s = 0; s < n; s++) { phaseMul(re, im, vr, vi); QM.fft(re, im, false); phaseMul(re, im, kr, ki); QM.fft(re, im, true); phaseMul(re, im, vr, vi); }
    };
  }

  /** Hψ for real ψ with the exact (spectral) kinetic energy */
  function applyH(g, V, m, v) {
    const re = Float64Array.from(v), im = new Float64Array(g.N);
    QM.fft(re, im, false);
    for (let j = 0; j < g.N; j++) { const f = (g.k[j] * g.k[j]) / (2 * m); re[j] *= f; im[j] *= f; }
    QM.fft(re, im, true);
    for (let j = 0; j < g.N; j++) re[j] += V[j] * v[j];
    return re;
  }
  const dotR = (a, b, dx) => { let s = 0; for (let j = 0; j < a.length; j++) s += a[j] * b[j]; return s * dx; };
  function fixSign(v) { // rightmost significant lobe positive (Hermite-function convention)
    let mx = 0; for (const a of v) mx = Math.max(mx, Math.abs(a));
    for (let j = v.length - 1; j >= 0; j--) if (Math.abs(v[j]) > 0.05 * mx) { if (v[j] < 0) for (let i = 0; i < v.length; i++) v[i] = -v[i]; break; }
    return v;
  }

  /** lowest K eigenpairs {E[], vec[]} (vec normalised, ∫v²dx = 1). FD tridiagonal H → Sturm bisection + inverse iteration,
      then Rayleigh–Ritz in the span of K+8 FD vectors with the spectral kinetic energy. */
  function eigenstates(g, V, m, K = NLEV, raw = false) {
    const I = V.inner || { i0: 0, i1: g.N - 1 }, N = I.i1 - I.i0 + 1, h = 1 / (m * g.dx * g.dx), e = -h / 2, M = Math.min(N, K + 8);
    const d = Float64Array.from(V.subarray(I.i0, I.i1 + 1), (v) => v + h);
    let lo = Infinity, hi = -Infinity; for (const v of d) { lo = Math.min(lo, v); hi = Math.max(hi, v); } lo -= h; hi += h;
    const count = (x) => { let c = 0, q = 1; for (let j = 0; j < N; j++) { q = d[j] - x - (j ? (e * e) / q : 0); if (q === 0) q = -1e-300; if (q < 0) c++; } return c; };
    const solve = (lam, r) => { // (T − λ) y = r, Thomas algorithm
      const c = new Float64Array(N), y = new Float64Array(N);
      for (let j = 0; j < N; j++) {
        let piv = d[j] - lam - (j ? e * c[j - 1] : 0); if (Math.abs(piv) < 1e-300) piv = 1e-300;
        c[j] = e / piv; y[j] = (r[j] - (j ? e * y[j - 1] : 0)) / piv;
      }
      for (let j = N - 2; j >= 0; j--) y[j] -= c[j] * y[j + 1];
      return y;
    };
    const B = [];
    for (let n = 0; n < M; n++) {
      let a = lo, b = hi;
      for (let it = 0; it < 200 && b - a > 1e-13 * Math.max(1, Math.abs(a), Math.abs(b)); it++) { const c = 0.5 * (a + b); if (count(c) > n) b = c; else a = c; }
      const lam = 0.5 * (a + b);
      let v = Float64Array.from({ length: N }, (_, j) => 1 + 0.5 * Math.sin(0.37 * j + n));
      for (let it = 0; it < 4; it++) {
        v = solve(lam, v);
        for (const u of B) { const s = dotR(u, v, g.dx); for (let j = 0; j < N; j++) v[j] -= s * u[j]; }
        const s = Math.sqrt(dotR(v, v, g.dx)); for (let j = 0; j < N; j++) v[j] /= s;
      }
      B.push(v);
    }
    if (V.inner) { // walls: FD eigenvectors embedded in the full grid (box: ideal energies, matching the sine propagator)
      const E = [], vec = [];
      for (let n = 0; n < Math.min(K, M); n++) {
        let lam = 0; for (let j = 0; j < N; j++) { const l = j ? B[n][j - 1] : 0, r = j < N - 1 ? B[n][j + 1] : 0; lam += B[n][j] * (d[j] * B[n][j] + e * (l + r)); }
        E.push(V.inner.box && !raw ? boxEnergy(g, V, m, n) : lam * g.dx);
        const u = new Float64Array(g.N); u.set(B[n], I.i0); vec.push(fixSign(u));
      }
      return { E, vec };
    }
    // Rayleigh–Ritz with the exact kinetic energy
    const HB = B.map((v) => applyH(g, V, m, v)), H = { n: M, re: new Float64Array(M * M), im: new Float64Array(M * M) };
    for (let i = 0; i < M; i++) for (let j = 0; j <= i; j++) H.re[i * M + j] = H.re[j * M + i] = 0.5 * (dotR(B[i], HB[j], g.dx) + dotR(B[j], HB[i], g.dx));
    const R = QM.eigh(H), E = R.values.slice(0, K);
    const vec = R.vectors.slice(0, K).map((c) => {
      const u = new Float64Array(N); for (let i = 0; i < M; i++) { const w = c.re[i]; for (let j = 0; j < N; j++) u[j] += w * B[i][j]; }
      const s = Math.sqrt(dotR(u, u, g.dx)); for (let j = 0; j < N; j++) u[j] /= s;
      return fixSign(u);
    });
    return { E, vec };
  }

  /** HO eigenfunctions ψ_0..ψ_nMax(x − shift), x0 = 1/√(mω) (Hermite-function recurrence) */
  function hermite(xs, x0, nMax, shift = 0) {
    const N = xs.length, out = [], c0 = Math.pow(Math.PI, -0.25) / Math.sqrt(x0);
    let prev = new Float64Array(N), cur = new Float64Array(N);
    for (let j = 0; j < N; j++) { const xi = (xs[j] - shift) / x0; cur[j] = c0 * Math.exp((-xi * xi) / 2); }
    out.push(cur);
    for (let n = 0; n < nMax; n++) {
      const nxt = new Float64Array(N), a = Math.sqrt(2 / (n + 1)), b = Math.sqrt(n / (n + 1));
      for (let j = 0; j < N; j++) nxt[j] = a * ((xs[j] - shift) / x0) * cur[j] - b * prev[j];
      out.push(nxt); prev = cur; cur = nxt;
    }
    return out;
  }

  function initState(g, P, eig, V) {
    const s = initRaw(g, P, eig);
    if (V && V.inner) {
      for (let j = 0; j < g.N; j++) if (j < V.inner.i0 || j > V.inner.i1) s.re[j] = s.im[j] = 0;
      if (dotR(s.re, s.re, 1) + dotR(s.im, s.im, 1) < 1e-20) s.re.set(eig.vec[0]); // state entirely behind the wall
      normalize(g, s.re, s.im);
    }
    return s;
  }
  function initRaw(g, P, eig) {
    const N = g.N, G = gaussInfo(P), re = new Float64Array(N), im = new Float64Array(N), x = g.x;
    const setAmp = (j, amp, ph) => { re[j] = amp * Math.cos(ph); im[j] = amp * Math.sin(ph); };
    const st = P.state, nMax = eig ? eig.vec.length - 1 : 0, clampN = (n) => Math.max(0, Math.min(nMax, Math.round(n)));
    if (st === 'ground' || st === 'eigen') re.set(eig.vec[st === 'ground' ? 0 : clampN(P.n1)]);
    else if (st === 'eigen2') {
      const u = eig.vec[clampN(P.n1)], v = eig.vec[clampN(P.n2)], c = Math.cos(P.phi), s = Math.sin(P.phi);
      for (let j = 0; j < N; j++) { re[j] = u[j] + c * v[j]; im[j] = s * v[j]; }
    } else if (st === 'hoN') re.set(hermite(x, 1 / Math.sqrt(P.m * P.omega), Math.round(P.nHO), P.x0)[Math.round(P.nHO)]);
    else for (let j = 0; j < N; j++) {
      if (G) { const dd = x[j] - G.xm; setAmp(j, Math.exp(-(dd * dd) / (4 * G.s0 * G.s0)), G.pm * dd); }
      else if (st === 'cat') {
        const s0 = 1 / (2 * P.Delta), d1 = x[j] - P.x0 - P.sep / 2, d2 = x[j] - P.x0 + P.sep / 2;
        const g1 = Math.exp(-(d1 * d1) / (4 * s0 * s0)), g2 = Math.exp(-(d2 * d2) / (4 * s0 * s0)), ph = P.p0 * (x[j] - P.x0);
        re[j] = g1 * Math.cos(ph) + g2 * Math.cos(ph + P.phi); im[j] = g1 * Math.sin(ph) + g2 * Math.sin(ph + P.phi);
      } else if (st === 'box') re[j] = Math.abs(x[j]) < P.wBox / 2 ? Math.sin((Math.round(P.nBox) * Math.PI * (x[j] + P.wBox / 2)) / P.wBox) : 0;
      else if (st === 'exp') setAmp(j, Math.exp(-P.kappa * Math.abs(x[j] - P.x0)), P.p0 * (x[j] - P.x0));
      else re[j] = Math.exp(-P.m * P.lambda * Math.abs(x[j])); // δ-well bound state
    }
    normalize(g, re, im);
    return { re, im };
  }

  /** momentum representation, sorted by p: {p, re, im, rho} with ∫rho dp = 1 */
  function toP(g, re, im) {
    const N = g.N, fr = Float64Array.from(re), fi = Float64Array.from(im);
    QM.fft(fr, fi, false);
    const c = g.dx / Math.sqrt(2 * Math.PI), p = new Float64Array(N), pr = new Float64Array(N), pi = new Float64Array(N), rho = new Float64Array(N);
    for (let s = 0; s < N; s++) {
      const j = (s + N / 2) % N, pk = g.k[j], ph = -pk * g.x[0], cr = Math.cos(ph), ci = Math.sin(ph);
      p[s] = pk; pr[s] = c * (fr[j] * cr - fi[j] * ci); pi[s] = c * (fr[j] * ci + fi[j] * cr);
      rho[s] = pr[s] * pr[s] + pi[s] * pi[s];
    }
    return { p, re: pr, im: pi, rho };
  }

  function moments(xs, rho, d) {
    let n = 0, m1 = 0, m2 = 0;
    for (let i = 0; i < xs.length; i++) { const w = rho[i] * d; n += w; m1 += w * xs[i]; m2 += w * xs[i] * xs[i]; }
    m1 /= n; m2 /= n;
    return { norm: n, mean: m1, sd: Math.sqrt(Math.max(0, m2 - m1 * m1)), m2 };
  }
  /** P(var > t), each sample representing the bin [x − d/2, x + d/2] */
  function tail(xs, rho, d, t) {
    let s = 0;
    for (let i = 0; i < xs.length; i++) { const hi = xs[i] + d / 2; if (hi <= t) continue; s += rho[i] * Math.min(d, hi - t); }
    return s;
  }

  function analyse(g, re, im, V, m) {
    const rx = new Float64Array(g.N);
    let ev = 0;
    for (let j = 0; j < g.N; j++) { rx[j] = re[j] * re[j] + im[j] * im[j]; ev += rx[j] * V[j] * g.dx; }
    const P = toP(g, re, im), mx = moments(g.x, rx, g.dx), mp = moments(P.p, P.rho, g.dp);
    return { rx, P, mx, mp, E: mp.m2 * mp.norm / (2 * m) + ev };
  }

  /** |⟨n|ψ⟩|² for HO eigenstates n = 0..nMax */
  function fock(g, re, im, m, omega, nMax) {
    return Float64Array.from(hermite(g.x, 1 / Math.sqrt(m * omega), nMax), (h) => { const a = dotR(h, re, g.dx), b = dotR(h, im, g.dx); return a * a + b * b; });
  }

  // ---------- scattering: edges of the scatterer on the grid and plane-wave transmission ----------
  function scatterEdges(g, P) {
    if (P.pot === 'step') return { xl: -g.dx / 2, xr: -g.dx / 2 };
    const w = { delta: 0.5, barrier: P.bw, well: P.a }[P.pot]; if (w === undefined) return null; // grid sites with |x| < w/2 carry V
    let i0 = -1, i1 = -1; for (let j = 0; j < g.N; j++) if (Math.abs(g.x[j]) < (P.pot === 'delta' ? g.dx * w : w / 2)) { if (i0 < 0) i0 = j; i1 = j; }
    return i0 < 0 ? null : { xl: g.x[i0] - g.dx / 2, xr: g.x[i1] + g.dx / 2 };
  }
  /** plane-wave transmission probability T(p), p > 0, incident from the left */
  function transCoef(P, e, p) {
    const m = P.m, E = (p * p) / (2 * m);
    if (P.pot === 'delta') return 1 / (1 + (m * P.lambda / p) ** 2);
    if (P.pot === 'step') { if (E <= P.Vs) return 0; const q = Math.sqrt(2 * m * (E - P.Vs)); return (4 * p * q) / (p + q) ** 2; }
    const U = P.pot === 'barrier' ? P.Vb : -P.V0, w = e.xr - e.xl, D = E - U;
    if (Math.abs(D) < 1e-9) return 1 / (1 + (m * w * w * U) / 2);
    const s = D > 0 ? Math.sin(Math.sqrt(2 * m * D) * w) : Math.sinh(Math.sqrt(-2 * m * D) * w);
    return 1 / (1 + (U * U * s * s) / (4 * E * Math.abs(D)));
  }
  function sideProbs(g, rx, e) {
    let T = 0, R = 0; for (let j = 0; j < g.N; j++) { if (g.x[j] > e.xr) T += rx[j]; else if (g.x[j] < e.xl) R += rx[j]; }
    return { T: T * g.dx, R: R * g.dx };
  }
  /** expected final transmission of a packet incident from the left: ∫_{p>0} |ψ̃(p)|² T(p) dp */
  function tInf(g, P, e, Pp) { let s = 0; for (let i = 0; i < Pp.p.length; i++) if (Pp.p[i] > 0) s += Pp.rho[i] * transCoef(P, e, Pp.p[i]); return s * g.dp; }

  /** analytic description at time t (or null): Gaussian {kind:'g', xm, pm, sx, sp} or δ-state {kind:'d', k} */
  function analytic(P, t) {
    const G = gaussInfo(P);
    if (G) {
      const sp = 1 / (2 * G.s0);
      if (P.pot === 'free') return { kind: 'g', xm: G.xm + (G.pm * t) / P.m, pm: G.pm, sx: G.s0 * Math.hypot(1, t / (2 * P.m * G.s0 * G.s0)), sp };
      if (P.pot === 'harmonic') {
        const w = P.omega, c = Math.cos(w * t), s = Math.sin(w * t), mw = P.m * w;
        return { kind: 'g', xm: G.xm * c + (G.pm / mw) * s, pm: G.pm * c - mw * G.xm * s, sx: Math.hypot(G.s0 * c, (sp / mw) * s), sp: Math.hypot(sp * c, mw * G.s0 * s) };
      }
      return t < 1e-12 ? { kind: 'g', xm: G.xm, pm: G.pm, sx: G.s0, sp } : null;
    }
    if (P.state === 'delta' && (P.pot === 'delta' || t < 1e-12)) return { kind: 'd', k: P.m * P.lambda };
    return null;
  }
  function gaussPdf(x, mu, s) { const d = (x - mu) / s; return Math.exp(-d * d / 2) / (Math.sqrt(2 * Math.PI) * s); }
  function aDensX(A, x) { return A.kind === 'g' ? gaussPdf(x, A.xm, A.sx) : A.k * Math.exp(-2 * A.k * Math.abs(x)); }
  function aDensP(A, p) { if (A.kind === 'g') return gaussPdf(p, A.pm, A.sp); const u = p / A.k; return 2 / (Math.PI * A.k) / ((1 + u * u) * (1 + u * u)); }
  function aTailX(A, t) { if (A.kind === 'g') return 0.5 * erfc((t - A.xm) / (SQ2 * A.sx)); return t >= 0 ? 0.5 * Math.exp(-2 * A.k * t) : 1 - 0.5 * Math.exp(2 * A.k * t); }
  function aTailP(A, t) { if (A.kind === 'g') return 0.5 * erfc((t - A.pm) / (SQ2 * A.sp)); const u = t / A.k; return 0.5 - (u / (1 + u * u) + Math.atan(u)) / Math.PI; }
  /** analytic level n (or null) for comparison in the level panel */
  function aLevel(P, n) {
    if (P.pot === 'harmonic') return (n + 0.5) * P.omega;
    if (P.pot === 'box') return ((n + 1) * Math.PI / P.wBox) ** 2 / (2 * P.m);
    if (P.pot === 'delta' && n === 0) return -P.m * P.lambda * P.lambda / 2;
    return null;
  }

  // ---------- UI option lists ----------
  const POTS = [['free', 'free particle'], ['harmonic', 'harmonic ½mω²x²'], ['delta', 'δ-well −λδ(x)'], ['well', 'finite square well'], ['box', 'infinite square well (box)'],
    ['barrier', 'rectangular barrier (tunnelling)'], ['step', 'potential step'], ['double', 'double well a(x²−b²)²'], ['linear', 'linear Fx (gravity)'], ['cos', 'periodic V₀ sin²(πx/d)']];
  const STATES = [['gauss', 'Gaussian packet'], ['delta', 'δ-well bound state'], ['coherent', 'coherent state U|0⟩'], ['hoN', 'HO eigenstate |n⟩ (displaced)'],
    ['cat', 'two Gaussians (cat state)'], ['box', 'box eigenstate n'], ['ground', 'ground state of V (numerical)'], ['eigen', 'eigenstate n of V (numerical)'],
    ['eigen2', 'superposition of two eigenstates of V'], ['exp', 'exponential e^{−κ|x−x₀|}']];
  const opts = (l) => l.map(([value, label]) => ({ value, label }));
  const R_ = (key, label, min, max, step, hint) => ({ key, label, type: 'range', min, max, step, hint });

  // ---------- engine ----------
  QM.registerEngine('wave1d', {
    title: '1D wavefunction lab',
    mount(host, preset) {
      const P = Object.assign({}, DEF, preset || {});
      const { controls, views } = QM.ui.layout(host);
      let g, V, step, psi0, re, im, tCur = 0, S = null, hist = [], yMaxX = 0, yMaxP = 0, tMax = 10, eig = null, eigKey = '', edges = null, Tinf = null;

      // --- controls ---
      const sysPanel = QM.ui.panel('Hamiltonian &amp; initial state');
      function sysSpecs() {
        const s = [{ key: 'pot', label: 'Potential $V(x)$', type: 'select', options: opts(POTS) }, R_('m', 'mass $m$', 0.25, 4, 0.05)], st = P.state, pot = P.pot;
        if (pot === 'harmonic' || ['coherent', 'hoN'].includes(st)) s.push(R_('omega', 'frequency $\\omega$', 0.1, 4, 0.05));
        if (pot === 'delta' || st === 'delta') s.push(R_('lambda', 'δ strength $\\lambda$', 0.1, 4, 0.05, pot === 'delta' ? 'On the grid the δ is one site of depth $\\lambda/\\Delta x$.' : null));
        if (pot === 'well') s.push(R_('V0', 'depth $V_0$', 0.5, 30, 0.5), R_('a', 'width $a$', 0.2, 10, 0.1));
        if (pot === 'box' || st === 'box') s.push(R_('wBox', 'box width $w$', 1, 20, 0.1, pot === 'box' ? 'Hard walls (Dirichlet); evolved exactly as a sine series.' : null));
        if (pot === 'barrier') s.push(R_('Vb', 'barrier height $V_b$', -10, 20, 0.1), R_('bw', 'barrier width $b$', 0.05, 6, 0.05));
        if (pot === 'step') s.push(R_('Vs', 'step height $V_s$ (x > 0)', -10, 20, 0.1));
        if (pot === 'double') s.push(R_('dwA', 'strength $a$', 0.02, 3, 0.01), R_('dwB', 'minima at $\\pm b$', 0.5, 5, 0.05, 'Barrier height $ab^4$ at $x=0$.'));
        if (pot === 'linear') s.push(R_('F', 'force $F$ (slope)', -3, 3, 0.05), { key: 'floor', label: 'hard floor for $x<0$ (quantum bouncer)', type: 'check' });
        if (pot === 'cos') s.push(R_('Vc', 'amplitude $V_0$', 0, 20, 0.1), R_('per', 'period $d$', 0.5, 10, 0.1, '$d$ is rounded so the grid holds an integer number of periods.'));
        s.push({ type: 'heading', label: 'Initial state' }, { key: 'state', label: 'State $\\psi(t=0)$', type: 'select', options: opts(STATES) });
        const has = new Set(s.map((q) => q.key)), add = (...q) => q.forEach((x) => { if (!has.has(x.key)) { has.add(x.key); s.push(x); } });
        const P0 = R_('p0', 'mean momentum $p_0$', -10, 10, 0.01), X0 = R_('x0', 'shift $x_0$', -20, 20, 0.05);
        if (st === 'gauss') add(P0, R_('Delta', 'momentum width $\\Delta$', 0.1, 4, 0.01, '$\\tilde\\psi(p)\\propto e^{-(p-p_0)^2/4\\Delta^2}e^{-ipx_0}$'), X0);
        if (st === 'delta') s.push({ type: 'heading', label: '$\\psi=\\sqrt\\kappa\\,e^{-\\kappa|x|}$, $\\kappa=m\\lambda$' });
        if (st === 'coherent') add(R_('aRe', 'Re $\\alpha$', -3, 3, 0.01), R_('aIm', 'Im $\\alpha$', -3, 3, 0.01, '$U=e^{\\alpha(\\hat a-\\hat a^\\dagger)}=D(-\\alpha)$, so $U|0\\rangle=|\\beta{=}{-\\alpha}\\rangle$, $x_0=1/\\sqrt{m\\omega}$.'));
        if (st === 'hoN') add(R_('nHO', 'level $n$', 0, 10, 1, 'Hermite function of the HO with the current $m,\\omega$, centred at $x_0$.'), X0);
        if (st === 'cat') add(R_('sep', 'separation $s$', 0, 20, 0.1), R_('phi', 'relative phase $\\varphi$', 0, 6.2832, 0.01), R_('Delta', 'momentum width $\\Delta$', 0.1, 4, 0.01), P0, X0);
        if (st === 'box') add(R_('nBox', 'level $n$', 1, 12, 1, '$\\psi\\propto\\sin(n\\pi(x+w/2)/w)$ for $|x|<w/2$.'), R_('wBox', 'box width $w$', 1, 20, 0.1));
        if (st === 'ground') s.push({ type: 'heading', label: 'lowest eigenvector of the current $H$' });
        if (st === 'eigen' || st === 'eigen2') add(R_('n1', 'level $n_1$', 0, NLEV - 1, 1));
        if (st === 'eigen2') add(R_('n2', 'level $n_2$', 0, NLEV - 1, 1), R_('phi', 'relative phase $\\varphi$', 0, 6.2832, 0.01, '$(|n_1\\rangle+e^{i\\varphi}|n_2\\rangle)/\\sqrt2$ oscillates with $\\omega=E_{n_2}-E_{n_1}$.'));
        if (st === 'exp') add(R_('kappa', 'decay $\\kappa$', 0.1, 5, 0.05), X0, P0);
        s.push({ type: 'heading', label: 'Grid' }, R_('L', 'box length $L$ (periodic)', 10, 200, 1, `N = ${P.N} points`));
        return s;
      }
      function buildSys() {
        sysPanel.body.innerHTML = '';
        sysPanel.body.append(QM.ui.params(sysSpecs(), P, (k) => { if (k === 'pot' || k === 'state') buildSys(); reset(k); }).el);
      }
      buildSys();

      const tr = QM.ui.transport({ tMax: 10, onTime: (t) => { advance(t); compute(); record(); redraw(); } });
      const evoPanel = QM.ui.panel('Time evolution');
      const tParams = QM.ui.params([R_('tMax', 'duration $t_{\\max}$', 1, 200, 0.5)], P, () => { tMax = P.tMax; tr.setMax(tMax); redraw(); });
      evoPanel.body.append(tr.el, tParams.el);

      const tailPanel = QM.ui.panel('Tail probabilities');
      const tailOut = el('div');
      tailPanel.body.append(QM.ui.params([
        R_('pStar', 'threshold $p^*$', -10, 10, 0.01), R_('xStar', 'threshold $x^*$', -20, 20, 0.01),
      ], P, () => redraw()).el, tailOut);

      const viewPanel = QM.ui.panel('Display');
      viewPanel.body.append(QM.ui.params([
        R_('xView', 'x range $\\pm$', 1, 100, 0.5), R_('pView', 'p range $\\pm$', 1, 60, 0.5),
        { key: 'showReIm', label: 'Re/Im of $\\psi(x)$', type: 'check' },
        { key: 'showPReIm', label: 'Re/Im of $\\tilde\\psi(p)$', type: 'check' },
        { key: 'analytic', label: 'analytic curves (dashed)', type: 'check' },
      ], P, () => redraw()).el);
      controls.append(sysPanel.el, evoPanel.el, tailPanel.el, viewPanel.el);

      // --- views ---
      const cx = QM.ui.canvas(270), cp = QM.ui.canvas(270), ch = QM.ui.canvas(190), cs = QM.ui.canvas(190), cl = QM.ui.canvas(300), cf = QM.ui.canvas(180);
      const px = QM.ui.panel('Position space $|\\psi(x)|^2$'), pp = QM.ui.panel('Momentum space $|\\tilde\\psi(p)|^2$');
      px.body.append(cx); pp.body.append(cp);
      const st = QM.ui.stats([{ key: 't', label: '$t$' }, { key: 'xm', label: '$\\langle x\\rangle$' }, { key: 'dx', label: '$\\Delta x$' }, { key: 'pm', label: '$\\langle p\\rangle$' },
        { key: 'dp', label: '$\\Delta p$' }, { key: 'unc', label: '$\\Delta x\\,\\Delta p$' }, { key: 'E', label: '$\\langle H\\rangle$' }, { key: 'n', label: 'norm' }]);
      const stPanel = QM.ui.panel('Expectation values');
      stPanel.body.append(st.el);
      const hPanel = QM.ui.panel('$\\langle x\\rangle(t)$ and $\\Delta x(t)$');
      hPanel.body.append(ch, el('div', { class: 'note' }, 'Solid: simulation; dashed: analytic (Gaussian states in free / harmonic potential).'));
      const sPanel = QM.ui.panel('Scattering: transmitted / reflected probability');
      const sNote = el('div', { class: 'note' });
      sPanel.body.append(cs, sNote);
      const lPanel = QM.ui.panel('Energy levels of $H$ and weights $|\\langle E_n|\\psi\\rangle|^2$');
      const lNote = el('div', { class: 'note' });
      lPanel.body.append(cl, lNote);
      const fPanel = QM.ui.panel('Fock-number distribution $|\\langle n|\\psi\\rangle|^2$ (HO with current $m,\\omega$)');
      const fNote = el('div', { class: 'note' });
      fPanel.body.append(cf, fNote);
      views.append(el('div', { class: 'row2' }, px.el, pp.el), stPanel.el, hPanel.el, sPanel.el, lPanel.el, fPanel.el);

      // --- state handling ---
      function reset(key = 'init') {
        if (!g || g.L !== P.L || g.N !== P.N) g = makeGrid(P.N, P.L);
        const kk = JSON.stringify(POT_KEYS.map((k) => P[k]));
        if (kk !== eigKey) { V = potential(g, P); step = propagator(g, V, P.m); eig = eigenstates(g, V, P.m); eigKey = kk; }
        evoPanel.head.firstChild.textContent = `Time evolution (${!V.inner ? 'split-step Fourier' : V.inner.box ? 'exact sine series' : 'Crank–Nicolson'})`;
        psi0 = initState(g, P, eig, V); re = Float64Array.from(psi0.re); im = Float64Array.from(psi0.im);
        tCur = 0; hist = []; yMaxX = 0; yMaxP = 0;
        // automatic duration for periodic motion
        if (['init', 'pot', 'state', 'omega', 'm', 'n1', 'n2', 'wBox'].includes(key)) {
          let auto = null;
          if (P.pot === 'harmonic') auto = (4 * Math.PI) / P.omega;
          else if (P.state === 'eigen2' && Math.round(P.n1) !== Math.round(P.n2)) auto = (4 * Math.PI) / Math.abs(eig.E[Math.round(P.n1)] - eig.E[Math.round(P.n2)]);
          else if (P.pot === 'box') auto = (4 * P.m * ((V.inner.i1 - V.inner.i0 + 2) * g.dx) ** 2) / Math.PI; // revival time
          if (auto) tParams.set('tMax', +Math.min(200, Math.max(1, auto)).toFixed(3), true);
        }
        tMax = P.tMax;
        edges = scatterEdges(g, P);
        const A0 = analyse(g, re, im, V, P.m);
        Tinf = edges && A0.mx.mean + 3 * A0.mx.sd < edges.xl ? tInf(g, P, edges, A0.P) : null;
        tr.pause(); tr.set(0, true); tr.setMax(tMax);
        compute(); record(); redraw();
      }
      function advance(t) {
        if (Math.abs(t - tCur) > t + 1e-12) { re.set(psi0.re); im.set(psi0.im); tCur = 0; }
        const d = t - tCur; if (Math.abs(d) < 1e-14) return;
        const n = Math.ceil(Math.abs(d) / DT); step(re, im, d / n, n); tCur = t;
      }
      function compute() { S = analyse(g, re, im, V, P.m); }
      function record() {
        if (hist.length && hist[hist.length - 1].t >= tCur - 1e-12) hist = hist.filter((h) => h.t < tCur - 1e-12);
        const h = { t: tCur, xm: S.mx.mean, sx: S.mx.sd };
        if (edges) Object.assign(h, sideProbs(g, S.rx, edges));
        hist.push(h);
      }

      // --- drawing ---
      function windowIdx(arr, half) { let i0 = 0, i1 = arr.length - 1; while (i0 < i1 && arr[i0] < -half) i0++; while (i1 > i0 && arr[i1] > half) i1--; return [i0, i1 + 1]; }
      function vRange(a0, a1) { // V range over the window, walls and runaway tails clipped
        let vmin = Infinity, vmax = -Infinity; for (let j = a0; j < a1; j++) { vmin = Math.min(vmin, V[j]); if (V[j] < WALL / 2) vmax = Math.max(vmax, V[j]); }
        if (!(vmax > vmin)) vmax = vmin + 1;
        return [vmin, Math.min(vmax, vmin + 6 * Math.max(S.E - vmin, 1))];
      }
      function redraw() {
        if (!S) return;
        const th = QM.plot.theme(), col = th.series, A = P.analytic ? analytic(P, tCur) : null;
        // x-space
        const [a0, a1] = windowIdx(g.x, P.xView), xs = g.x.subarray(a0, a1), rx = S.rx.subarray(a0, a1);
        let rmax = 0; for (const v of rx) rmax = Math.max(rmax, v);
        const amax = P.showReIm ? Math.sqrt(rmax) : 0;
        yMaxX = Math.max(yMaxX, rmax, amax);
        const sx = [{ x: xs, y: rx, color: col[0], label: '|ψ|²', fillRange: [P.xStar, Infinity] }];
        if (P.showReIm) sx.push({ x: xs, y: re.subarray(a0, a1), color: col[2], width: 1.3, label: 'Re ψ' }, { x: xs, y: im.subarray(a0, a1), color: col[3], width: 1.3, label: 'Im ψ' });
        if (A) sx.push({ x: xs, y: Array.from(xs, (x) => aDensX(A, x)), color: th.ink, width: 1.2, dash: [5, 4], label: 'analytic' });
        if (P.pot !== 'free' && P.pot !== 'delta') {
          const [vmin, vmax] = vRange(a0, a1);
          sx.push({ x: xs, y: Array.from(V.subarray(a0, a1), (v) => yMaxX * (0.03 + (0.85 * (Math.min(v, vmax) - vmin)) / (vmax - vmin))), color: th.muted, width: 1, dash: [2, 3], label: 'V (scaled)' });
        }
        const vl = [{ x: P.xStar, color: col[1] }];
        if (P.pot === 'delta') vl.push({ x: 0, color: th.muted, dash: [2, 3] });
        QM.plot.line(cx, { xMin: -P.xView, xMax: P.xView, yMin: P.showReIm ? -1.1 * Math.max(amax, 1e-9) : 0, yMax: 1.1 * yMaxX, series: sx, vlines: vl, xlabel: 'x' });
        // p-space
        const [b0, b1] = windowIdx(S.P.p, P.pView), ps = S.P.p.subarray(b0, b1), rp = S.P.rho.subarray(b0, b1);
        let pmax = 0; for (const v of rp) pmax = Math.max(pmax, v);
        const pamax = P.showPReIm ? Math.sqrt(pmax) : 0;
        yMaxP = Math.max(yMaxP, pmax, pamax);
        const sp = [{ x: ps, y: rp, color: col[0], label: '|ψ̃|²', fillRange: [P.pStar, Infinity] }];
        if (P.showPReIm) sp.push({ x: ps, y: S.P.re.subarray(b0, b1), color: col[2], width: 1.3, label: 'Re ψ̃' }, { x: ps, y: S.P.im.subarray(b0, b1), color: col[3], width: 1.3, label: 'Im ψ̃' });
        if (A) sp.push({ x: ps, y: Array.from(ps, (p) => aDensP(A, p)), color: th.ink, width: 1.2, dash: [5, 4], label: 'analytic' });
        QM.plot.line(cp, { xMin: -P.pView, xMax: P.pView, yMin: P.showPReIm ? -1.1 * Math.max(pamax, 1e-9) : 0, yMax: 1.1 * yMaxP, series: sp, vlines: [{ x: P.pStar, color: col[1] }], xlabel: 'p' });
        // stats
        const f = (v) => QM.fmt(v, 4);
        st.set('t', f(tCur)); st.set('xm', f(S.mx.mean)); st.set('dx', f(S.mx.sd)); st.set('pm', f(S.mp.mean)); st.set('dp', f(S.mp.sd));
        st.set('unc', f(S.mx.sd * S.mp.sd) + (S.mx.sd * S.mp.sd < 0.5005 ? ' (min ½)' : '')); st.set('E', f(S.E)); st.set('n', QM.fmt(S.mx.norm, 8));
        // tails
        const nP = tail(S.P.p, S.P.rho, g.dp, P.pStar), nX = tail(g.x, S.rx, g.dx, P.xStar);
        const row = (lab, n, a) => `<tr><td>${lab}</td><td class="mono">${QM.fmt(n, 5)}</td><td class="mono">${a === null ? '–' : QM.fmt(a, 5)}</td></tr>`;
        const AA = analytic(P, tCur);
        tailOut.innerHTML = `<table class="meas-table"><tr><th></th><th>numerical</th><th>analytic</th></tr>${row('P(p &gt; p*)', nP, AA ? aTailP(AA, P.pStar) : null)}${row('P(x &gt; x*)', nX, AA ? aTailX(AA, P.xStar) : null)}</table>` +
          (AA && AA.kind === 'g' ? '<div class="note">Gaussian: P(p&gt;p*) = ½ erfc((p*−⟨p⟩)/(√2 Δp)).</div>' : AA ? '<div class="note">δ-state: P(p&gt;p*) = ½ − [u/(1+u²) + arctan u]/π, u = p*/κ.</div>' : '');
        // history
        const ts = hist.map((h) => h.t);
        const ser = [{ x: ts, y: hist.map((h) => h.xm), color: col[0], label: '⟨x⟩' }, { x: ts, y: hist.map((h) => h.sx), color: col[1], label: 'Δx' }];
        if (P.analytic && gaussInfo(P) && (P.pot === 'free' || P.pot === 'harmonic')) {
          const ta = Array.from({ length: 201 }, (_, i) => (i * tMax) / 200), an = ta.map((t) => analytic(P, t));
          ser.push({ x: ta, y: an.map((a) => a.xm), color: col[0], width: 1, dash: [5, 4] }, { x: ta, y: an.map((a) => a.sx), color: col[1], width: 1, dash: [5, 4] });
        }
        QM.plot.line(ch, { xMin: 0, xMax: tMax, series: ser, xlabel: 't' });
        // scattering
        const showS = edges && (P.pot === 'barrier' || P.pot === 'step' || !BOUND.includes(P.state));
        sPanel.el.style.display = showS ? '' : 'none';
        if (showS) {
          const last = hist[hist.length - 1];
          QM.plot.line(cs, { xMin: 0, xMax: tMax, yMin: 0, yMax: 1.05, xlabel: 't', series: [{ x: ts, y: hist.map((h) => h.T), color: col[0], label: 'T = P(x > right edge)' }, { x: ts, y: hist.map((h) => h.R), color: col[1], label: 'R = P(x < left edge)' }],
            hlines: Tinf !== null ? [{ y: Tinf, color: col[0] }, { y: 1 - Tinf, color: col[1] }] : [] });
          sNote.textContent = `T = ${QM.fmt(last.T, 5)},  R = ${QM.fmt(last.R, 5)},  inside = ${QM.fmt(S.mx.norm - last.T - last.R, 5)},  T + R + inside = ${QM.fmt(S.mx.norm, 8)}.` +
            (Tinf !== null ? `  Dashed: plane-wave prediction for t → ∞, T∞ = ∫_{p>0}|ψ̃(p)|²T(p)dp = ${QM.fmt(Tinf, 4)} (R∞ = 1 − T∞).` : '') +
            '  (Finite grid: the transmitted part eventually reflects from the grid edge.)';
        }
        // energy levels
        const showL = !CONTINUUM(P);
        lPanel.el.style.display = showL ? '' : 'none';
        if (showL) drawLevels(th, col, a0, a1, xs);
        // Fock distribution
        const showF = P.pot === 'harmonic' || P.state === 'coherent' || P.state === 'hoN';
        fPanel.el.style.display = showF ? '' : 'none';
        if (showF) {
          const b2 = P.aRe * P.aRe + P.aIm * P.aIm, nMax = P.state === 'coherent' ? Math.min(40, Math.max(12, Math.ceil(b2 + 5 * Math.sqrt(b2) + 6))) : 20;
          const pn = fock(g, re, im, P.m, P.omega, nMax);
          let sum = 0, nm = 0; pn.forEach((v, n) => { sum += v; nm += n * v; });
          const r = QM.plot.bars(cf, { labels: Array.from(pn, (_, n) => String(n)), values: Array.from(pn), yMax: Math.min(1, Math.max(...pn) * 1.2 + 0.02), xlabel: 'n', ylabel: P.state === 'coherent' ? 'bars: numerical · ▬ Poisson, mean |α|²' : null });
          if (P.state === 'coherent') {
            const ctx = cf.getContext('2d'); ctx.strokeStyle = th.ink; ctx.lineWidth = 2; let pois = Math.exp(-b2);
            for (let n = 0; n <= nMax; n++) { if (n > 0) pois *= b2 / n; ctx.beginPath(); ctx.moveTo(r.X(n + 0.1), r.Y(pois)); ctx.lineTo(r.X(n + 0.9), r.Y(pois)); ctx.stroke(); }
          }
          fNote.textContent = `Σ P(n) = ${QM.fmt(sum, 5)} (n ≤ ${nMax}),  ⟨n⟩ = ${QM.fmt(nm, 4)}` + (P.state === 'coherent' ? `,  |α|² = ${QM.fmt(b2, 4)}` : '') +
            (P.pot === 'harmonic' ? '. The distribution is conserved by the HO evolution.' : '.');
        }
      }
      function drawLevels(th, col, a0, a1, xs) {
        const K = eig.E.length, E = eig.E;
        const w = E.map((_, n) => { const u = eig.vec[n]; return dotR(u, re, g.dx) ** 2 + dotR(u, im, g.dx) ** 2; });
        let vmin = Infinity; for (let j = a0; j < a1; j++) vmin = Math.min(vmin, V[j]);
        const lo = Math.min(vmin, E[0]), top = Math.min(E[K - 1], Math.max(E[Math.min(7, K - 1)], S.E + 0.3 * (S.E - E[0]))), yMin = lo - 0.05 * (top - lo), yMax = top + 0.12 * (top - lo);
        const r = QM.plot.line(cl, { xMin: -P.xView, xMax: P.xView, yMin, yMax, xlabel: 'x', ylabel: 'E', legend: false,
          series: [{ x: xs, y: Array.from(V.subarray(a0, a1), (v) => Math.min(v, yMax + (yMax - yMin))), color: th.muted, width: 1.5 }], hlines: [{ y: S.E, color: col[1], dash: [6, 4] }] });
        const ctx = cl.getContext('2d'), wmax = Math.max(...w, 1e-12);
        ctx.save(); ctx.font = '11px system-ui, sans-serif'; ctx.textBaseline = 'bottom'; ctx.textAlign = 'right';
        let lastY = Infinity, shown = 0;
        E.forEach((e, n) => {
          const y = r.Y(e); if (y < r.box.y) return; shown++;
          const lw = 1 + 6 * Math.sqrt(w[n] / wmax) * (w[n] > 1e-4 ? 1 : 0);
          ctx.strokeStyle = w[n] > 1e-4 ? col[0] : th.axis; ctx.lineWidth = lw; ctx.beginPath(); ctx.moveTo(r.box.x, y); ctx.lineTo(r.box.x + r.box.w, y); ctx.stroke();
          if (lastY - y < 13 || y - 13 < r.box.y) return; // skip labels that would overlap
          lastY = y; ctx.fillStyle = th.ink; ctx.fillText(`n=${n}  E=${QM.fmt(e, 4)}  |c|²=${QM.fmt(w[n], 3)}`, r.box.x + r.box.w - 4, y - 2);
        });
        ctx.restore();
        const wsum = w.reduce((a, b) => a + b, 0), an = E.slice(0, 8).map((_, n) => aLevel(P, n));
        lNote.textContent = `Grey: V(x); lines: the ${shown} lowest eigenvalues in range (thickness ∝ weight); dashed: ⟨H⟩ = ${QM.fmt(S.E, 4)}.  Σ_{n<${K}} |c_n|² = ${QM.fmt(wsum, 5)} (conserved in time).` +
          (an[0] !== null ? '  Analytic: ' + an.filter((v) => v !== null).map((v) => QM.fmt(v, 4)).join(', ') + (P.pot === 'box' ? ' (nominal width w; the grid box is slightly wider)' : '') + '.' : '');
      }

      reset();
      const off = QM.ui.onResize(views, redraw);
      return { destroy() { tr.destroy(); off(); }, redraw };
    },
  });

  // ---------- self-tests ----------
  const T = (P0) => Object.assign({}, DEF, P0);
  function prep(P) { const g = makeGrid(P.N, P.L), V = potential(g, P), eig = eigenstates(g, V, P.m), s = initState(g, P, eig, V); return { g, V, s, eig, A: analyse(g, s.re, s.im, V, P.m) }; }
  QM.addTest('wave1d: Gaussian packet is minimum-uncertainty, Δx·Δp = ½, ⟨p⟩ = p0, Δp = Δ', () => {
    const { A } = prep(T({ p0: 2, Delta: 1 })), u = A.mx.sd * A.mp.sd;
    return { pass: Math.abs(u - 0.5) < 1e-4 && Math.abs(A.mp.mean - 2) < 1e-4 && Math.abs(A.mp.sd - 1) < 1e-4 && Math.abs(A.mx.sd - 0.5) < 1e-4, detail: `Δx·Δp = ${u.toFixed(6)}, ⟨p⟩ = ${A.mp.mean.toFixed(5)}, Δp = ${A.mp.sd.toFixed(5)}` };
  });
  QM.addTest('wave1d: P(p > 0) ≈ 0.99 for p0 = 2.33 Δ', () => {
    const { g, A } = prep(T({ p0: 2.33, Delta: 1 })), n = tail(A.P.p, A.P.rho, g.dp, 0), a = 0.5 * (1 + erf(2.33 / SQ2));
    return { pass: Math.abs(n - a) < 5e-4 && Math.abs(a - 0.99) < 1e-3, detail: `numerical ${n.toFixed(5)}, ½[1+erf(p0/√2Δ)] = ${a.toFixed(5)}` };
  });
  QM.addTest('wave1d: δ-well bound state P(p > p0) = 1/4 − 1/(2π) ≈ 0.0908', () => {
    const { g, A } = prep(T({ pot: 'delta', state: 'delta', lambda: 1 })), n = tail(A.P.p, A.P.rho, g.dp, 1), a = 0.25 - 1 / (2 * Math.PI);
    return { pass: Math.abs(n - a) < 1e-3, detail: `numerical ${n.toFixed(5)}, analytic ${a.toFixed(5)}` };
  });
  QM.addTest('wave1d: δ-well bound state is stationary (Crank–Nicolson, t = 5), E = −mλ²/2', () => {
    const { g, V, s, eig, A } = prep(T({ pot: 'delta', state: 'delta', lambda: 1 })), r0 = Float64Array.from(A.rx);
    propagator(g, V, 1)(s.re, s.im, DT, 500); const B = analyse(g, s.re, s.im, V, 1); let dev = 0; for (let j = 0; j < g.N; j++) dev = Math.max(dev, Math.abs(B.rx[j] - r0[j]));
    return { pass: dev < 2e-3 && Math.abs(B.E - A.E) < 1e-3 && Math.abs(eig.E[0] + 0.5) < 2e-3, detail: `max |Δρ| = ${dev.toExponential(2)}, E_0 (FD) = ${eig.E[0].toFixed(5)}` };
  });
  QM.addTest('wave1d: U = e^{α(a−a†)}: ⟨x⟩ = −√2 x0 α, ⟨n⟩ = |α|², after half period ⟨x⟩ = +√2 x0 α', () => {
    const P = T({ pot: 'harmonic', omega: 1, state: 'coherent', aRe: 1, aIm: 0 }), { g, V, s, A } = prep(P);
    const pn = fock(g, s.re, s.im, 1, 1, 30); let nm = 0; pn.forEach((v, n) => (nm += n * v));
    stepper(g, V, 1)(s.re, s.im, Math.PI / 400, 400);
    const B = analyse(g, s.re, s.im, V, 1);
    const ok = Math.abs(A.mx.mean + SQ2) < 1e-6 && Math.abs(nm - 1) < 1e-6 && Math.abs(B.mx.mean - SQ2) < 1e-3 && Math.abs(B.mx.sd - 1 / SQ2) < 1e-3;
    return { pass: ok, detail: `⟨x⟩(0) = ${A.mx.mean.toFixed(5)}, ⟨n⟩ = ${nm.toFixed(5)}, ⟨x⟩(π/ω) = ${B.mx.mean.toFixed(5)}, Δx(π/ω) = ${B.mx.sd.toFixed(5)}` };
  });
  QM.addTest('wave1d: split-step evolution conserves the norm (finite well, 1000 steps)', () => {
    const { g, V, s } = prep(T({ pot: 'well', V0: 5, a: 2, p0: 3 }));
    stepper(g, V, 1)(s.re, s.im, DT, 1000);
    let n = 0; for (let j = 0; j < g.N; j++) n += (s.re[j] * s.re[j] + s.im[j] * s.im[j]) * g.dx;
    return { pass: Math.abs(n - 1) < 1e-9, detail: `‖ψ‖² − 1 = ${(n - 1).toExponential(2)}` };
  });
  QM.addTest('wave1d: numerical ground state of the harmonic potential has E = ω/2', () => {
    const P = T({ pot: 'harmonic', omega: 1.5, state: 'ground' }), { A } = prep(P);
    return { pass: Math.abs(A.E - 0.75) < 1e-3, detail: `E0 = ${A.E.toFixed(6)}` };
  });
  QM.addTest('wave1d: HO eigenstates n = 0…10 have E = n + ½ (Hermite functions ⟨H⟩ and numerical eigenvalues)', () => {
    const P = T({ pot: 'harmonic', omega: 1 }), g = makeGrid(P.N, P.L), V = potential(g, P), eig = eigenstates(g, V, 1);
    let err = 0;
    for (let n = 0; n <= 10; n++) {
      const s = initState(g, T({ pot: 'harmonic', state: 'hoN', nHO: n, x0: 0 }), null), E = analyse(g, s.re, s.im, V, 1).E;
      err = Math.max(err, Math.abs(E - n - 0.5), Math.abs(eig.E[n] - n - 0.5));
    }
    return { pass: err < 1e-3, detail: `max |E_n − (n+½)| = ${err.toExponential(2)}` };
  });
  QM.addTest('wave1d: infinite well (finite-difference solver) levels E_n ∝ n² (n = 1…8), E_1 ≈ π²/(2mw²)', () => {
    const P = T({ pot: 'box', wBox: 6 }), g = makeGrid(P.N, P.L), eig = eigenstates(g, potential(g, P), 1, NLEV, true), E1 = eig.E[0], a = Math.PI ** 2 / 72;
    let err = 0; for (let n = 1; n <= 8; n++) err = Math.max(err, Math.abs(eig.E[n - 1] / E1 / (n * n) - 1));
    return { pass: err < 5e-3 && Math.abs(E1 / a - 1) < 0.03, detail: `max |E_n/(n²E_1) − 1| = ${err.toExponential(2)}, E_1 = ${E1.toFixed(5)} (ideal ${a.toFixed(5)})` };
  });
  QM.addTest('wave1d: box – Gaussian revives at T_rev = 4mw²/π (|⟨ψ(0)|ψ(T_rev)⟩|² = 1), norm kept', () => {
    const P = T({ pot: 'box', wBox: 6, state: 'gauss', p0: 3, x0: -1 }), { g, V, s } = prep(P), r0 = Float64Array.from(s.re), i0 = Float64Array.from(s.im);
    const W = (V.inner.i1 - V.inner.i0 + 2) * g.dx, prop = propagator(g, V, 1);
    prop(s.re, s.im, (4 * W * W) / Math.PI / 3, 1); const A = analyse(g, s.re, s.im, V, 1); prop(s.re, s.im, (8 * W * W) / Math.PI / 3, 1);
    const a = dotR(r0, s.re, g.dx) + dotR(i0, s.im, g.dx), b = dotR(r0, s.im, g.dx) - dotR(i0, s.re, g.dx), ov = a * a + b * b;
    return { pass: Math.abs(ov - 1) < 1e-9 && Math.abs(A.mx.norm - 1) < 1e-9, detail: `overlap ${ov.toFixed(10)}, norm at T_rev/3 = ${A.mx.norm.toFixed(10)}` };
  });
  QM.addTest('wave1d: quantum bouncer (Crank–Nicolson, hard floor) – eigenstate stationary, norm conserved, E_1 ≈ 2.338 (F²/2m)^{1/3}', () => {
    const P = T({ pot: 'linear', floor: true, F: 0.5, state: 'eigen', n1: 0 }), { g, V, s, eig, A } = prep(P), r0 = Float64Array.from(A.rx);
    propagator(g, V, 1)(s.re, s.im, DT, 500); const B = analyse(g, s.re, s.im, V, 1);
    let dev = 0; for (let j = 0; j < g.N; j++) dev = Math.max(dev, Math.abs(B.rx[j] - r0[j]));
    const a = 2.33811 * Math.cbrt(0.125);
    return { pass: dev < 1e-8 && Math.abs(B.mx.norm - 1) < 1e-9 && Math.abs(eig.E[0] / a - 1) < 0.03, detail: `max |Δρ| = ${dev.toExponential(2)}, ‖ψ‖² − 1 = ${(B.mx.norm - 1).toExponential(2)}, E_1 = ${eig.E[0].toFixed(5)} (Airy ${a.toFixed(5)})` };
  });
  QM.addTest('wave1d: barrier scattering conserves norm, T + R = 1, T → ∫|ψ̃(p)|²T(p)dp', () => {
    const P = T({ pot: 'barrier', Vb: 2, bw: 1, state: 'gauss', p0: 2, Delta: 0.25, x0: -15, L: 100 }), { g, V, s, A } = prep(P), e = scatterEdges(g, P), Ti = tInf(g, P, e, A.P);
    propagator(g, V, 1)(s.re, s.im, DT, 1500);
    const B = analyse(g, s.re, s.im, V, 1), { T: Tn, R } = sideProbs(g, B.rx, e);
    return { pass: Math.abs(B.mx.norm - 1) < 1e-9 && Math.abs(Tn + R - 1) < 1e-3 && Math.abs(Tn - Ti) < 0.01, detail: `‖ψ‖² − 1 = ${(B.mx.norm - 1).toExponential(2)}, T = ${Tn.toFixed(5)}, R = ${R.toFixed(5)}, T∞ (plane waves) = ${Ti.toFixed(5)}` };
  });
})();
