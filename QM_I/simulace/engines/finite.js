/* ENGINE 'finite' – finite-dimensional quantum systems (ħ = 1).
   One general lab; the problems are only presets of it.

   Builders (choosable on the page via preset.builderChoices):
     'network'     – N quantum dots: topology chain / ring / complete / star / grid, on-site ε, hopping J,
                     next-nearest γ, tilt F, disorder W; optionally ⊗ spin ½ with λ(2Ŷ−1)⊗σ_x and Zeeman b σ_z
     'spins'       – chain of Ns = 1…6 spins/qubits: level energies ε (+gradient), transverse field,
                     nearest-neighbour Jxx, Jyy, Jzz, antisymmetric ω(σxσy − σyσx), optional switch-off time T
     'orbitalSpin' – orbital l (1…3) ⊗ spin ½, spin–orbit ξ L·S, field b(Lz + 2Sz)

   PRESET SCHEMA
   { builder, builderChoices?: ['network', …], values: {…builder params},
     extraParams?: [paramSpec], labels?: [...], levelChars?: ['0','1'], subLabels?: [...], blochLabels?,
     initial: { mode, ...data }   modes: amp, site, eigen, eig2, packet, twosite, uniform, plane, standing, cat,
                                  random, product, ghz, w, cg, fn   (+ spin direction ts, ps for network⊗spin)
     initialFn?: (values, QM) => Vector, initialLabel?, statePresets?: [{label, vector}],
     observables?: [keys first], track?: [keys], tMax?, showBand? }                                     */
(function () {
  'use strict';
  const { el } = QM.ui;
  const { V, M, C } = QM;
  const builders = {};
  QM.finite = { builders };

  // ------------------------------------------------------------------ helpers
  const fmtE = (v) => QM.fmt(v, 4);
  const piF = (x) => QM.fmt(x / Math.PI, 3) + 'π';
  const LET = 'ABCDEFGHIJKLMNOP';
  function proj(n, k) { const m = M.zeros(n); m.re[k * n + k] = 1; return m; }
  function perm(N, f) { const m = M.zeros(N); for (let i = 0; i < N; i++) m.re[f(i) * N + i] = 1; return m; }
  function rng(seed) { let a = (seed * 2654435761) >>> 0 || 1; return () => { a = (a + 0x6D2B79F5) >>> 0; let t = a; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }
  const qubit = (t, p) => V.from([[Math.cos(t / 2), 0], [Math.sin(t / 2) * Math.cos(p), Math.sin(t / 2) * Math.sin(p)]]);

  // ------------------------------------------------------------------ builder: network of dots
  builders.network = {
    title: 'Quantum-dot network',
    params: () => [
      { key: 'topology', label: 'Topology', type: 'select', options: [{ value: 'chain', label: 'open chain' }, { value: 'ring', label: 'ring' }, { value: 'complete', label: 'complete graph (all-to-all)' }, { value: 'star', label: 'star (centre + N−1 leaves)' }, { value: 'grid', label: '2D square grid' }] },
      { key: 'N', label: 'Number of dots $N$', min: 2, max: 121, step: 1 },
      { key: 'eps', label: 'On-site energy $\\varepsilon$', min: -3, max: 3, step: 0.01 },
      { key: 'J', label: 'Hopping $J$ (nearest neighbours)', min: -2, max: 2, step: 0.01 },
      { key: 'gamma', label: 'Next-nearest hopping $\\gamma$ (ring: diagonals $n\\leftrightarrow n+2$)', min: -2, max: 2, step: 0.01 },
      { key: 'F', label: 'Tilt $F$: $\\varepsilon_n \\mathrel{+}= F\\,n$ (electric field)', min: -1, max: 1, step: 0.005 },
      { key: 'W', label: 'Disorder $W$: random $\\varepsilon_n\\in[-W/2,W/2]$', min: 0, max: 4, step: 0.01 },
      { key: 'seed', label: 'Disorder seed', min: 1, max: 50, step: 1 },
      { type: 'heading', label: 'Spin ½ on the particle' },
      { key: 'spin', label: 'Attach spin ½ (space $\\mathbb C^N\\otimes\\mathbb C^2$)', type: 'check' },
      { key: 'lambda', label: 'Coupling $\\lambda(2\\hat Y-1)\\otimes\\hat\\sigma_x$ ($\\hat Y$ = mirror)', min: -1, max: 1, step: 0.005 },
      { key: 'bz', label: 'Zeeman field $b\\,\\hat\\sigma_z$', min: -1, max: 1, step: 0.01 },
    ],
    defaults: { topology: 'chain', N: 3, eps: 0, J: -1, gamma: 0, F: 0, W: 0, seed: 1, spin: false, lambda: 0.1, bz: 0, labelStyle: 'numbers' },
    build(v) {
      const N = Math.max(2, Math.round(v.N)), topo = v.topology;
      const st = v.labelStyle === 'letters' && N <= 16 ? 'letters' : v.labelStyle === 'centered' ? 'centered' : 'numbers';
      const siteLab = (i) => (st === 'letters' ? LET[i] : st === 'centered' ? String(i - Math.floor(N / 2)) : String(i + 1));
      const siteVal = (i) => (st === 'centered' ? i - Math.floor(N / 2) : i + 1);
      // geometry
      const Lx = Math.ceil(Math.sqrt(N));
      const pos = [...Array(N)].map((_, i) => {
        if (topo === 'chain') return [N === 1 ? 0 : -1 + (2 * i) / (N - 1), 0];
        if (topo === 'star') return i === 0 ? [0, 0] : [Math.cos((2 * Math.PI * (i - 1)) / (N - 1)), Math.sin((2 * Math.PI * (i - 1)) / (N - 1))];
        if (topo === 'grid') { const r = Math.floor(i / Lx), c = i % Lx, s = Math.max(1, Lx - 1); return [-1 + (2 * c) / s, 1 - (2 * r) / s]; }
        return [Math.cos((2 * Math.PI * i) / N), Math.sin((2 * Math.PI * i) / N)];
      });
      // orbital Hamiltonian
      const Hc = M.zeros(N), seen = new Set();
      const hop = (i, j, a) => { if (i === j || !a || i < 0 || j < 0 || i >= N || j >= N) return; const k = Math.min(i, j) + ',' + Math.max(i, j); if (seen.has(k)) return; seen.add(k); Hc.re[i * N + j] += a; Hc.re[j * N + i] += a; };
      const r = rng(Math.round(v.seed || 1));
      for (let i = 0; i < N; i++) Hc.re[i * N + i] = v.eps + (v.F || 0) * (i - (N - 1) / 2) + (v.W ? v.W * (r() - 0.5) : 0);
      if (topo === 'chain') { for (let i = 0; i < N - 1; i++) hop(i, i + 1, v.J); for (let i = 0; i < N - 2; i++) hop(i, i + 2, v.gamma); }
      else if (topo === 'ring') { for (let i = 0; i < N; i++) hop(i, (i + 1) % N, v.J); if (N >= 4) for (let i = 0; i < N; i++) hop(i, (i + 2) % N, v.gamma); }
      else if (topo === 'complete') { for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) hop(i, j, v.J); }
      else if (topo === 'star') { for (let i = 1; i < N; i++) hop(0, i, v.J); for (let i = 1; i < N; i++) hop(i, i === N - 1 ? 1 : i + 1, v.gamma); }
      else if (topo === 'grid') { for (let i = 0; i < N; i++) { const c = i % Lx; if (c < Lx - 1) hop(i, i + 1, v.J); hop(i, i + Lx, v.J); if (c < Lx - 1) hop(i, i + Lx + 1, v.gamma); if (c > 0) hop(i, i + Lx - 1, v.gamma); } }
      const edges = []; for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) { const w = Hc.re[i * N + j]; if (Math.abs(w) > 1e-12) edges.push([i, j, w]); }
      // orbital observables
      const Q = M.diag([...Array(N)].map((_, i) => siteVal(i)));
      const mirror = topo === 'ring' ? perm(N, (i) => (N - i) % N) : perm(N, (i) => N - 1 - i);
      const orb = [];
      if (st === 'letters') orb.push({ key: 'site', label: 'Which dot', op: M.diag([...Array(N)].map((_, i) => i)), fmtValue: (x) => LET[Math.round(x)] });
      else orb.push({ key: 'Q', label: 'Position Q̂ = Σ n|n⟩⟨n|', op: Q });
      if (topo === 'ring' || topo === 'complete' || topo === 'star' || topo === 'grid') {
        orb.push({ key: 'X', label: 'x-coordinate X̂ = Σ x_n|n⟩⟨n|', op: M.diag(pos.map((p) => +p[0].toFixed(12))) });
        orb.push({ key: 'Yc', label: 'y-coordinate Ŷ = Σ y_n|n⟩⟨n|', op: M.diag(pos.map((p) => +p[1].toFixed(12))) });
      }
      if (topo === 'chain' || topo === 'ring') orb.push({ key: 'Y', label: topo === 'ring' ? 'Mirror (n ↔ −n)' : 'Mirror Ŷ (n ↔ N+1−n)', op: mirror });
      if (topo === 'chain') orb.push({ key: 'V', label: 'Velocity V̂ = i[Ĥ, Q̂]', op: M.scale(M.comm(Hc, Q), { re: 0, im: 1 }), fmtValue: fmtE });
      if (topo === 'ring' && N === 4) {
        orb.push({ key: 'PX', label: 'P̂_X = |1⟩⟨3| + |3⟩⟨1| (A↔C)', op: M.from([[0, 0, 1, 0], [0, 0, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0]]) });
        orb.push({ key: 'PY', label: 'P̂_Y = |2⟩⟨4| + |4⟩⟨2| (B↔D)', op: M.from([[0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0], [0, 1, 0, 0]]) });
      }
      const nodes = pos.map((p, i) => ({ x: p[0], y: p[1], label: siteLab(i), idx: [i] }));
      const sys = { n: N, dims: [N], labels: nodes.map((x) => x.label), layout: { nodes, edges }, orbitalDim: N };
      if (!v.spin) {
        sys.H = Hc; sys.observables = [{ key: 'H', label: 'Energy Ĥ', op: Hc, fmtValue: fmtE }].concat(orb);
        if (topo === 'chain' && N > 30 && !v.gamma && !v.F && !v.W) {
          const pairs = [...Array(N)].map((_, k) => ({ E: v.eps + 2 * v.J * Math.cos(((k + 1) * Math.PI) / (N + 1)), k: k + 1 })).sort((a, b) => a.E - b.E);
          sys.eig = { values: pairs.map((p) => p.E), vectors: pairs.map((p) => { const vec = V.zeros(N); for (let n = 0; n < N; n++) vec.re[n] = Math.sqrt(2 / (N + 1)) * Math.sin((p.k * (n + 1) * Math.PI) / (N + 1)); return vec; }) };
        }
      } else {
        const I2 = M.eye(2), IN = M.eye(N), sx = QM.ops.sx(), sy = QM.ops.sy(), sz = QM.ops.sz();
        sys.H = M.add(M.kron(Hc, I2), M.scale(M.kron(M.sub(M.scale(mirror, 2), IN), sx), v.lambda), M.scale(M.kron(IN, sz), v.bz));
        sys.n = 2 * N; sys.dims = [N, 2];
        sys.labels = []; nodes.forEach((nd) => { sys.labels.push(nd.label + '↑', nd.label + '↓'); });
        nodes.forEach((nd, i) => (nd.idx = [2 * i, 2 * i + 1]));
        sys.observables = [{ key: 'H', label: 'Energy Ĥ', op: sys.H, fmtValue: fmtE },
          { key: 'sx', label: 'Spin ŝ_x', op: M.kron(IN, M.scale(sx, 0.5)) }, { key: 'sy', label: 'Spin ŝ_y', op: M.kron(IN, M.scale(sy, 0.5)) }, { key: 'sz', label: 'Spin ŝ_z', op: M.kron(IN, M.scale(sz, 0.5)) }]
          .concat(orb.map((o) => Object.assign({}, o, { op: M.kron(o.op, I2), label: o.label + ' ⊗ 1' })));
        sys.qubits = [{ sub: 1, label: 'spin (reduced)' }]; sys.entSub = 1; sys.spin = true;
      }
      if (topo === 'chain' || topo === 'ring') sys.band = { eps: v.eps, J: v.J, gamma: v.gamma };
      return sys;
    },
  };

  // ------------------------------------------------------------------ builder: spin chain / qubits
  builders.spins = {
    title: 'Spins / qubits',
    params: () => [
      { key: 'Ns', label: 'Number of spins (qubits) $N_s$', min: 1, max: 6, step: 1 },
      { type: 'heading', label: 'Single-spin terms' },
      { key: 'e', label: 'Energy of upper level $\\varepsilon\\,|1\\rangle\\langle1|$', min: -3, max: 3, step: 0.01 },
      { key: 'de', label: 'Gradient: $\\varepsilon_k=\\varepsilon+k\\,\\delta\\varepsilon$', min: -2, max: 2, step: 0.01 },
      { key: 'hx', label: 'Transverse field $h_x\\,\\sigma_x$', min: -2, max: 2, step: 0.01 },
      { type: 'heading', label: 'Nearest-neighbour couplings' },
      { key: 'Jxx', label: '$J_{xx}\\,\\sigma_x\\sigma_x$', min: -2, max: 2, step: 0.01 },
      { key: 'Jyy', label: '$J_{yy}\\,\\sigma_y\\sigma_y$', min: -2, max: 2, step: 0.01 },
      { key: 'Jzz', label: '$J_{zz}\\,\\sigma_z\\sigma_z$', min: -2, max: 2, step: 0.01 },
      { key: 'D', label: '$\\omega\\,(\\sigma_x\\sigma_y-\\sigma_y\\sigma_x)$', min: -2, max: 2, step: 0.01 },
      { key: 'closed', label: 'Close into a ring (last–first bond)', type: 'check' },
      { key: 'switchOn', label: 'Couplings only for $0\\le t\\le T$', type: 'check' },
      { key: 'T', label: 'Switch-off time $T$', min: 0, max: 6, step: 0.001 },
      { type: 'heading', label: 'Measurement axis $\\vec n$' },
      { key: 'thN', label: '$\\theta_n$', min: 0, max: Math.PI, step: 0.01, fmt: piF },
      { key: 'phN', label: '$\\varphi_n$', min: 0, max: 2 * Math.PI, step: 0.01, fmt: piF },
    ],
    defaults: { Ns: 2, e: 0, de: 0, hx: 0, Jxx: 0, Jyy: 0, Jzz: 0, D: 0, closed: false, switchOn: false, T: Math.PI / 8, thN: 0, phN: 0 },
    build(v, preset) {
      const Ns = Math.max(1, Math.round(v.Ns)), dims = Array(Ns).fill(2), n = 1 << Ns;
      const I = M.eye(2), sx = QM.ops.sx(), sy = QM.ops.sy(), sz = QM.ops.sz(), P1 = proj(2, 1);
      const loc = (op, k) => QM.ops.embed(op, dims, k);
      let H0 = M.zeros(n), Vc = M.zeros(n);
      for (let k = 0; k < Ns; k++) { H0 = M.add(H0, M.scale(loc(P1, k), v.e + k * v.de), M.scale(loc(sx, k), v.hx)); }
      const bonds = []; for (let k = 0; k < Ns - 1; k++) bonds.push([k, k + 1]); if (v.closed && Ns > 2) bonds.push([Ns - 1, 0]);
      const two = (A, a, B, b) => M.mul(loc(A, a), loc(B, b));
      bonds.forEach(([a, b]) => { Vc = M.add(Vc, M.scale(two(sx, a, sx, b), v.Jxx), M.scale(two(sy, a, sy, b), v.Jyy), M.scale(two(sz, a, sz, b), v.Jzz), M.scale(M.sub(two(sx, a, sy, b), two(sy, a, sx, b)), v.D)); });
      const H = M.add(H0, Vc);
      const sn = QM.ops.sn(v.thN, v.phN);
      const tot = (op) => { let s = M.zeros(n); for (let k = 0; k < Ns; k++) s = M.add(s, loc(M.scale(op, 0.5), k)); return s; };
      const Sx = tot(sx), Sy = tot(sy), Sz = tot(sz);
      const ch = (preset && preset.levelChars) || ['0', '1'];
      const labels = preset && preset.labels && preset.labels.length === n ? preset.labels : [...Array(n)].map((_, i) => [...Array(Ns)].map((_, k) => ch[(i >> (Ns - 1 - k)) & 1]).join(''));
      const sub = (k) => (preset && preset.subLabels && preset.subLabels[k]) || 'spin ' + (k + 1);
      const obs = [{ key: 'H', label: 'Energy Ĥ', op: H, fmtValue: fmtE }];
      if (v.switchOn) obs.push({ key: 'H0', label: 'Ĥ₀ (without couplings)', op: H0, fmtValue: fmtE });
      for (let k = 0; k < Ns; k++) {
        obs.push({ key: 'E' + (k + 1), label: `Energy of ${sub(k)}: ε_k|1⟩⟨1|`, op: M.scale(loc(P1, k), v.e + k * v.de), fmtValue: fmtE });
        obs.push({ key: 'Z' + (k + 1), label: `Level of ${sub(k)} (|0⟩→0, |1⟩→1)`, op: loc(P1, k) });
        obs.push({ key: 'Sn' + (k + 1), label: `Ŝ_n of ${sub(k)} (½σ_n)`, op: loc(M.scale(sn, 0.5), k) });
      }
      obs.push({ key: 'Sn', label: 'Total Ŝ_n = Ŝ·n', op: tot(sn) }, { key: 'Sz', label: 'Total Ŝ_z', op: Sz },
        { key: 'S2', label: 'Total Ŝ² (= S(S+1))', op: M.add(M.mul(Sx, Sx), M.mul(Sy, Sy), M.mul(Sz, Sz)) },
        { key: 'A', label: 'Â = Ŝ_x² − Ŝ_y² (×ω/ħ, ω=1)', op: M.sub(M.mul(Sx, Sx), M.mul(Sy, Sy)), fmtValue: (x) => QM.fmt(x, 4) + ' ω' });
      const sys = { n, dims, labels, H, observables: obs, qubits: dims.map((_, k) => ({ sub: k, label: sub(k) })), entSub: 0, Ns,
        popOps: Ns > 2 ? dims.map((_, k) => ({ label: `P(${sub(k)} in |${ch[1]}⟩)`, op: loc(P1, k) })) : null,
        axis: [Math.sin(v.thN) * Math.cos(v.phN), Math.sin(v.thN) * Math.sin(v.phN), Math.cos(v.thN)] };
      if (Ns <= 3) sys.layout = { nodes: [...Array(Ns)].map((_, k) => ({ x: Ns === 1 ? 0 : -1 + (2 * k) / (Ns - 1), y: 0, label: sub(k), idx: [] })), edges: bonds.map(([a, b]) => [a, b, 1]) };
      if (v.switchOn) { sys.Hoff = H0; sys.tSwitch = v.T; }
      return sys;
    },
  };

  // ------------------------------------------------------------------ builder: orbital ⊗ spin
  /** CG state |j mj> in (orbital l) ⊗ (spin ½), basis m_l = l..−l, spin ↑,↓ */
  function cgState(l, j, mj) {
    const n = 2 * (2 * l + 1), v = V.zeros(n), idx = (ml, up) => (l - ml) * 2 + (up ? 0 : 1);
    const a = Math.sqrt((l + mj + 0.5) / (2 * l + 1)), b = Math.sqrt((l - mj + 0.5) / (2 * l + 1));
    const mUp = mj - 0.5, mDn = mj + 0.5;
    if (Math.abs(j - (l + 0.5)) < 1e-9) { if (Math.abs(mUp) <= l) v.re[idx(mUp, true)] = a; if (Math.abs(mDn) <= l) v.re[idx(mDn, false)] = b; }
    else { if (Math.abs(mUp) <= l) v.re[idx(mUp, true)] = -b; if (Math.abs(mDn) <= l) v.re[idx(mDn, false)] = a; }
    return V.normalize(v);
  }
  QM.finite.cgState = cgState;

  builders.orbitalSpin = {
    title: 'Orbital l ⊗ spin ½',
    params: () => [
      { key: 'l', label: 'Orbital quantum number $l$', min: 1, max: 3, step: 1 },
      { key: 'xi', label: 'Spin–orbit $\\xi\\,\\hat{\\vec L}\\cdot\\hat{\\vec S}$', min: -1, max: 1, step: 0.01 },
      { key: 'bz', label: 'Magnetic field $b\\,(\\hat L_z+2\\hat S_z)$', min: -1, max: 1, step: 0.01 },
      { type: 'heading', label: 'Spin measurement axis $\\vec n$' },
      { key: 'thN', label: '$\\theta_n$', min: 0, max: Math.PI, step: 0.01, fmt: piF },
      { key: 'phN', label: '$\\varphi_n$', min: 0, max: 2 * Math.PI, step: 0.01, fmt: piF },
    ],
    defaults: { l: 1, xi: 0.2, bz: 0.1, thN: Math.PI / 4, phN: 0 },
    build(v) {
      const l = Math.round(v.l), L = QM.ops.spin(l), S = QM.ops.spin(0.5), d = 2 * l + 1, Il = M.eye(d), Is = M.eye(2);
      const Jx = M.add(M.kron(L.Jx, Is), M.kron(Il, S.Jx)), Jy = M.add(M.kron(L.Jy, Is), M.kron(Il, S.Jy)), Jz = M.add(M.kron(L.Jz, Is), M.kron(Il, S.Jz));
      const J2 = M.add(M.mul(Jx, Jx), M.mul(Jy, Jy), M.mul(Jz, Jz));
      const LS = M.add(M.kron(L.Jx, S.Jx), M.kron(L.Jy, S.Jy), M.kron(L.Jz, S.Jz));
      const H = M.add(M.scale(LS, v.xi), M.scale(M.add(M.kron(L.Jz, Is), M.scale(M.kron(Il, S.Jz), 2)), v.bz));
      const labels = []; for (let k = 0; k < d; k++) for (const s of ['↑', '↓']) labels.push(`${QM.fmt(l - k, 0)}${s}`);
      return {
        n: 2 * d, dims: [d, 2], labels, H, l,
        observables: [
          { key: 'Sn', label: 'Spin Ŝ_n (1 ⊗ ½σ_n)', op: M.kron(Il, M.scale(QM.ops.sn(v.thN, v.phN), 0.5)) },
          { key: 'J2', label: 'Ĵ² (= j(j+1))', op: J2 }, { key: 'Jz', label: 'Ĵ_z', op: Jz },
          { key: 'Lz', label: 'L̂_z', op: M.kron(L.Jz, Is) }, { key: 'Sz', label: 'Ŝ_z', op: M.kron(Il, S.Jz) },
          { key: 'LS', label: 'L̂·Ŝ', op: LS }, { key: 'H', label: 'Energy Ĥ', op: H, fmtValue: fmtE },
        ],
        qubits: [{ sub: 1, label: 'spin (reduced)' }], entSub: 1,
        axis: [Math.sin(v.thN) * Math.cos(v.phN), Math.sin(v.thN) * Math.sin(v.phN), Math.cos(v.thN)],
      };
    },
  };

  // ------------------------------------------------------------------ initial states
  const MODE_NAMES = {
    amp: 'amplitudes (Re/Im)', site: 'basis state', eigen: 'energy eigenstate', eig2: 'superposition of two eigenstates',
    packet: 'Gaussian wave packet', twosite: 'a|n₀⟩ + b|n₀+1⟩', uniform: 'uniform superposition', plane: 'plane wave e^{ikn}',
    standing: 'standing wave sin(kπn/(N+1))', cat: 'two packets ("cat")', random: 'random state',
    product: 'product state (Bloch vectors)', ghz: 'GHZ (|0…0⟩ + e^{iφ}|1…1⟩)/√2', w: 'W state', cg: 'coupled state |j m_j⟩',
  };
  function modesFor(builderName, sys, preset) {
    let m;
    if (builderName === 'network') m = ['site', 'packet', 'twosite', 'uniform', 'plane', 'standing', 'cat', 'eigen', 'eig2', 'random', 'amp'];
    else if (builderName === 'spins') m = ['product', 'site', 'ghz', 'w', 'eigen', 'eig2', 'random', 'amp'];
    else m = ['cg', 'site', 'eigen', 'eig2', 'random', 'amp'];
    if (sys.n > 16) m = m.filter((x) => x !== 'amp');
    if (sys.Ns === 1) m = m.filter((x) => x !== 'ghz' && x !== 'w');
    return (preset.initialFn ? ['fn'] : []).concat(m);
  }
  function initSpecs(mode, sys) {
    const No = sys.orbitalDim || sys.n;
    const k0 = { key: 'k0', label: 'Quasi-momentum $k_0$', min: -Math.PI, max: Math.PI, step: 0.01, fmt: piF };
    const sp = {
      packet: [{ key: 'center', label: 'Centre $n_0$', min: 0, max: No - 1, step: 1 }, { key: 'width', label: 'Width $\\sigma$ (sites)', min: 0.3, max: Math.max(2, No / 3), step: 0.1 }, k0],
      twosite: [{ key: 'center', label: 'Site $n_0$', min: 0, max: No - 2, step: 1 }, { key: 'theta', label: 'Mixing: $|a|=\\cos\\frac\\theta2,\\ |b|=\\sin\\frac\\theta2$', min: 0, max: Math.PI, step: 0.01, fmt: piF }, { key: 'phi', label: 'Relative phase $\\arg b$', min: 0, max: 2 * Math.PI, step: 0.01, fmt: piF }],
      plane: [k0], standing: [{ key: 'kk', label: 'Mode number $k$', min: 1, max: No, step: 1 }],
      cat: [{ key: 'center', label: 'Packet 1 centre', min: 0, max: No - 1, step: 1 }, { key: 'center2', label: 'Packet 2 centre', min: 0, max: No - 1, step: 1 }, { key: 'width', label: 'Width $\\sigma$', min: 0.3, max: Math.max(2, No / 4), step: 0.1 }, k0, { key: 'phi', label: 'Relative phase', min: 0, max: 2 * Math.PI, step: 0.01, fmt: piF }],
      random: [{ key: 'seed', label: 'Random seed', min: 1, max: 100, step: 1 }],
      eig2: [{ key: 'k1', label: 'Eigenstate $k_1$ (0 = ground)', min: 0, max: sys.n - 1, step: 1 }, { key: 'k2', label: 'Eigenstate $k_2$', min: 0, max: sys.n - 1, step: 1 }, { key: 'phi', label: 'Relative phase', min: 0, max: 2 * Math.PI, step: 0.01, fmt: piF }],
      ghz: [{ key: 'phi', label: 'Phase $\\varphi$', min: 0, max: 2 * Math.PI, step: 0.01, fmt: piF }],
      product: [...Array(sys.Ns || 0)].flatMap((_, k) => [{ key: 't' + (k + 1), label: `${sys.qubits[k].label}: $\\theta$`, min: 0, max: Math.PI, step: 0.01, fmt: piF }, { key: 'p' + (k + 1), label: `${sys.qubits[k].label}: $\\varphi$`, min: 0, max: 2 * Math.PI, step: 0.01, fmt: piF }]),
      cg: sys.l ? [{ key: 'j', label: '$j$', type: 'select', options: [{ value: sys.l + 0.5, label: `l + ½ = ${sys.l + 0.5}` }, { value: sys.l - 0.5, label: `l − ½ = ${sys.l - 0.5}` }] }, { key: 'mj', label: '$m_j$', min: -(sys.l + 0.5), max: sys.l + 0.5, step: 1 }] : [],
    }[mode] || [];
    const orbital = ['packet', 'twosite', 'plane', 'standing', 'cat', 'uniform'].includes(mode) || (mode === 'random' && sys.spin);
    if (sys.spin && orbital) sp.push({ type: 'heading', label: 'Spin direction' }, { key: 'ts', label: 'Spin $\\theta$', min: 0, max: Math.PI, step: 0.01, fmt: piF }, { key: 'ps', label: 'Spin $\\varphi$', min: 0, max: 2 * Math.PI, step: 0.01, fmt: piF });
    return sp;
  }
  function makeInitial(initial, sys, eig, values, preset) {
    const n = sys.n, No = sys.orbitalDim || n, m = initial.mode, g = (k, d) => (initial[k] === undefined ? d : +initial[k]);
    const withSpin = (orb) => (sys.spin ? V.kron(V.normalize(orb), qubit(g('ts', 0), g('ps', 0))) : V.normalize(orb));
    const pk = (c, w, k) => { const v = V.zeros(No); for (let i = 0; i < No; i++) { const a = Math.exp(-((i - c) ** 2) / (4 * w * w)); v.re[i] = a * Math.cos(k * i); v.im[i] = a * Math.sin(k * i); } return v; };
    switch (m) {
      case 'amp': return initial.vector && initial.vector.length === n ? V.normalize(V.from(initial.vector)) : V.basis(n, 0);
      case 'site': return V.basis(n, Math.min(n - 1, g('k', 0)));
      case 'eigen': return V.clone(eig.vectors[Math.min(n - 1, g('k', 0))]);
      case 'eig2': { const v = V.clone(eig.vectors[Math.min(n - 1, g('k1', 0))]); V.addTo(v, eig.vectors[Math.min(n - 1, g('k2', 1))], C.expi(g('phi', 0))); return V.normalize(v); }
      case 'packet': return withSpin(pk(g('center', Math.floor(No / 2)), g('width', 1.5), g('k0', 0)));
      case 'cat': { const v = pk(g('center', Math.floor(No / 4)), g('width', 1.5), g('k0', 0)); V.addTo(v, pk(g('center2', Math.floor((3 * No) / 4)), g('width', 1.5), -g('k0', 0)), C.expi(g('phi', 0))); return withSpin(v); }
      case 'twosite': { const c = Math.min(No - 2, g('center', Math.floor(No / 2))), v = V.zeros(No), th = g('theta', 0), ph = g('phi', 0); v.re[c] = Math.cos(th / 2); v.re[c + 1] = Math.sin(th / 2) * Math.cos(ph); v.im[c + 1] = Math.sin(th / 2) * Math.sin(ph); return withSpin(v); }
      case 'uniform': { const v = V.zeros(No); v.re.fill(1); return withSpin(v); }
      case 'plane': { const v = V.zeros(No); for (let i = 0; i < No; i++) { v.re[i] = Math.cos(g('k0', 0) * i); v.im[i] = Math.sin(g('k0', 0) * i); } return withSpin(v); }
      case 'standing': { const v = V.zeros(No); for (let i = 0; i < No; i++) v.re[i] = Math.sin((g('kk', 1) * Math.PI * (i + 1)) / (No + 1)); return withSpin(v); }
      case 'random': { const r = rng(g('seed', 1)), v = V.zeros(n); for (let i = 0; i < n; i++) { const a = Math.sqrt(-2 * Math.log(r() + 1e-12)), b = 2 * Math.PI * r(); v.re[i] = a * Math.cos(b); v.im[i] = a * Math.sin(b); } return V.normalize(v); }
      case 'product': { let v = null; for (let k = 0; k < sys.Ns; k++) { const q = qubit(g('t' + (k + 1), 0), g('p' + (k + 1), 0)); v = v ? V.kron(v, q) : q; } return v; }
      case 'ghz': { const v = V.zeros(n); v.re[0] = Math.SQRT1_2; v.re[n - 1] = Math.SQRT1_2 * Math.cos(g('phi', 0)); v.im[n - 1] = Math.SQRT1_2 * Math.sin(g('phi', 0)); return v; }
      case 'w': { const v = V.zeros(n); for (let k = 0; k < sys.Ns; k++) v.re[1 << k] = 1; return V.normalize(v); }
      case 'cg': { const j = g('j', sys.l + 0.5), mj = Math.max(-j, Math.min(j, Math.round(g('mj', 0.5) - 0.5) + 0.5)); return cgState(sys.l, j, mj); }
      case 'fn': return V.normalize(preset.initialFn(values, QM));
      default: return V.basis(n, 0);
    }
  }

  // ------------------------------------------------------------------ engine
  function mount(host, preset0) {
    let inner = null;
    const switchTo = (b) => {
      if (inner) inner.destroy();
      host.innerHTML = '';
      const p = b === preset0.builder ? preset0 : { builder: b, builderChoices: preset0.builderChoices, values: {}, blochLabels: preset0.blochLabels };
      inner = mountInner(host, p, switchTo);
      QM.ui.tex(host);
    };
    switchTo(preset0.builder);
    return { destroy: () => inner && inner.destroy(), redraw: () => inner && inner.redraw() };
  }

  function mountInner(host, preset, switchTo) {
    const builderName = preset.builder, builder = builders[builderName];
    const values = Object.assign({}, builder.defaults, preset.values || {});
    const hide = new Set(preset.hide || []);
    const specs = builder.params().filter((s) => !hide.has(s.key))
      .filter((s, i, arr) => s.type !== 'heading' || (arr[i + 1] && arr[i + 1].type !== 'heading'))
      .concat(preset.extraParams || []);
    let sys, eigOn, eigOff, psiRef, tRef = 0, psiAt, traj = null;
    const initial = Object.assign({}, preset.initial || {});
    const cleanups = [];

    const { controls, views } = QM.ui.layout(host);
    // --- controls
    const sysPanel = QM.ui.panel('Hamiltonian');
    if (preset.builderChoices && preset.builderChoices.length > 1) {
      const bs = el('select', { class: 'param-select' }, preset.builderChoices.map((b) => el('option', { value: b }, builders[b].title)));
      bs.value = builderName;
      bs.addEventListener('change', () => switchTo(bs.value));
      sysPanel.body.append(el('div', { class: 'param', style: { marginBottom: '10px' } }, el('label', { class: 'param-label' }, 'System'), bs));
    }
    const paramUI = QM.ui.params(specs, values, onParam);
    sysPanel.body.append(paramUI.el);
    const initPanel = QM.ui.panel('Initial state $|\\psi(0)\\rangle$');
    const initBody = el('div'); initPanel.body.append(initBody);
    const timePanel = QM.ui.panel('Time evolution');
    const transport = QM.ui.transport({ tMax: preset.tMax || 10, onTime: () => drawAll(false), label: 'Time $t$ ($\\hbar=1$)' });
    const tMaxUI = QM.ui.params([{ key: 'tMax', label: 'Time window $t_{\\max}$', min: 1, max: 100, step: 0.5 }], { tMax: preset.tMax || 10 }, () => { transport.setMax(tMaxUI.values.tMax); recomputeTraj(); drawAll(true); });
    timePanel.body.append(transport.el, tMaxUI.el);
    const measPanel = QM.ui.panel('Measurement');
    const meas = QM.measurePanel({
      observables: () => { const keys = preset.observables || []; return keys.map((k) => sys.observables.find((o) => o.key === k)).filter(Boolean).concat(sys.observables.filter((o) => !keys.includes(o.key))); },
      getState: () => psiAt(transport.t),
      onCollapse: (psi) => { psiRef = psi; tRef = transport.t; recomputeTraj(); drawAll(true); },
      get basisLabels() { return sys.labels; },
    });
    measPanel.body.append(meas.el, el('p', { class: 'note' }, 'After a measurement the state is replaced by the collapsed state at the current time; evolution continues from it. Changing the initial state undoes it.'));
    controls.append(sysPanel.el, initPanel.el, timePanel.el, measPanel.el);

    // --- views
    const stateCanvas = QM.ui.canvas(250), barCanvas = QM.ui.canvas(190);
    const ketBox = el('div', { class: 'mono note', style: { wordBreak: 'break-word' } });
    const statsUI = QM.ui.stats([{ key: 't', label: 'time $t$' }, { key: 'E', label: '$\\langle H\\rangle$' }, { key: 'dE', label: '$\\Delta H$' }, { key: 'norm', label: '$\\|\\psi\\|^2$' }, { key: 'ent', label: 'entanglement' }]);
    const statePanel = QM.ui.panel('State at time $t$');
    statePanel.body.append(stateCanvas, barCanvas, ketBox, statsUI.el);
    const popCanvas = QM.ui.canvas(240), trackCanvas = QM.ui.canvas(190);
    const trackSel = el('select', { class: 'param-select', style: { width: 'auto' } });
    const popPanel = QM.ui.panel('Probabilities and expectation values over time', { right: trackSel });
    popPanel.body.append(popCanvas, trackCanvas);
    const levelCanvas = QM.ui.canvas(230), eigTable = el('div', { class: 'matrix-wrap' });
    const levelPanel = QM.ui.panel('Spectrum of $\\hat H$ and weights $|\\langle E|\\psi\\rangle|^2$');
    levelPanel.body.append(levelCanvas, eigTable);
    const blochWrap = el('div', { style: { display: 'grid', gap: '10px' } });
    const blochPanel = QM.ui.panel('Subsystems – reduced density matrices');
    blochPanel.body.append(blochWrap);
    const bandCanvas = QM.ui.canvas(200);
    const bandPanel = QM.ui.panel('Band structure $E(k)$ and group velocity $v(k)$ (infinite, clean chain)');
    bandPanel.body.append(bandCanvas);
    views.append(statePanel.el, popPanel.el, blochPanel.el, levelPanel.el, bandPanel.el);
    let blochs = [], trackKey = (preset.track || [])[0] || '';
    trackSel.addEventListener('change', () => { trackKey = trackSel.value; recomputeTraj(); drawAll(false); });

    // ------------------------------------------------------------------ logic
    function rebuild() {
      sys = builder.build(values, preset);
      eigOn = sys.eig || QM.eigh(sys.H);
      eigOff = sys.Hoff ? QM.eigh(sys.Hoff) : null;
      blochPanel.el.style.display = sys.qubits ? '' : 'none';
      bandPanel.el.style.display = sys.band ? '' : 'none';
      if (sys.qubits) {
        const cols = Math.min(3, sys.qubits.length);
        blochWrap.style.gridTemplateColumns = `repeat(${cols}, minmax(0, 1fr))`;
        if (blochs.length !== sys.qubits.length || blochs.some((b, i) => b.q.label !== sys.qubits[i].label)) {
          blochWrap.innerHTML = '';
          blochs = sys.qubits.map((q) => { const c = QM.ui.canvas(sys.qubits.length > 2 ? 190 : 230); const info = el('div', { class: 'note mono' }); blochWrap.append(el('div', null, c, info)); return { q, b: QM.bloch(c, { labels: preset.blochLabels }), info }; });
        }
      }
      // tracked observable selector
      const opts = [{ key: '', label: '⟨O⟩(t): none' }].concat(sys.observables.map((o) => ({ key: o.key, label: '⟨O⟩(t): ' + o.label })));
      trackSel.innerHTML = ''; opts.forEach((o) => trackSel.append(el('option', { value: o.key }, o.label)));
      if (!sys.observables.some((o) => o.key === trackKey)) trackKey = '';
      trackSel.value = trackKey;
    }
    function onParam(key) {
      const structural = ['N', 'topology', 'l', 'Ns', 'spin'].includes(key);
      rebuild();
      if (structural || initial.mode === 'fn' || initial.mode === 'eigen' || initial.mode === 'eig2' || (preset.resetOn || []).includes(key)) { buildInitUI(); setInitial(); }
      else recomputeTraj();
      meas.refresh();
      drawAll(true);
    }

    function buildInitUI() {
      initBody.innerHTML = '';
      const modes = modesFor(builderName, sys, preset);
      if (!modes.includes(initial.mode)) initial.mode = modes[0];
      const names = Object.assign({}, MODE_NAMES, { fn: preset.initialLabel || 'problem state' });
      const sel = el('select', { class: 'param-select' }, modes.map((m) => el('option', { value: m }, names[m])));
      sel.value = initial.mode;
      sel.addEventListener('change', () => { initial.mode = sel.value; buildInitUI(); setInitial(); meas.update(); drawAll(true); });
      initBody.append(el('div', { class: 'param' }, el('label', { class: 'param-label' }, 'Type'), sel));
      const sub = el('div', { style: { marginTop: '8px' } });
      initBody.append(sub);
      const change = () => { setInitial(); meas.update(); drawAll(true); };
      const m = initial.mode;
      if (m === 'amp') {
        const vec = initial.vector && initial.vector.length === sys.n ? V.from(initial.vector) : V.basis(sys.n, 0);
        sub.append(QM.ui.amplitudes({ labels: sys.labels.map((l) => '|' + l + '⟩'), vector: vec, presets: preset.statePresets, onChange: (v) => { initial.vector = [...Array(v.n)].map((_, i) => [v.re[i], v.im[i]]); change(); } }).el);
      } else if (m === 'site' || m === 'eigen') {
        const s2 = el('select', { class: 'param-select' }, [...Array(sys.n)].map((_, i) => el('option', { value: i }, m === 'site' ? '|' + sys.labels[i] + '⟩' : `#${i}: E = ${fmtE(eigOn.values[i])}`)));
        s2.value = Math.min(initial.k || 0, sys.n - 1);
        s2.addEventListener('change', () => { initial.k = +s2.value; change(); });
        sub.append(s2);
      } else if (m === 'fn') {
        sub.append(el('p', { class: 'note' }, 'Built from the problem parameters in the Hamiltonian panel.'));
      }
      const sp = initSpecs(m, sys);
      if (sp.length) sub.append(QM.ui.params(sp, initial, change).el);
      if (preset.statePresets && m !== 'amp' && sys.n <= 16) {
        const ps = el('select', { class: 'param-select', style: { marginTop: '8px' } }, el('option', { value: '' }, 'Problem states…'), preset.statePresets.map((p, i) => el('option', { value: i }, p.label)));
        ps.addEventListener('change', () => { if (ps.value === '') return; initial.mode = 'amp'; initial.vector = preset.statePresets[+ps.value].vector.map((x) => (typeof x === 'number' ? [x, 0] : x)); buildInitUI(); change(); });
        sub.append(ps);
      }
      QM.ui.tex(initBody);
    }
    function setInitial() { psiRef = makeInitial(initial, sys, eigOn, values, preset); tRef = 0; transport.set(0, true); recomputeTraj(); }

    // evolution ---------------------------------------------------------
    function recomputeTraj() {
      if (!eigOff) { const f = QM.evolver(eigOn, psiRef); psiAt = (t) => f(t - tRef); }
      else {
        const T = sys.tSwitch;
        if (tRef <= T) { const fOn = QM.evolver(eigOn, psiRef), psiT = fOn(T - tRef), fOff = QM.evolver(eigOff, psiT); psiAt = (t) => (t <= T ? fOn(t - tRef) : fOff(t - T)); }
        else { const fOff = QM.evolver(eigOff, psiRef), psiT = fOff(T - tRef), fOn = QM.evolver(eigOn, psiT); psiAt = (t) => (t >= T ? fOff(t - tRef) : fOn(t - T)); }
      }
      const tMax = tMaxUI.values.tMax, S = sys.n > 60 ? 160 : 300;
      const ts = [], probs = [], track = [], ent = [], pops = [];
      const tOp = trackKey ? sys.observables.find((o) => o.key === trackKey) : null;
      for (let s = 0; s <= S; s++) {
        const t = (tMax * s) / S, psi = psiAt(t); ts.push(t); probs.push(V.probs(psi));
        if (tOp) track.push(M.expect(tOp.op, psi).re);
        if (sys.popOps) pops.push(sys.popOps.map((o) => M.expect(o.op, psi).re));
        if (sys.qubits && sys.n <= 64) ent.push(entanglement(psi));
      }
      traj = { ts, probs, track, tOp, ent, tMax, pops };
    }
    function entanglement(psi) {
      const r = QM.rho.reducedFromPure(psi, sys.dims, sys.entSub || 0);
      return { C: sys.dims.length === 2 && sys.dims[0] === 2 && sys.dims[1] === 2 ? QM.concurrencePure(psi) : NaN, S: QM.rho.entropy(r) };
    }

    // drawing -------------------------------------------------------------
    function drawAll(full) {
      const t = transport.t, psi = psiAt(t), p = V.probs(psi);
      drawState(psi, p); drawPops(t); drawLevels(psi, t); drawSubsystems(psi);
      if (full) drawBand();
      meas.update();
      const Hn = sys.Hoff && t > sys.tSwitch ? sys.Hoff : sys.H;
      const E = M.expect(Hn, psi).re, E2 = V.norm2(M.apply(Hn, psi));
      statsUI.set('t', QM.fmt(t, 3)); statsUI.set('E', fmtE(E)); statsUI.set('dE', fmtE(Math.sqrt(Math.max(0, E2 - E * E))));
      statsUI.set('norm', QM.fmt(V.norm2(psi), 6));
      if (sys.qubits) { const e = entanglement(psi); statsUI.set('ent', isNaN(e.C) ? `S = ${QM.fmt(e.S, 3)} bit` : `C = ${QM.fmt(e.C, 3)}, S = ${QM.fmt(e.S, 3)}`); } else statsUI.set('ent', '–');
      ketBox.textContent = sys.n <= 16 ? '|ψ(t)⟩ = ' + QM.ui.ket(psi, sys.labels) : '';
    }
    function drawState(psi, p) {
      const phases = [...Array(psi.n)].map((_, i) => Math.atan2(psi.im[i], psi.re[i]));
      QM.plot.bars(barCanvas, { labels: sys.labels.map((l) => (sys.n <= 16 ? '|' + l + '⟩' : l)), values: Array.from(p), phases, yMax: 1.12, valueLabels: true, ylabel: '|⟨n|ψ⟩|²  (colour = phase)' });
      const lay = sys.layout, show = lay && lay.nodes.length <= 36 && lay.nodes.some((n) => n.idx.length);
      stateCanvas.style.display = show ? '' : 'none';
      if (!show) return;
      const { ctx, w, h } = QM.plot.setup(stateCanvas), th = QM.plot.theme();
      const xs = lay.nodes.map((n) => n.x), ys = lay.nodes.map((n) => n.y);
      const flat = Math.max(...ys) - Math.min(...ys) < 1e-9;
      const R = Math.min(w / 2 - 40, flat ? 1e9 : h / 2 - 40);
      const X = (x) => w / 2 + x * R, Y = (y) => h / 2 - y * R;
      lay.edges.forEach(([a, b, wgt]) => { ctx.strokeStyle = QM.plot.hexA('#6b7280', 0.5); ctx.lineWidth = 1 + 2 * Math.min(1, Math.abs(wgt)); ctx.setLineDash(wgt > 0 ? [] : [5, 3]); ctx.beginPath(); ctx.moveTo(X(xs[a]), Y(ys[a])); ctx.lineTo(X(xs[b]), Y(ys[b])); ctx.stroke(); });
      ctx.setLineDash([]);
      const rMax = Math.max(6, Math.min(34, (flat ? (2 * R) / Math.max(1, lay.nodes.length - 1) : R / Math.max(1.5, lay.nodes.length / 3)) * 0.42));
      lay.nodes.forEach((nd) => {
        let pr = 0, best = -1, ph = 0; nd.idx.forEach((i) => { pr += p[i]; if (p[i] > best) { best = p[i]; ph = Math.atan2(psi.im[i], psi.re[i]); } });
        ctx.strokeStyle = th.axis; ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(X(nd.x), Y(nd.y), rMax, 0, 7); ctx.stroke();
        ctx.fillStyle = QM.plot.phaseColor(ph); ctx.beginPath(); ctx.arc(X(nd.x), Y(nd.y), rMax * Math.sqrt(pr), 0, 7); ctx.fill();
        if (lay.nodes.length <= 16) {
          ctx.fillStyle = th.ink; ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.font = '600 13px system-ui';
          ctx.fillText(nd.label, X(nd.x), Y(nd.y) + rMax + 3);
          ctx.font = '11px system-ui'; ctx.fillStyle = th.muted; ctx.fillText(QM.fmt(pr, 3), X(nd.x), Y(nd.y) + rMax + 18);
        }
      });
      ctx.fillStyle = th.muted; ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.font = '11.5px system-ui';
      ctx.fillText('disc area = probability · colour = phase · solid edge J>0, dashed J<0', 6, 4);
    }
    function drawPops(t) {
      if (!traj) return;
      const big = sys.n > 12;
      const vl = [{ x: t }].concat(sys.tSwitch !== undefined ? [{ x: sys.tSwitch, color: '#c92a2a', dash: [2, 3] }] : []);
      if (sys.popOps) {
        QM.plot.line(popCanvas, { xMin: 0, xMax: traj.tMax, yMin: 0, yMax: 1.02, series: sys.popOps.map((o, k) => ({ x: traj.ts, y: traj.pops.map((r) => r[k]), label: o.label })), vlines: vl, xlabel: 't' });
      } else if (big) {
        const nx = sys.n, ny = traj.ts.length, data = new Float32Array(nx * ny);
        traj.probs.forEach((pr, j) => { for (let i = 0; i < nx; i++) data[j * nx + i] = pr[i]; });
        let mx = 0; traj.probs.forEach((pr) => pr.forEach((x) => (mx = Math.max(mx, x))));
        QM.plot.heatmap(popCanvas, { data, nx, ny, xMin: -0.5, xMax: nx - 0.5, yMin: 0, yMax: traj.tMax, xlabel: 'basis index', ylabel: 't', hline: t, vMax: Math.max(0.02, mx) * 0.5 });
      } else {
        const series = sys.labels.map((l, i) => ({ x: traj.ts, y: traj.probs.map((pr) => pr[i]), label: 'P(' + l + ')' }));
        QM.plot.line(popCanvas, { xMin: 0, xMax: traj.tMax, yMin: 0, yMax: 1.02, series, vlines: vl, xlabel: 't', legend: sys.n <= 8 });
      }
      const series = [];
      if (traj.tOp) series.push({ x: traj.ts, y: traj.track, label: '⟨' + traj.tOp.label.split(/[ (=]/)[0] + '⟩' });
      if (traj.ent.length) {
        if (!isNaN(traj.ent[0].C)) series.push({ x: traj.ts, y: traj.ent.map((e) => e.C), label: 'concurrence C', dash: [6, 3] });
        series.push({ x: traj.ts, y: traj.ent.map((e) => e.S), label: 'entropy S (bit)', dash: [2, 2] });
      }
      trackCanvas.style.display = series.length ? '' : 'none';
      if (series.length) QM.plot.line(trackCanvas, { xMin: 0, xMax: traj.tMax, series, vlines: vl, xlabel: 't' });
    }
    function drawLevels(psi, t) {
      const eig = sys.Hoff && t > sys.tSwitch ? eigOff : eigOn;
      const w8 = eig.vectors.map((v) => C.abs2(V.dot(v, psi)));
      eigTable.innerHTML = '';
      if (eig.values.length > 8) {
        const sx = [], sy = []; eig.values.forEach((E, i) => { sx.push(E, E, E); sy.push(0, w8[i], NaN); });
        QM.plot.line(levelCanvas, { series: [{ x: sx, y: sy, label: '|⟨E|ψ⟩|² (stems at eigenvalues)', width: 2 }], xMin: Math.min(...eig.values) - 1e-6, xMax: Math.max(...eig.values) + 1e-6, yMin: 0, xlabel: 'E' });
        return;
      }
      const { ctx, w, h } = QM.plot.setup(levelCanvas), th = QM.plot.theme();
      const Es = eig.values, lo = Math.min(...Es), hi = Math.max(...Es), span = hi - lo || 1;
      const Y = (E) => h - 18 - ((E - lo) / span) * (h - 36);
      const groups = []; Es.forEach((E, i) => { const g = groups.find((g) => Math.abs(g.E - E) < 1e-6 * Math.max(1, Math.abs(E))); if (g) { g.ids.push(i); g.w += w8[i]; } else groups.push({ E, ids: [i], w: w8[i] }); });
      ctx.font = '12px system-ui'; ctx.textBaseline = 'middle';
      groups.forEach((g) => {
        const y = Y(g.E), x0 = 70, x1 = w - 170;
        ctx.strokeStyle = th.ink; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(x0, y); ctx.lineTo(x1, y); ctx.stroke();
        ctx.fillStyle = QM.plot.hexA('#175ec5', 0.35); ctx.fillRect(x0, y - 4, (x1 - x0) * g.w, 8);
        ctx.fillStyle = th.ink; ctx.textAlign = 'right'; ctx.fillText(fmtE(g.E), x0 - 8, y);
        ctx.textAlign = 'left'; ctx.fillText(`p = ${QM.fmt(g.w, 4)}` + (g.ids.length > 1 ? `  (×${g.ids.length})` : ''), x1 + 8, y);
      });
      ctx.fillStyle = th.muted; ctx.textAlign = 'left'; ctx.textBaseline = 'top';
      ctx.fillText(sys.Hoff ? (t > sys.tSwitch ? 'H₀ (couplings switched off)' : 'H = H₀ + V (couplings on)') : 'E', 6, 2);
      if (sys.n <= 8) {
        const tbl = el('table', { class: 'meas-table' }, el('tr', null, el('th', null, 'E'), el('th', null, 'eigenvector')));
        eig.values.forEach((E, i) => tbl.append(el('tr', null, el('td', { class: 'mono' }, fmtE(E)), el('td', { class: 'mono' }, QM.ui.ket(eig.vectors[i], sys.labels)))));
        eigTable.append(tbl);
      }
    }
    function drawSubsystems(psi) {
      if (!sys.qubits) return;
      blochs.forEach(({ q, b, info }) => {
        const r = QM.rho.reducedFromPure(psi, sys.dims, q.sub), bv = QM.rho.bloch(r);
        b.draw({ vectors: [{ v: bv, label: '' }], axis: sys.axis ? { v: sys.axis, label: 'n' } : null, title: q.label, caption: `|r| = ${QM.fmt(Math.hypot(...bv), 3)} · Tr ρ² = ${QM.fmt(QM.rho.purity(r), 3)}` });
        info.textContent = `ρ = [[${QM.C.fmt(M.get(r, 0, 0), 3)}, ${QM.C.fmt(M.get(r, 0, 1), 3)}], [${QM.C.fmt(M.get(r, 1, 0), 3)}, ${QM.C.fmt(M.get(r, 1, 1), 3)}]]`;
      });
    }
    function drawBand() {
      if (!sys.band) return;
      const { eps, J, gamma } = sys.band, ks = [], E = [], v = [];
      for (let i = 0; i <= 200; i++) { const k = -Math.PI + (2 * Math.PI * i) / 200; ks.push(k); E.push(eps + 2 * J * Math.cos(k) + 2 * gamma * Math.cos(2 * k)); v.push(-2 * J * Math.sin(k) - 4 * gamma * Math.sin(2 * k)); }
      const vl = ['packet', 'plane', 'cat'].includes(initial.mode) ? [{ x: initial.k0 || 0 }] : [];
      QM.plot.line(bandCanvas, { xMin: -Math.PI, xMax: Math.PI, piTicks: true, series: [{ x: ks, y: E, label: 'E(k)' }, { x: ks, y: v, label: 'v(k)', dash: [6, 3] }], vlines: vl, xlabel: 'k' });
    }

    // init
    rebuild(); buildInitUI(); setInitial(); meas.refresh(); drawAll(true);
    cleanups.push(QM.ui.onResize(views, () => drawAll(true)));
    return { destroy() { transport.destroy(); cleanups.forEach((f) => f()); }, redraw: () => drawAll(true) };
  }

  QM.registerEngine('finite', { title: 'Finite-dimensional lab', mount });

  // ------------------------------------------------------------------ tests
  const near = (a, b, tol = 1e-6) => Math.abs(a - b) < tol;
  const NET = (o) => builders.network.build(Object.assign({}, builders.network.defaults, o));
  const SP = (o) => builders.spins.build(Object.assign({}, builders.spins.defaults, o));
  QM.addTest('eigh: random Hermitian 6×6 satisfies H v = λ v and orthonormality', () => {
    const n = 6, H = M.zeros(n);
    for (let i = 0; i < n; i++) for (let j = i; j < n; j++) { const a = Math.random() - 0.5, b = i === j ? 0 : Math.random() - 0.5; M.set(H, i, j, C.of(a, b)); M.set(H, j, i, C.of(a, -b)); }
    const { values, vectors } = QM.eigh(H); let err = 0;
    vectors.forEach((v, k) => { const r = M.apply(H, v); for (let i = 0; i < n; i++) err = Math.max(err, Math.hypot(r.re[i] - values[k] * v.re[i], r.im[i] - values[k] * v.im[i])); vectors.forEach((u, l) => { const d = V.dot(u, v); err = Math.max(err, Math.abs(d.re - (k === l ? 1 : 0)), Math.abs(d.im)); }); });
    return { pass: err < 1e-9, detail: 'max error ' + err.toExponential(2) };
  });
  QM.addTest('Úloha 1 (2018): 3-dot chain from |1⟩ → P(E) = 1/4, 1/2, 1/4; E = ∓√2β, 0', () => {
    const g = QM.spectral(NET({ topology: 'chain', N: 3, J: -1 }).H), p = QM.measureProbs(g, V.basis(3, 0));
    return { pass: g.length === 3 && near(g[0].value, -Math.SQRT2) && near(p[0], 0.25) && near(p[1], 0.5) && near(p[2], 0.25), detail: p.map((x) => QM.fmt(x, 4)).join(', ') };
  });
  QM.addTest('DÚ1 2020: ring α,β,γ – from |B⟩ P = 1/4 (α+γ+2β), 1/4 (α+γ−2β), 1/2 (α−γ); X̂ = |A⟩⟨A|−|C⟩⟨C|', () => {
    const a = 0.3, b = 0.7, c = 0.2, sys = NET({ topology: 'ring', N: 4, eps: a, J: b, gamma: c, labelStyle: 'letters' });
    const g = QM.spectral(sys.H), p = QM.measureProbs(g, V.basis(4, 1));
    const find = (E) => p[g.findIndex((x) => near(x.value, E))];
    const X = sys.observables.find((o) => o.key === 'X').op;
    return { pass: near(find(a + c + 2 * b), 0.25) && near(find(a + c - 2 * b), 0.25) && near(find(a - c), 0.5) && near(X.re[0], 1) && near(X.re[10], -1) && near(X.re[5], 0), detail: g.map((x, i) => `${QM.fmt(x.value, 3)}:${QM.fmt(p[i], 3)}`).join(' ') };
  });
  QM.addTest('Test 2014 Ú2: complete graph K₄ eigenvalues ε−3t (×1), ε+t (×3)', () => {
    const g = QM.spectral(NET({ topology: 'complete', N: 4, eps: 0.5, J: -1 }).H);
    return { pass: g.length === 2 && near(g[0].value, -2.5) && g[1].vectors.length === 3 && near(g[1].value, 1.5), detail: g.map((x) => `${QM.fmt(x.value, 3)}×${x.vectors.length}`).join(', ') };
  });
  QM.addTest('Úloha 3 (2016): ⟨V⟩ = 2β Im(a*b) = β for a=1/√2, b=i/√2', () => {
    const sys = NET({ topology: 'chain', N: 9, J: -1, labelStyle: 'centered' });
    const psi = V.zeros(9); psi.re[4] = Math.SQRT1_2; psi.im[5] = Math.SQRT1_2;
    const v = M.expect(sys.observables.find((o) => o.key === 'V').op, psi).re;
    return { pass: near(v, 1), detail: '⟨V⟩ = ' + QM.fmt(v, 6) };
  });
  QM.addTest('Úloha 2 (2016): chain⊗spin levels ∓√2β ± λ, ±3λ; ground state P(s_y=±½)=½', () => {
    const sys = NET({ topology: 'chain', N: 3, J: -1, spin: true, lambda: 0.1 });
    const eig = QM.eigh(sys.H), want = [-Math.SQRT2 - 0.1, -Math.SQRT2 + 0.1, -0.3, 0.3, Math.SQRT2 - 0.1, Math.SQRT2 + 0.1];
    const g = QM.spectral(sys.observables.find((o) => o.key === 'sy').op), p = QM.measureProbs(g, eig.vectors[0]);
    return { pass: want.every((E, i) => near(E, eig.values[i])) && near(p[0], 0.5) && near(p[1], 0.5), detail: eig.values.map((x) => QM.fmt(x, 3)).join(', ') };
  });
  QM.addTest('Network: long clean chain analytic eigenbasis equals numerical (N=40 vs Jacobi)', () => {
    const a = NET({ topology: 'chain', N: 40, J: -1 }), b = QM.eigh(a.H);
    return { pass: a.eig.values.every((E, i) => near(E, b.values[i], 1e-9)), detail: 'E₀ = ' + QM.fmt(a.eig.values[0], 6) };
  });
  QM.addTest('Entanglement Úloha 3: |10⟩ → concurrence 1 at T = π/(8ω); full swap at π/(4ω)', () => {
    const sys = SP({ e: 1, D: 1 });
    const f = QM.evolver(QM.eigh(sys.H), V.basis(4, 2));
    const C1 = QM.concurrencePure(f(Math.PI / 8)), p01 = V.probs(f(Math.PI / 4))[1];
    return { pass: near(C1, 1) && near(p01, 1), detail: `C(π/8) = ${QM.fmt(C1, 6)}, P(01, π/4) = ${QM.fmt(p01, 6)}` };
  });
  QM.addTest('Úloha 2: A = Sx² − Sy² on |x+x+⟩ → P(+1)=½, P(0)=½, P(−1)=0; singlet → P(0)=1', () => {
    const g = QM.spectral(SP({}).observables.find((o) => o.key === 'A').op);
    const p = QM.measureProbs(g, V.from([0.5, 0.5, 0.5, 0.5])), ps = QM.measureProbs(g, V.from([0, Math.SQRT1_2, -Math.SQRT1_2, 0]));
    const at = (pp, val) => pp[g.findIndex((x) => near(x.value, val))] || 0;
    return { pass: near(at(p, 1), 0.5) && near(at(p, 0), 0.5) && near(at(p, -1), 0) && near(at(ps, 0), 1), detail: g.map((x, i) => `${QM.fmt(x.value, 2)}:${QM.fmt(p[i], 3)}`).join(' ') };
  });
  QM.addTest('Test 2014 Ú3: |++⟩, S along (1,1,1)/√3 → P(0) = 1/3', () => {
    const sys = SP({ thN: Math.acos(1 / Math.sqrt(3)), phN: Math.PI / 4 });
    const g = QM.spectral(sys.observables.find((o) => o.key === 'Sn').op), p = QM.measureProbs(g, V.basis(4, 0));
    const p0 = p[g.findIndex((x) => near(x.value, 0))];
    return { pass: near(p0, 1 / 3), detail: 'P(0) = ' + QM.fmt(p0, 6) };
  });
  QM.addTest('Spins: 3-spin GHZ has reduced entropy 1 bit; Heisenberg chain conserves total S_z', () => {
    const sys = SP({ Ns: 3, Jxx: 1, Jyy: 1, Jzz: 0.5 }), ghz = V.zeros(8); ghz.re[0] = ghz.re[7] = Math.SQRT1_2;
    const S = QM.rho.entropy(QM.rho.reducedFromPure(ghz, sys.dims, 0));
    const Sz = sys.observables.find((o) => o.key === 'Sz').op, cm = M.maxAbs(M.comm(sys.H, Sz));
    return { pass: near(S, 1) && cm < 1e-12, detail: `S = ${QM.fmt(S, 6)}, ‖[H,S_z]‖ = ${cm.toExponential(1)}` };
  });
  QM.addTest('Úloha 5: 2p, j=½, m=½ → ρ_spin = diag(1/3, 2/3), P(+ along (1,0,1)/√2) ≈ 0.382', () => {
    const sys = builders.orbitalSpin.build({ l: 1, xi: 0, bz: 0, thN: Math.PI / 4, phN: 0 }), psi = cgState(1, 0.5, 0.5);
    const r = QM.rho.reducedFromPure(psi, sys.dims, 1);
    const g = QM.spectral(sys.observables.find((o) => o.key === 'Sn').op), p = QM.measureProbs(g, psi);
    const J2 = M.expect(sys.observables.find((o) => o.key === 'J2').op, psi).re;
    return { pass: near(r.re[0], 1 / 3) && near(r.re[3], 2 / 3) && near(p[1], 0.5 * (1 - 1 / (3 * Math.SQRT2))) && near(J2, 0.75), detail: `ρ↑↑=${QM.fmt(r.re[0], 4)}, P(+)=${QM.fmt(p[1], 4)}, ⟨J²⟩=${QM.fmt(J2, 4)}` };
  });
  QM.addTest('U1 (Nobel): p₊ = ½ + (√3/4)(|α|²−|β|²) for θ = π/3', () => {
    const sys = SP({ thN: Math.PI / 3 }), a = Math.cos(0.4), b = Math.sin(0.4), s = Math.SQRT1_2;
    const psi = V.from([a * s, b * s, a * s, -b * s]);
    const g = QM.spectral(sys.observables.find((o) => o.key === 'Sn1').op), p = QM.measureProbs(g, psi);
    const want = 0.5 + (Math.sqrt(3) / 4) * (a * a - b * b);
    return { pass: near(p[1], want), detail: `p₊ = ${QM.fmt(p[1], 5)} (expected ${QM.fmt(want, 5)})` };
  });
})();
