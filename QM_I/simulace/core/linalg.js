/* Complex numbers & linear algebra for small Hilbert spaces (ħ = 1 everywhere).
   Scalar: {re, im}.  Vector: {n, re: Float64Array, im: Float64Array}.
   Matrix: {n, re: Float64Array(n*n), im: Float64Array(n*n)}, row-major (i*n + j). */
(function () {
  'use strict';
  const QM = (window.QM = window.QM || {});

  // ---------- scalars ----------
  const C = {
    of: (re, im = 0) => ({ re, im }),
    add: (a, b) => ({ re: a.re + b.re, im: a.im + b.im }),
    sub: (a, b) => ({ re: a.re - b.re, im: a.im - b.im }),
    mul: (a, b) => ({ re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re }),
    conj: (a) => ({ re: a.re, im: -a.im }),
    scale: (a, s) => ({ re: a.re * s, im: a.im * s }),
    abs2: (a) => a.re * a.re + a.im * a.im,
    abs: (a) => Math.hypot(a.re, a.im),
    arg: (a) => Math.atan2(a.im, a.re),
    expi: (t) => ({ re: Math.cos(t), im: Math.sin(t) }),
    polar: (r, t) => ({ re: r * Math.cos(t), im: r * Math.sin(t) }),
    div: (a, b) => { const d = b.re * b.re + b.im * b.im; return { re: (a.re * b.re + a.im * b.im) / d, im: (a.im * b.re - a.re * b.im) / d }; },
    fmt(a, digits = 3) {
      const r = +a.re.toFixed(digits), i = +a.im.toFixed(digits);
      if (Math.abs(i) < Math.pow(10, -digits)) return fmtNum(r, digits);
      if (Math.abs(r) < Math.pow(10, -digits)) return fmtNum(i, digits) + 'i';
      return fmtNum(r, digits) + (i < 0 ? ' − ' : ' + ') + fmtNum(Math.abs(i), digits) + 'i';
    },
  };
  function fmtNum(x, digits = 3) {
    if (!isFinite(x)) return String(x);
    let s = (+x.toFixed(digits)).toString();
    if (s === '-0') s = '0';
    return s.replace('-', '−');
  }
  QM.C = C;
  QM.fmt = fmtNum;

  // ---------- vectors ----------
  const V = {
    zeros: (n) => ({ n, re: new Float64Array(n), im: new Float64Array(n) }),
    /** from array of numbers or {re,im} or [re,im] */
    from(arr) {
      const v = V.zeros(arr.length);
      arr.forEach((x, i) => {
        if (typeof x === 'number') v.re[i] = x;
        else if (Array.isArray(x)) { v.re[i] = x[0]; v.im[i] = x[1] || 0; }
        else { v.re[i] = x.re; v.im[i] = x.im || 0; }
      });
      return v;
    },
    basis(n, k) { const v = V.zeros(n); v.re[k] = 1; return v; },
    clone: (v) => ({ n: v.n, re: Float64Array.from(v.re), im: Float64Array.from(v.im) }),
    get: (v, i) => ({ re: v.re[i], im: v.im[i] }),
    set(v, i, z) { v.re[i] = z.re; v.im[i] = z.im; },
    norm2(v) { let s = 0; for (let i = 0; i < v.n; i++) s += v.re[i] * v.re[i] + v.im[i] * v.im[i]; return s; },
    normalize(v) { const s = Math.sqrt(V.norm2(v)); if (s > 0) for (let i = 0; i < v.n; i++) { v.re[i] /= s; v.im[i] /= s; } return v; },
    /** <a|b> */
    dot(a, b) { let re = 0, im = 0; for (let i = 0; i < a.n; i++) { re += a.re[i] * b.re[i] + a.im[i] * b.im[i]; im += a.re[i] * b.im[i] - a.im[i] * b.re[i]; } return { re, im }; },
    scale(v, z) { const w = V.zeros(v.n); for (let i = 0; i < v.n; i++) { w.re[i] = v.re[i] * z.re - v.im[i] * z.im; w.im[i] = v.re[i] * z.im + v.im[i] * z.re; } return w; },
    addTo(acc, v, z) { // acc += z v
      for (let i = 0; i < v.n; i++) { acc.re[i] += v.re[i] * z.re - v.im[i] * z.im; acc.im[i] += v.re[i] * z.im + v.im[i] * z.re; }
      return acc;
    },
    probs(v) { const p = new Float64Array(v.n); for (let i = 0; i < v.n; i++) p[i] = v.re[i] * v.re[i] + v.im[i] * v.im[i]; return p; },
    kron(a, b) { const v = V.zeros(a.n * b.n); for (let i = 0; i < a.n; i++) for (let j = 0; j < b.n; j++) { const k = i * b.n + j; v.re[k] = a.re[i] * b.re[j] - a.im[i] * b.im[j]; v.im[k] = a.re[i] * b.im[j] + a.im[i] * b.re[j]; } return v; },
  };
  QM.V = V;

  // ---------- matrices ----------
  const M = {
    zeros: (n) => ({ n, re: new Float64Array(n * n), im: new Float64Array(n * n) }),
    eye(n) { const m = M.zeros(n); for (let i = 0; i < n; i++) m.re[i * n + i] = 1; return m; },
    /** rows: array of arrays of number | [re,im] | {re,im} */
    from(rows) {
      const n = rows.length, m = M.zeros(n);
      rows.forEach((row, i) => row.forEach((x, j) => {
        const k = i * n + j;
        if (typeof x === 'number') m.re[k] = x;
        else if (Array.isArray(x)) { m.re[k] = x[0]; m.im[k] = x[1] || 0; }
        else { m.re[k] = x.re; m.im[k] = x.im || 0; }
      }));
      return m;
    },
    diag(arr) { const n = arr.length, m = M.zeros(n); arr.forEach((x, i) => (m.re[i * n + i] = x)); return m; },
    clone: (m) => ({ n: m.n, re: Float64Array.from(m.re), im: Float64Array.from(m.im) }),
    get: (m, i, j) => ({ re: m.re[i * m.n + j], im: m.im[i * m.n + j] }),
    set(m, i, j, z) { m.re[i * m.n + j] = typeof z === 'number' ? z : z.re; m.im[i * m.n + j] = typeof z === 'number' ? 0 : z.im; },
    add(...ms) { const r = M.zeros(ms[0].n); for (const m of ms) for (let k = 0; k < r.re.length; k++) { r.re[k] += m.re[k]; r.im[k] += m.im[k]; } return r; },
    sub(a, b) { const r = M.clone(a); for (let k = 0; k < r.re.length; k++) { r.re[k] -= b.re[k]; r.im[k] -= b.im[k]; } return r; },
    /** scale by real number or complex {re,im} */
    scale(m, z) {
      if (typeof z === 'number') z = { re: z, im: 0 };
      const r = M.zeros(m.n);
      for (let k = 0; k < r.re.length; k++) { r.re[k] = m.re[k] * z.re - m.im[k] * z.im; r.im[k] = m.re[k] * z.im + m.im[k] * z.re; }
      return r;
    },
    mul(a, b) {
      const n = a.n, r = M.zeros(n);
      for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) {
        const ar = a.re[i * n + k], ai = a.im[i * n + k];
        if (ar === 0 && ai === 0) continue;
        for (let j = 0; j < n; j++) {
          const br = b.re[k * n + j], bi = b.im[k * n + j];
          r.re[i * n + j] += ar * br - ai * bi; r.im[i * n + j] += ar * bi + ai * br;
        }
      }
      return r;
    },
    dagger(m) { const n = m.n, r = M.zeros(n); for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { r.re[j * n + i] = m.re[i * n + j]; r.im[j * n + i] = -m.im[i * n + j]; } return r; },
    comm: (a, b) => M.sub(M.mul(a, b), M.mul(b, a)),
    apply(m, v) {
      const n = m.n, w = V.zeros(n);
      for (let i = 0; i < n; i++) { let re = 0, im = 0; for (let j = 0; j < n; j++) { const ar = m.re[i * n + j], ai = m.im[i * n + j]; re += ar * v.re[j] - ai * v.im[j]; im += ar * v.im[j] + ai * v.re[j]; } w.re[i] = re; w.im[i] = im; }
      return w;
    },
    /** <v|M|v> (real part if Hermitian) */
    expect(m, v) { return V.dot(v, M.apply(m, v)); },
    trace(m) { let re = 0, im = 0; for (let i = 0; i < m.n; i++) { re += m.re[i * m.n + i]; im += m.im[i * m.n + i]; } return { re, im }; },
    kron(a, b) {
      const n = a.n * b.n, r = M.zeros(n);
      for (let i1 = 0; i1 < a.n; i1++) for (let j1 = 0; j1 < a.n; j1++) {
        const ar = a.re[i1 * a.n + j1], ai = a.im[i1 * a.n + j1];
        if (ar === 0 && ai === 0) continue;
        for (let i2 = 0; i2 < b.n; i2++) for (let j2 = 0; j2 < b.n; j2++) {
          const br = b.re[i2 * b.n + j2], bi = b.im[i2 * b.n + j2];
          const k = (i1 * b.n + i2) * n + (j1 * b.n + j2);
          r.re[k] = ar * br - ai * bi; r.im[k] = ar * bi + ai * br;
        }
      }
      return r;
    },
    /** kron of a list */
    kronAll(list) { return list.reduce((acc, m) => M.kron(acc, m)); },
    outer(a, b) { // |a><b|
      const n = a.n, r = M.zeros(n);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { r.re[i * n + j] = a.re[i] * b.re[j] + a.im[i] * b.im[j]; r.im[i * n + j] = a.im[i] * b.re[j] - a.re[i] * b.im[j]; }
      return r;
    },
    isHermitian(m, tol = 1e-9) { const n = m.n; for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { if (Math.abs(m.re[i * n + j] - m.re[j * n + i]) > tol || Math.abs(m.im[i * n + j] + m.im[j * n + i]) > tol) return false; } return true; },
    maxAbs(m) { let s = 0; for (let k = 0; k < m.re.length; k++) s = Math.max(s, Math.hypot(m.re[k], m.im[k])); return s; },
  };
  QM.M = M;

  // ---------- Hermitian eigensolver: cyclic complex Jacobi ----------
  /** returns {values: number[] ascending, vectors: Vector[]} */
  function eigh(H) {
    const n = H.n;
    const are = Float64Array.from(H.re), aim = Float64Array.from(H.im);
    // symmetrize (guard against round-off)
    for (let i = 0; i < n; i++) { aim[i * n + i] = 0; for (let j = i + 1; j < n; j++) { const r = 0.5 * (are[i * n + j] + are[j * n + i]), m = 0.5 * (aim[i * n + j] - aim[j * n + i]); are[i * n + j] = r; are[j * n + i] = r; aim[i * n + j] = m; aim[j * n + i] = -m; } }
    const vre = new Float64Array(n * n), vim = new Float64Array(n * n);
    for (let i = 0; i < n; i++) vre[i * n + i] = 1;
    let scale = 0; for (let k = 0; k < n * n; k++) scale = Math.max(scale, Math.abs(are[k]), Math.abs(aim[k]));
    const eps = 1e-14 * (scale || 1);
    for (let sweep = 0; sweep < 100; sweep++) {
      let off = 0;
      for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += are[p * n + q] ** 2 + aim[p * n + q] ** 2;
      if (Math.sqrt(off) < eps) break;
      for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
        const xr = are[p * n + q], xi = aim[p * n + q];
        const r = Math.hypot(xr, xi);
        if (r < eps * 1e-3) continue;
        const er = xr / r, ei = xi / r; // phase e
        const app = are[p * n + p], aqq = are[q * n + q];
        const theta = (aqq - app) / (2 * r);
        const t = (theta >= 0 ? 1 : -1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        // U block (rows p,q ; cols p,q): [[c, s], [-s*conj(e), c*conj(e)]]
        const upp_r = c, upp_i = 0, upq_r = s, upq_i = 0;
        const uqp_r = -s * er, uqp_i = s * ei;
        const uqq_r = c * er, uqq_i = -c * ei;
        // A <- A U (columns p,q)
        for (let k = 0; k < n; k++) {
          const kp = k * n + p, kq = k * n + q;
          const apr = are[kp], api = aim[kp], aqr = are[kq], aqi = aim[kq];
          are[kp] = apr * upp_r - api * upp_i + aqr * uqp_r - aqi * uqp_i;
          aim[kp] = apr * upp_i + api * upp_r + aqr * uqp_i + aqi * uqp_r;
          are[kq] = apr * upq_r - api * upq_i + aqr * uqq_r - aqi * uqq_i;
          aim[kq] = apr * upq_i + api * upq_r + aqr * uqq_i + aqi * uqq_r;
          const vpr = vre[kp], vpi = vim[kp], vqr = vre[kq], vqi = vim[kq];
          vre[kp] = vpr * upp_r - vpi * upp_i + vqr * uqp_r - vqi * uqp_i;
          vim[kp] = vpr * upp_i + vpi * upp_r + vqr * uqp_i + vqi * uqp_r;
          vre[kq] = vpr * upq_r - vpi * upq_i + vqr * uqq_r - vqi * uqq_i;
          vim[kq] = vpr * upq_i + vpi * upq_r + vqr * uqq_i + vqi * uqq_r;
        }
        // A <- U† A (rows p,q)
        for (let k = 0; k < n; k++) {
          const pk = p * n + k, qk = q * n + k;
          const apr = are[pk], api = aim[pk], aqr = are[qk], aqi = aim[qk];
          // conj(u_pp) A_pk + conj(u_qp) A_qk
          are[pk] = upp_r * apr + upp_i * api + uqp_r * aqr + uqp_i * aqi;
          aim[pk] = upp_r * api - upp_i * apr + uqp_r * aqi - uqp_i * aqr;
          are[qk] = upq_r * apr + upq_i * api + uqq_r * aqr + uqq_i * aqi;
          aim[qk] = upq_r * api - upq_i * apr + uqq_r * aqi - uqq_i * aqr;
        }
        are[p * n + q] = are[q * n + p] = 0; aim[p * n + q] = aim[q * n + p] = 0;
        aim[p * n + p] = aim[q * n + q] = 0;
      }
    }
    const idx = [...Array(n).keys()].sort((a, b) => are[a * n + a] - are[b * n + b]);
    const values = idx.map((i) => are[i * n + i]);
    const vectors = idx.map((j) => {
      const v = V.zeros(n);
      for (let i = 0; i < n; i++) { v.re[i] = vre[i * n + j]; v.im[i] = vim[i * n + j]; }
      // fix global phase: make largest component real positive (nicer display)
      let best = 0, bk = 0; for (let i = 0; i < n; i++) { const a = v.re[i] ** 2 + v.im[i] ** 2; if (a > best + 1e-12) { best = a; bk = i; } }
      const ph = Math.atan2(v.im[bk], v.re[bk]);
      return V.scale(v, C.expi(-ph));
    });
    return { values, vectors };
  }

  /** group eigenpairs into eigenspaces: [{value, vectors, P}] */
  function spectral(A, tol = 1e-7) {
    const { values, vectors } = eigh(A);
    const groups = [];
    values.forEach((val, i) => {
      const g = groups.length && Math.abs(groups[groups.length - 1].value - val) < tol * Math.max(1, Math.abs(val)) ? groups[groups.length - 1] : null;
      if (g) { g.vectors.push(vectors[i]); g._sum += val; g.value = g._sum / g.vectors.length; }
      else groups.push({ value: val, _sum: val, vectors: [vectors[i]] });
    });
    groups.forEach((g) => { g.P = g.vectors.reduce((acc, v) => M.add(acc, M.outer(v, v)), M.zeros(A.n)); delete g._sum; });
    return groups;
  }

  /** state at time t from eigendecomposition {values, vectors} (ħ=1) */
  function evolve(eig, psi0, t) {
    const out = V.zeros(psi0.n);
    eig.vectors.forEach((v, k) => { const c = C.mul(C.expi(-eig.values[k] * t), V.dot(v, psi0)); V.addTo(out, v, c); });
    return out;
  }
  /** precompute coefficients for fast repeated evolution */
  function evolver(eig, psi0) {
    const coeffs = eig.vectors.map((v) => V.dot(v, psi0));
    return (t) => { const out = V.zeros(psi0.n); eig.vectors.forEach((v, k) => V.addTo(out, v, C.mul(C.expi(-eig.values[k] * t), coeffs[k]))); return out; };
  }

  // ---------- density matrices ----------
  const rho = {
    pure: (v) => M.outer(v, v),
    /** partial trace keeping subsystem `keep` of dims array */
    partial(r, dims, keep) {
      const D = dims.reduce((a, b) => a * b, 1), dk = dims[keep];
      const before = dims.slice(0, keep).reduce((a, b) => a * b, 1), after = dims.slice(keep + 1).reduce((a, b) => a * b, 1);
      const out = M.zeros(dk);
      for (let a = 0; a < dk; a++) for (let b = 0; b < dk; b++) {
        let re = 0, im = 0;
        for (let x = 0; x < before; x++) for (let y = 0; y < after; y++) {
          const i = (x * dk + a) * after + y, j = (x * dk + b) * after + y;
          re += r.re[i * D + j]; im += r.im[i * D + j];
        }
        out.re[a * dk + b] = re; out.im[a * dk + b] = im;
      }
      return out;
    },
    /** reduced density matrix of subsystem `keep` directly from a pure state */
    reducedFromPure(v, dims, keep) { return rho.partial(M.outer(v, v), dims, keep); },
    purity(r) { return M.trace(M.mul(r, r)).re; },
    entropy(r) { // von Neumann, log2
      const { values } = eigh(r); let s = 0; values.forEach((p) => { if (p > 1e-12) s -= p * Math.log2(p); }); return s;
    },
    bloch(r2) { // 2x2 -> [x,y,z]
      return [2 * r2.re[1], -2 * r2.im[1], r2.re[0] - r2.re[3]];
    },
  };
  QM.rho = rho;

  /** concurrence of a pure two-qubit state (a|00>+b|01>+c|10>+d|11>) */
  function concurrencePure(v) {
    const ad = C.mul(V.get(v, 0), V.get(v, 3)), bc = C.mul(V.get(v, 1), V.get(v, 2));
    return 2 * C.abs(C.sub(ad, bc));
  }

  // ---------- standard operators ----------
  const ops = {
    sx: () => M.from([[0, 1], [1, 0]]),
    sy: () => M.from([[0, [0, -1]], [[0, 1], 0]]),
    sz: () => M.from([[1, 0], [0, -1]]),
    /** σ·n with n=(sinθcosφ, sinθsinφ, cosθ) */
    sn(theta, phi) {
      return M.add(M.scale(ops.sx(), Math.sin(theta) * Math.cos(phi)), M.scale(ops.sy(), Math.sin(theta) * Math.sin(phi)), M.scale(ops.sz(), Math.cos(theta)));
    },
    /** spin-j matrices (ħ=1), basis m = j, j-1, ..., -j */
    spin(j) {
      const d = Math.round(2 * j + 1);
      const Jz = M.zeros(d), Jp = M.zeros(d);
      for (let k = 0; k < d; k++) { const m = j - k; Jz.re[k * d + k] = m; if (k > 0) Jp.re[(k - 1) * d + k] = Math.sqrt(j * (j + 1) - m * (m + 1)); }
      const Jm = M.dagger(Jp);
      const Jx = M.scale(M.add(Jp, Jm), 0.5);
      const Jy = M.scale(M.sub(Jp, Jm), { re: 0, im: -0.5 });
      return { Jx, Jy, Jz, Jp, Jm, d };
    },
    /** embed single-site operator op on subsystem k of dims */
    embed(op, dims, k) { return M.kronAll(dims.map((d, i) => (i === k ? op : M.eye(d)))); },
    proj: (v) => M.outer(v, v),
  };
  QM.ops = ops;

  /** Born rule & collapse: groups from spectral(); returns probs array */
  function measureProbs(groups, psi) { return groups.map((g) => g.vectors.reduce((s, v) => s + C.abs2(V.dot(v, psi)), 0)); }
  function collapse(group, psi) { return V.normalize(M.apply(group.P, psi)); }
  function sample(probs) { let r = Math.random(), acc = 0; for (let i = 0; i < probs.length; i++) { acc += probs[i]; if (r < acc) return i; } return probs.length - 1; }

  Object.assign(QM, { eigh, spectral, evolve, evolver, concurrencePure, measureProbs, collapse, sample });
})();
