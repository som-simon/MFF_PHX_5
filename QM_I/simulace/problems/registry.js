/* Catalogue: topics, source files (tracking which PDFs were processed) and the
   problem list. Problems themselves are added per engine in problems/<engine>.js
   via QM.addProblems([...]).

   Problem schema:
   { id, topic, title, source (file in ../Problems), sheet (label in the PDF),
     engine (key in QM.engines), statement (HTML, $tex$ ok), solution (HTML),
     preset (engine-specific object), sandbox?: true }                       */
(function () {
  'use strict';
  const QM = (window.QM = window.QM || {});

  const topics = [
    { id: 'spin', title: 'Spin ½ & measurement', desc: 'Single qubit on the Bloch sphere: observables $\\hat S_{\\vec n}$, Born rule, precession, mixed states.' },
    { id: 'composite', title: 'Spins, qubits & entanglement', desc: 'Chain of 1–6 spins/qubits with fields and couplings, or orbital $l$ ⊗ spin ½: joint measurements, reduced density matrices, concurrence, entanglement entropy.' },
    { id: 'dots', title: 'Quantum-dot networks', desc: 'N dots as a chain, ring, complete graph, star or 2D grid; hopping, next-nearest hopping, tilt, disorder, optional spin ½. Symmetry operators, velocity, band structure.' },
    { id: 'ring', title: 'Particle on a ring', desc: 'Free rotor, rotor in a potential (pendulum) or with magnetic flux; $L$ eigenbasis, ladder operators $e^{\\pm i\\varphi}$, revivals, quantum carpet, Dirichlet kernel.' },
    { id: 'wave1d', title: '1D wavefunctions', desc: 'General $H=p^2/2m+V(x)$: free, harmonic, δ-well, finite well, box, barrier, step, double well, linear, periodic. Packets, cat states, eigenstates, coherent states; $x$ and $p$ pictures, levels, scattering.' },
    { id: 'angular', title: 'Angular momentum', desc: 'Decomposition of angular functions into $Y_{lm}$, measurement of $\\hat L^2$, $\\hat L_z$.' },
  ];

  /* Every file in QM_I/Problems must be listed here (status: done | partial | no-sim | todo).
     Keep PROCESSED.md in sync. */
  const sources = [
    { file: '3-Uloha_3_entanglement_Qbitu.pdf', title: 'Úloha 3 – entanglement of two qubits', status: 'done', problems: ['entangle-qubits'] },
    { file: 'DU1_2020.pdf', title: 'DÚ 1 (2020) – quantum dots (4-dot ring)', status: 'done', problems: ['dots-ring4'] },
    { file: 'QM-II-2015-du1.pdf', title: 'QM II 2015 DÚ 1 – spin polarisation in a magnetic field', status: 'done', problems: ['spin-precession'] },
    { file: 'QMIcv_09.30_zadani.pdf', title: 'Cvičení 1 – linear algebra, Pauli matrices', status: 'partial', problems: ['pauli-decomp', 'pauli-eigen'], note: 'Ú1–4, Ú7 are proofs – no simulation.' },
    { file: 'QMIcv_10.07_zadani.pdf', title: 'Cvičení 2 – spin ½ along a general axis', status: 'done', problems: ['spin-axis'] },
    { file: 'reseni_kvantovy_rotor.pdf', title: 'Solution – simple quantum rotor (Cvičení 5)', status: 'done', problems: ['rotor-basics', 'dirichlet'] },
    { file: 'U1.pdf', title: 'Úloha 1 – tribute to the 2012 Nobel prize (spin ⊗ cavity)', status: 'done', problems: ['nobel-cavity'] },
    { file: 'Uloha1_2018.pdf', title: 'Úloha 1 (2018) – three-dot chain', status: 'done', problems: ['dots-chain3'] },
    { file: 'Uloha2.pdf', title: 'Úloha 2 – two spin-½ particles, $A=(\\omega/\\hbar)(S_x^2-S_y^2)$', status: 'done', problems: ['two-spins-A'] },
    { file: 'Uloha2_2016.pdf', title: 'Úloha 2 (2016) – three-dot chain with spin', status: 'done', problems: ['dots-chain3-spin'] },
    { file: 'Uloha3.pdf', title: 'Úloha 3 – Gaussian wave packets', status: 'done', problems: ['gauss-packet'] },
    { file: 'Uloha3_2016.pdf', title: 'Úloha 3 (2016) – infinite quantum chain, velocity operator', status: 'done', problems: ['chain-infinite'] },
    { file: 'Uloha3_2018.pdf', title: 'Úloha 3 (2018) – quantum rotor, ladder operators', status: 'done', problems: ['rotor-ladder'] },
    { file: 'UlohaV.pdf', title: 'Úloha 5 – reduced density matrix (2p, j=½)', status: 'done', problems: ['reduced-2p'] },
    { file: 'Vzorove_reseni_2014.pdf', title: 'Sample solutions 2014 (5 problems)', status: 'done', problems: ['displaced-ho', 'dots-k4', 'two-spins-e', 'spherical-well-L2', 'delta-well-p'] },
    { file: 'Dirichletovo jádro a periodická Diracova delta distribuce.htm', title: 'Old standalone sim – Dirichlet kernel', status: 'done', problems: ['dirichlet'], note: 'Ported into the ring engine.' },
    { file: 'Kvantové tečky – časový vývoj.htm', title: 'Old standalone sim – 4 quantum dots', status: 'done', problems: ['dots-ring4'], note: 'Replaced by the finite engine.' },
    { file: 'Kvantový rotor — vlastní stavy a časový vývoj.htm', title: 'Old standalone sim – quantum rotor', status: 'done', problems: ['rotor-basics'], note: 'Ported into the ring engine.' },
  ];

  const problems = [];
  QM.registry = { topics, sources, problems };
  QM.addProblems = (list) => problems.push(...list);
  QM.engines = QM.engines || {};
  QM.registerEngine = (name, def) => { QM.engines[name] = def; };
  QM.tests = QM.tests || [];
  QM.addTest = (name, fn) => QM.tests.push({ name, fn });
})();
