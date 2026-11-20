/* Problems using the 'finite' engine (quantum dots, chains, qubits, orbital⊗spin). */
(function () {
  'use strict';
  const S = Math.SQRT1_2;
  const pi = (x) => QM.fmt(x / Math.PI, 3) + 'π';

  QM.addProblems([
    // ---------------------------------------------------------------- dots
    {
      id: 'dots-ring4', topic: 'dots', engine: 'finite', source: 'DU1_2020.pdf', sheet: 'DÚ 1 (2020) – Úloha 1',
      title: 'Four quantum dots on a ring',
      statement: `A particle tunnels between four dots $A,B,C,D$ placed at $\\pm x$, $\\pm y$. In the basis $(|A\\rangle,|B\\rangle,|C\\rangle,|D\\rangle)$
        $$\\hat H=\\begin{pmatrix}\\alpha&\\beta&\\gamma&\\beta\\\\\\beta&\\alpha&\\beta&\\gamma\\\\\\gamma&\\beta&\\alpha&\\beta\\\\\\beta&\\gamma&\\beta&\\alpha\\end{pmatrix},\\quad
        \\hat X=|A\\rangle\\langle A|-|C\\rangle\\langle C|,\\ \\hat Y=|B\\rangle\\langle B|-|D\\rangle\\langle D|,$$
        $\\hat P_X=|A\\rangle\\langle C|+|C\\rangle\\langle A|$, $\\hat P_Y=|B\\rangle\\langle D|+|D\\rangle\\langle B|$.
        <ol><li>Find common eigenvectors of $\\hat P_X$, $\\hat P_Y$.</li><li>Write $\\hat H$ in that basis and diagonalise it.</li>
        <li>$\\hat X$ is measured with result 0, then $\\hat Y$ with result 1. Which energies can a subsequent energy measurement give, with which probabilities?</li></ol>
        <p class="note">Mapping: $\\varepsilon=\\alpha$, $J=\\beta$, $\\gamma=\\gamma$. Try: evolve, measure X̂ until you get 0, then Ŷ until 1, then Ĥ.</p>`,
      solution: `<ol><li>$(|A\\rangle\\pm|C\\rangle)/\\sqrt2$ with $(p_X,p_Y)=(\\pm1,0)$; $(|B\\rangle\\pm|D\\rangle)/\\sqrt2$ with $(0,\\pm1)$.</li>
        <li>Block form: $\\begin{pmatrix}\\alpha+\\gamma&2\\beta\\\\2\\beta&\\alpha+\\gamma\\end{pmatrix}$ on the symmetric pair, $\\alpha-\\gamma$ on both antisymmetric vectors.
        $E_1=\\alpha+\\gamma+2\\beta$: $\\tfrac12(1,1,1,1)$; $E_2=\\alpha+\\gamma-2\\beta$: $\\tfrac12(1,-1,1,-1)$; $E_3=\\alpha-\\gamma$ (twice): $(|A\\rangle-|C\\rangle)/\\sqrt2$, $(|B\\rangle-|D\\rangle)/\\sqrt2$.</li>
        <li>After the measurements the state is $|B\\rangle$: $P(\\alpha+\\gamma+2\\beta)=\\tfrac14$, $P(\\alpha+\\gamma-2\\beta)=\\tfrac14$, $P(\\alpha-\\gamma)=\\tfrac12$.</li></ol>`,
      preset: {
        builder: 'network', values: { topology: 'ring', N: 4, eps: 0, J: 1, gamma: 0.5, labelStyle: 'letters' },
        initial: { mode: 'site', k: 0 }, observables: ['X', 'Yc', 'H', 'PX', 'PY', 'site'], track: ['X'], tMax: 2 * Math.PI,
        statePresets: [{ label: '|A⟩', vector: [1, 0, 0, 0] }, { label: '(|A⟩ + i|B⟩)/√2', vector: [S, [0, S], 0, 0] }, { label: 'uniform ½(1,1,1,1)', vector: [0.5, 0.5, 0.5, 0.5] }],
      },
    },
    {
      id: 'dots-chain3', topic: 'dots', engine: 'finite', source: 'Uloha1_2018.pdf', sheet: 'Úloha 1 (2018)',
      title: 'Three-dot chain and mirror symmetry',
      statement: `$\\mathcal H=\\mathrm{span}(|1\\rangle,|2\\rangle,|3\\rangle)$, $\\hat H=-\\beta(|1\\rangle\\langle2|+|2\\rangle\\langle1|+|2\\rangle\\langle3|+|3\\rangle\\langle2|)$, $\\beta>0$; $\\hat Y=|2\\rangle\\langle2|+|1\\rangle\\langle3|+|3\\rangle\\langle1|$.
        <ol><li>Spectral decomposition of $\\hat H$.</li><li>Is $\\hat Y$ an observable? Meaning, eigenvectors, $[\\hat H,\\hat Y]$?</li><li>The particle starts in $|1\\rangle$: possible energies and their probabilities.</li></ol>
        <p class="note">Mapping: $J=-\\beta$ (default β = 1).</p>`,
      solution: `<ol><li>$E=-\\sqrt2\\beta$: $\\tfrac12(1,\\sqrt2,1)$; $E=0$: $\\tfrac1{\\sqrt2}(1,0,-1)$; $E=+\\sqrt2\\beta$: $\\tfrac12(1,-\\sqrt2,1)$.</li>
        <li>$\\hat Y$ is Hermitian – the mirror reflection $1\\leftrightarrow3$. Eigenvalue $+1$: $|2\\rangle$, $(|1\\rangle+|3\\rangle)/\\sqrt2$; $-1$: $(|1\\rangle-|3\\rangle)/\\sqrt2$. $[\\hat H,\\hat Y]=0$.</li>
        <li>$P(-\\sqrt2\\beta)=\\tfrac14$, $P(0)=\\tfrac12$, $P(\\sqrt2\\beta)=\\tfrac14$. (Dynamics: full transfer $1\\to3$ at $t=\\pi/(\\sqrt2\\beta)$.)</li></ol>`,
      preset: { builder: 'network', values: { topology: 'chain', N: 3, eps: 0, J: -1, gamma: 0 }, initial: { mode: 'site', k: 0 }, observables: ['H', 'Y', 'Q'], track: ['Q'], tMax: 2 * Math.PI / Math.SQRT2 * 2 },
    },
    {
      id: 'dots-chain3-spin', topic: 'dots', engine: 'finite', source: 'Uloha2_2016.pdf', sheet: 'Úloha 2 (2016)',
      title: 'Three-dot chain with spin',
      statement: `The three-dot chain ($\\hat H$, $\\hat Y$ as in the previous problem) now carries spin ½: $\\mathcal H_3\\otimes\\mathcal H_s$,
        $$\\hat H'=\\hat H\\otimes\\hat 1+\\hat V,\\qquad \\hat V=\\frac{2\\lambda}{\\hbar}(2\\hat Y-\\hat 1)\\otimes\\hat s_x,\\quad 0<\\lambda\\ll\\beta.$$
        <ol><li>Possible energy values of $\\hat H'$.</li><li>The system is in its ground state; $\\hat s_y$ is measured. Outcomes, probabilities, post-measurement state.</li></ol>`,
      solution: `<ol><li>$\\hat V=\\lambda(2\\hat Y-1)\\otimes\\hat\\sigma_x$ commutes with $\\hat H\\otimes1$; On the $\\hat H$ eigenvectors $\\hat Y=+1$ (for $E=\\mp\\sqrt2\\beta$) or $\\hat Y=-1$ (for $E=0$), so $2\\hat Y-1=+1$ or $-3$. Energies: $-\\sqrt2\\beta\\pm\\lambda$, $\\pm3\\lambda$, $\\sqrt2\\beta\\pm\\lambda$ – all non-degenerate (exact, not perturbative).</li>
        <li>Ground state $\\tfrac12(1,\\sqrt2,1)\\otimes|x-\\rangle$, $E=-\\sqrt2\\beta-\\lambda$. $s_y=\\pm\\hbar/2$ each with probability ½; afterwards $\\tfrac12(1,\\sqrt2,1)\\otimes|y\\pm\\rangle$.</li></ol>`,
      preset: { builder: 'network', values: { topology: 'chain', N: 3, J: -1, spin: true, lambda: 0.1 }, initial: { mode: 'eigen', k: 0 }, observables: ['sy', 'H', 'sx', 'sz', 'Y', 'Q'], track: ['sy'], tMax: 40, blochLabels: { z: ['↑ z', '↓'] } },
    },
    {
      id: 'dots-k4', topic: 'dots', engine: 'finite', source: 'Vzorove_reseni_2014.pdf', sheet: 'Test 2014 – Úloha 2',
      title: 'Four dots, all-to-all tunnelling',
      statement: `$\\hat H=\\varepsilon\\sum_m|m\\rangle\\langle m|-t\\sum_{m\\ne n}|m\\rangle\\langle n|$ on four dots (complete graph).
        Verify that $(1,1,1,1)$ is an eigenvector, construct the projector $\\hat P_0$ onto it and $\\hat P_1=1-\\hat P_0$, and write $\\hat H=\\alpha\\hat P_0+\\beta\\hat P_1$.
        <p class="note">Mapping: $J=-t$. The complete-graph topology also works for other $N$ (try other N).</p>`,
      solution: `$\\hat P_0=\\tfrac14\\sum_{mn}|m\\rangle\\langle n|$, $\\hat H=(\\varepsilon-3t)\\hat P_0+(\\varepsilon+t)\\hat P_1$: eigenvalue $\\varepsilon-3t$ (non-degenerate) and $\\varepsilon+t$ (three-fold), e.g. $(3,-1,-1,-1)/\\sqrt{12}$, $(0,2,-1,-1)/\\sqrt6$, $(0,0,1,-1)/\\sqrt2$.
        Starting from $|1\\rangle$: $P_1(t)=\\tfrac58+\\tfrac38\\cos(4t\\,\\tau/\\hbar)$ – a single frequency.`,
      preset: { builder: 'network', values: { topology: 'complete', N: 4, eps: 0, J: -1, gamma: 0 }, initial: { mode: 'site', k: 0 }, observables: ['H', 'Q'], tMax: 2 * Math.PI },
    },
    {
      id: 'chain-infinite', topic: 'dots', engine: 'finite', source: 'Uloha3_2016.pdf', sheet: 'Úloha 3 (2016)',
      title: 'Infinite quantum chain – velocity operator',
      statement: `$\\mathcal H=\\ell^2(\\mathbb Z)$, $\\hat H=\\hbar\\omega\\sum_n|n\\rangle\\langle n|-\\hbar\\beta\\sum_n(|n\\rangle\\langle n+1|+|n+1\\rangle\\langle n|)$, $\\hat Q=\\sum_n n|n\\rangle\\langle n|$.
        <ol><li>Velocity $\\hat V=\\tfrac i\\hbar[\\hat H,\\hat Q]$.</li><li>Show $\\tfrac{d}{dt}\\langle\\hat Q\\rangle=\\langle\\hat V\\rangle$.</li><li>Spectrum of $\\hat V$.</li><li>Normalised (generalised) eigenvectors (hint $\\psi_n=e^{ikn}$).</li><li>Continuous or discrete? Degeneracy? Completeness.</li><li>For $|\\psi\\rangle=a|0\\rangle+b|1\\rangle$ find $a,b$ maximising $\\langle\\hat V\\rangle$.</li></ol>
        <p class="note">Simulated as a long chain (101 sites, $\\varepsilon=\\omega$, $J=-\\beta$). Watch $\\langle Q\\rangle(t)$ grow linearly with slope $\\langle V\\rangle$; try the Gaussian packet with different $k_0$.</p>`,
      solution: `<ol><li>$\\hat V=i\\beta\\sum_n(|n+1\\rangle\\langle n|-|n\\rangle\\langle n+1|)$.</li><li>Heisenberg: $\\tfrac{d}{dt}\\langle Q\\rangle=\\tfrac i\\hbar\\langle[H,Q]\\rangle=\\langle V\\rangle$.</li>
        <li>On $|k\\rangle=(2\\pi)^{-1/2}\\sum_n e^{ikn}|n\\rangle$, $k\\in(-\\pi,\\pi]$: $\\hat V|k\\rangle=2\\beta\\sin k\\,|k\\rangle$ – spectrum $[-2\\beta,2\\beta]$; consistent with $E(k)=\\hbar\\omega-2\\hbar\\beta\\cos k$, $v=\\tfrac1\\hbar\\,dE/dk$.</li>
        <li>$\\langle k|k'\\rangle=\\delta(k-k')$.</li><li>Purely continuous, doubly degenerate ($k$ and $\\pi-k$) except at $\\pm2\\beta$; $\\int_{-\\pi}^{\\pi}|k\\rangle\\langle k|\\,dk=1$.</li>
        <li>$\\langle V\\rangle=2\\beta\\,\\mathrm{Im}(a^*b)$, maximal ($=\\beta$) for $a=1/\\sqrt2$, $b=i/\\sqrt2$.</li></ol>`,
      preset: { builder: 'network', values: { topology: 'chain', N: 101, eps: 0, J: -1, gamma: 0, labelStyle: 'centered' }, initial: { mode: 'twosite', theta: Math.PI / 2, phi: Math.PI / 2 }, observables: ['Q', 'V', 'H'], track: ['Q'], tMax: 20, showBand: true },
    },

    // ---------------------------------------------------------------- composite
    {
      id: 'entangle-qubits', topic: 'composite', engine: 'finite', source: '3-Uloha_3_entanglement_Qbitu.pdf', sheet: 'Úloha 3 – entanglement',
      title: 'Entangling two qubits',
      statement: `Qubit: $\\hat h=\\varepsilon|1\\rangle\\langle1|$. Two qubits: $\\hat H_0=\\hat h\\otimes1+1\\otimes\\hat h$, interaction $\\hat V=\\hbar\\omega(\\hat\\sigma_x\\otimes\\hat\\sigma_y-\\hat\\sigma_y\\otimes\\hat\\sigma_x)$ switched on for $0\\le t\\le T$. Initially $|\\psi\\rangle=|1\\rangle\\otimes|0\\rangle$.
        <ol><li>Evolution under $\\hat H_0$ alone.</li><li>State at $T$ under $\\hat H_0+\\hat V$.</li><li>Show the state at $T=\\pi/(8\\omega)$ is entangled.</li><li>Then the energy of qubit 1 and immediately of qubit 2 is measured: outcomes and probabilities.</li></ol>
        <p class="note">Mapping: $\\varepsilon_1=\\varepsilon_2=\\varepsilon$, $\\omega$ term, switch-off time $T$. The Bloch balls show each qubit's reduced state – they shrink to the centre at maximal entanglement.</p>`,
      solution: `<ol><li>$|10\\rangle$ is an eigenstate of $\\hat H_0$: $|\\psi(t)\\rangle=e^{-i\\varepsilon t/\\hbar}|10\\rangle$.</li>
        <li>$\\hat V$ annihilates $|00\\rangle,|11\\rangle$; $\\hat V|10\\rangle=2i\\hbar\\omega|01\\rangle$, $\\hat V|01\\rangle=-2i\\hbar\\omega|10\\rangle$. Hence $|\\psi(T)\\rangle=e^{-i\\varepsilon T/\\hbar}[\\cos(2\\omega T)|10\\rangle+\\sin(2\\omega T)|01\\rangle]$ (full swap at $T=\\pi/4\\omega$).</li>
        <li>At $T=\\pi/8\\omega$: $(|10\\rangle+|01\\rangle)/\\sqrt2$ (Bell state) – not a product; concurrence $C=|\\sin4\\omega T|=1$, reduced states maximally mixed.</li>
        <li>$(\\varepsilon,0)$ with probability ½, $(0,\\varepsilon)$ with probability ½; never $(0,0)$ or $(\\varepsilon,\\varepsilon)$. In general $\\cos^2 2\\omega T$ and $\\sin^22\\omega T$.</li></ol>`,
      preset: {
        builder: 'spins', builderChoices: ['spins', 'orbitalSpin'], values: { e: 1, D: 1, switchOn: true, T: Math.PI / 8 },
        initial: { mode: 'site', k: 2 }, observables: ['E1', 'E2', 'Z1', 'Z2', 'H', 'H0'], tMax: 2, blochLabels: { z: ['|0⟩', '|1⟩'] },
        statePresets: [{ label: '|10⟩', vector: [0, 0, 1, 0] }, { label: 'Bell (|01⟩+|10⟩)/√2', vector: [0, S, S, 0] }],
      },
    },
    {
      id: 'two-spins-A', topic: 'composite', engine: 'finite', source: 'Uloha2.pdf', sheet: 'Úloha 2',
      title: 'Two spins: measuring A = (ω/ħ)(Sx² − Sy²)',
      statement: `Two spin-½ particles, $\\hat{\\vec S}=\\hat{\\vec s}^{(1)}\\otimes1+1\\otimes\\hat{\\vec s}^{(2)}$, observable $\\hat A=\\frac\\omega\\hbar(\\hat S_x^2-\\hat S_y^2)$.
        <ol><li>Possible values of $A$.</li><li>Probabilities and post-measurement states for $|x+\\rangle\\otimes|x+\\rangle$.</li><li>The same for $|x+\\rangle|x-\\rangle-|x-\\rangle|x+\\rangle$.</li></ol>
        <p class="note">Basis $|\\pm\\pm\\rangle$ = $S_z$ eigenstates. Choose the input state in the amplitude editor (presets included) or as a product of two Bloch vectors.</p>`,
      solution: `$\\hat S_x^2-\\hat S_y^2=\\tfrac12(\\hat S_+^2+\\hat S_-^2)$ ⇒ $\\hat A=\\omega\\hbar(|{+}{+}\\rangle\\langle{-}{-}|+|{-}{-}\\rangle\\langle{+}{+}|)$.
        <ol><li>$+\\omega\\hbar$: $(|{+}{+}\\rangle+|{-}{-}\\rangle)/\\sqrt2$; $-\\omega\\hbar$: $(|{+}{+}\\rangle-|{-}{-}\\rangle)/\\sqrt2$; $0$ (twice): $|{+}{-}\\rangle,|{-}{+}\\rangle$.</li>
        <li>$|x+x+\\rangle=\\tfrac12(1,1,1,1)$: $P(\\omega\\hbar)=\\tfrac12$ → $(|{+}{+}\\rangle+|{-}{-}\\rangle)/\\sqrt2$; $P(0)=\\tfrac12$ → $(|{+}{-}\\rangle+|{-}{+}\\rangle)/\\sqrt2$; $P(-\\omega\\hbar)=0$.</li>
        <li>The state is the singlet (rotation invariant): $P(0)=1$, state unchanged.</li></ol>`,
      preset: {
        builder: 'spins', builderChoices: ['spins', 'orbitalSpin'], values: {}, labels: ['++', '+−', '−+', '−−'], subLabels: ['spin 1', 'spin 2'],
        initial: { mode: 'amp', vector: [0.5, 0.5, 0.5, 0.5] }, observables: ['A', 'S2', 'Sz', 'Sn', 'Sn1', 'Sn2'], tMax: 1,
        statePresets: [{ label: '|x+⟩|x+⟩', vector: [0.5, 0.5, 0.5, 0.5] }, { label: '|x+⟩|x−⟩ − |x−⟩|x+⟩ (normalised)', vector: [0, -S, S, 0] }, { label: '|++⟩', vector: [1, 0, 0, 0] }],
      },
    },
    {
      id: 'two-spins-e', topic: 'composite', engine: 'finite', source: 'Vzorove_reseni_2014.pdf', sheet: 'Test 2014 – Úloha 3',
      title: 'Two spins: total spin along (1,1,1)/√3',
      statement: `Two spins in $|{+}{+}\\rangle$ (z basis). The total spin projection onto $\\vec e=(1,1,1)/\\sqrt3$, $\\hat S_e=\\frac\\hbar2(\\vec e\\cdot\\hat{\\vec\\sigma}\\otimes1+1\\otimes\\vec e\\cdot\\hat{\\vec\\sigma})$, is measured. Probability of $S_e=0$ and the state after that outcome.
        <p class="note">The axis is set by $\\theta_n=\\arccos(1/\\sqrt3)$, $\\varphi_n=\\pi/4$. Sweep $\\theta_n$: $P(0)=\\tfrac12\\sin^2\\theta$.</p>`,
      solution: `Zero eigenspace: $\\varphi_1=(0,1,-1,0)/\\sqrt2$ (singlet) and $\\varphi_2=(i-1,1,1,i+1)/\\sqrt6$. $P(0)=|\\langle\\varphi_2|{+}{+}\\rangle|^2=\\tfrac13$; post-measurement state $\\varphi_2$ (up to phase). Generally $P(0)=\\tfrac12\\sin^2\\theta$ with $\\cos\\theta=\\vec e\\cdot\\hat z=1/\\sqrt3$.`,
      preset: {
        builder: 'spins', builderChoices: ['spins', 'orbitalSpin'], values: { thN: Math.acos(1 / Math.sqrt(3)), phN: Math.PI / 4 }, labels: ['++', '+−', '−+', '−−'], subLabels: ['spin 1', 'spin 2'],
        initial: { mode: 'site', k: 0 }, observables: ['Sn', 'S2', 'Sz', 'Sn1', 'Sn2'], tMax: 1,
      },
    },
    {
      id: 'nobel-cavity', topic: 'composite', engine: 'finite', source: 'U1.pdf', sheet: 'Úloha 1 – Nobel prize 2012',
      title: 'Spin entangled with a cavity photon',
      statement: `A spin-½ passes through a microwave cavity (two states $|0\\rangle,|1\\rangle$ = photon number). Afterwards
        $$|\\psi\\rangle=\\frac{\\alpha}{\\sqrt2}(|+\\rangle+|-\\rangle)\\otimes|0\\rangle+\\frac{\\beta}{\\sqrt2}(|+\\rangle-|-\\rangle)\\otimes|1\\rangle=\\alpha|x+\\rangle|0\\rangle+\\beta|x-\\rangle|1\\rangle.$$
        The spin projection along $\\vec n$ in the $xz$ plane at $\\theta=\\pi/3$ from $z$ is measured.
        <ol><li>Projectors onto the measurement outcomes.</li><li>Probabilities of the outcomes and the state of the whole system after the measurement.</li></ol>
        <p class="note">Parameters: $|\\alpha|=\\cos\\chi$, $|\\beta|=\\sin\\chi$, relative phase $\\arg\\beta$.</p>`,
      solution: `<ol><li>$\\hat P_\\pm=\\tfrac12(1\\pm\\vec n\\cdot\\hat{\\vec\\sigma})\\otimes1_{\\rm cav}$, $\\vec n=(\\tfrac{\\sqrt3}2,0,\\tfrac12)$.</li>
        <li>$\\vec n$ is at angle $\\pi/6$ from $x$: $|\\langle n+|x+\\rangle|^2=\\cos^2\\tfrac\\pi{12}=\\tfrac{2+\\sqrt3}4$, $|\\langle n+|x-\\rangle|^2=\\tfrac{2-\\sqrt3}4$.
        $p_\\pm=\\tfrac12\\pm\\tfrac{\\sqrt3}4(|\\alpha|^2-|\\beta|^2)$. After outcome $+$: $|n+\\rangle\\otimes(\\alpha\\langle n+|x+\\rangle|0\\rangle+\\beta\\langle n+|x-\\rangle|1\\rangle)/\\sqrt{p_+}$ – a product state; the cavity is left in a superposition conditioned on the result.</li></ol>`,
      preset: {
        builder: 'spins', builderChoices: ['spins', 'orbitalSpin'], values: { thN: Math.PI / 3, phN: 0, chi: Math.PI / 4, phb: 0 }, labels: ['+,0', '+,1', '−,0', '−,1'], subLabels: ['spin', 'cavity'],
        extraParams: [{ type: 'heading', label: 'Amplitudes $\\alpha$, $\\beta$' }, { key: 'chi', label: '$\\chi$: $|\\alpha|=\\cos\\chi$', min: 0, max: Math.PI / 2, step: 0.01, fmt: pi }, { key: 'phb', label: '$\\arg\\beta$', min: 0, max: 2 * Math.PI, step: 0.01, fmt: pi }],
        initial: { mode: 'fn' }, initialLabel: 'α|x+⟩|0⟩ + β|x−⟩|1⟩',
        initialFn: (v, QM) => { const a = Math.cos(v.chi) * Math.SQRT1_2, b = Math.sin(v.chi) * Math.SQRT1_2, c = Math.cos(v.phb), s = Math.sin(v.phb); return QM.V.from([[a, 0], [b * c, b * s], [a, 0], [-b * c, -b * s]]); },
        observables: ['Sn1', 'Z2', 'Sz', 'Sn'], tMax: 1, blochLabels: { z: ['|+⟩ / |0⟩', '|−⟩ / |1⟩'] },
      },
    },
    {
      id: 'reduced-2p', topic: 'composite', engine: 'finite', source: 'UlohaV.pdf', sheet: 'Úloha 5 – reduced density matrix',
      title: 'Reduced spin density matrix of a 2p state',
      statement: `A hydrogen electron is in a $2p$ state ($l=1$) with $J^2=\\tfrac34\\hbar^2$ and $J_z=\\tfrac12\\hbar$, i.e. $|j=\\tfrac12,m_j=\\tfrac12\\rangle$ in $(l=1)\\otimes(s=\\tfrac12)$.
        Find the reduced density matrix of the spin, decide whether it is pure or mixed, and find the probabilities of the spin projection along $\\vec n=(1,0,1)/\\sqrt2$.
        <p class="note">Spin axis: $\\theta_n=\\pi/4$, $\\varphi_n=0$. Try other $j, m_j$ in "coupled state" and switch on spin–orbit or a field for dynamics.</p>`,
      solution: `$|\\tfrac12,\\tfrac12\\rangle=\\sqrt{\\tfrac23}|m_l{=}1\\rangle|\\downarrow\\rangle-\\sqrt{\\tfrac13}|m_l{=}0\\rangle|\\uparrow\\rangle$ (Condon–Shortley).
        $\\hat\\rho_s=\\mathrm{Tr}_{\\rm orb}|\\psi\\rangle\\langle\\psi|=\\mathrm{diag}(\\tfrac13,\\tfrac23)$ in $(\\uparrow,\\downarrow)$ – mixed ($\\mathrm{Tr}\\rho^2=\\tfrac59$), Bloch vector $(0,0,-\\tfrac13)$.
        $P(\\pm)=\\tfrac12(1\\pm\\vec n\\cdot\\vec r)=\\tfrac12(1\\mp\\tfrac1{3\\sqrt2})$: $P(+)\\approx0.382$, $P(-)\\approx0.618$.`,
      preset: { builder: 'orbitalSpin', builderChoices: ['spins', 'orbitalSpin'], values: { l: 1, thN: Math.PI / 4, phN: 0 }, initial: { mode: 'cg', j: 0.5, mj: 0.5 }, observables: ['Sn', 'J2', 'Jz', 'Lz', 'Sz', 'H'], track: ['Sz'], tMax: 20 },
    },

    // ---------------------------------------------------------------- sandboxes
    {
      id: 'lab-dots', topic: 'dots', engine: 'finite', sandbox: true, title: 'Free play: Bloch oscillations in a tilted chain',
      statement: 'A Gaussian packet in a chain with a linear tilt $F$ oscillates instead of running away (Bloch oscillations, period $2\pi/F$). Change topology, number of dots, disorder $W$ (Anderson localisation) or attach a spin.',
      preset: { builder: 'network', values: { topology: 'chain', N: 41, eps: 0, J: -1, gamma: 0, F: 0.25 }, initial: { mode: 'packet', center: 20, width: 3, k0: 0 }, track: ['Q'], tMax: 60 },
    },
    {
      id: 'lab-qubits', topic: 'composite', engine: 'finite', sandbox: true, title: 'Free play: XY spin chain',
      statement: 'An excitation (one flipped spin) hops along an XY chain of 4 spins – excitation transfer, entanglement spreads. Change the number of spins (1–6), couplings and fields, or switch the system to orbital $l$ ⊗ spin ½.',
      preset: { builder: 'spins', builderChoices: ['spins', 'orbitalSpin'], values: { Ns: 4, e: 1, de: 0.1, Jxx: 0.5, Jyy: 0.5 }, initial: { mode: 'product', t1: Math.PI, p1: 0 }, tMax: 15, blochLabels: { z: ['|0⟩', '|1⟩'] } },
    },
  ]);
})();
