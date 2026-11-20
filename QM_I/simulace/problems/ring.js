/* Problems for the 'ring' engine (particle on a ring / quantum rotor). Preset schema: see engines/ring.js. */
(function () {
  'use strict';
  const R = String.raw;
  QM.addProblems([
    {
      id: 'rotor-basics', topic: 'ring', title: 'Quantum rotor – observables, eigenbasis, ladder',
      source: 'reseni_kvantovy_rotor.pdf', sheet: 'Cvičení 5 – quantum rotor (solution)', engine: 'ring',
      statement: R`A particle moves on a circle of radius $d$; its configuration is the angle $\varphi$. Operators:
        $\hat\varphi$, $\hat x = d\cos\hat\varphi$, $\hat y = d\sin\hat\varphi$, $\hat L = -i\,\partial_\varphi$, $\hat H = \hat L^2/2I$ and parity $\hat P\psi(\varphi) = \psi(\varphi+\pi)$.
        <ol>
          <li>Find the spectra of $\hat\varphi$, $\hat x$, $\hat y$.</li>
          <li>Which of the operators are mutually compatible? Compute $[\hat L,\hat x]$, $[\hat L,\hat y]$ and the relation between $\hat P$ and $\hat x$.</li>
          <li>Find a complete set of commuting observables built from $\hat x,\hat y$.</li>
          <li>Find the common eigenbasis of $\hat L,\hat H,\hat P$ and the matrix elements of $\hat x$ in it.</li>
          <li>Show that $\hat A^\dagger = \hat x + i\hat y$ acts as a ladder operator for $\hat L$.</li>
        </ol>
        In the lab a Gaussian packet with a kick $n_0$ runs around the ring; watch it disperse, reassemble on the opposite side at $T/2$ and fully revive at $T = 4\pi I$. Try the $\hat A^\dagger$, $\hat A$ and $\hat P$ buttons.`,
      solution: R`<ol>
          <li>$\hat\varphi$: continuous spectrum $[0,2\pi)$ (eigenfunctions $\delta_{2\pi}(\varphi-\varphi_0)$). $\hat x$ and $\hat y$: continuous spectrum $[-d,d]$; every inner value of $\hat x$ is doubly degenerate ($\varphi$ and $2\pi-\varphi$ give the same $x$), similarly for $\hat y$ ($\varphi$ and $\pi-\varphi$).</li>
          <li>From $[\hat L, f(\hat\varphi)] = -i f'(\hat\varphi)$: $[\hat L,\hat x] = i\hat y$, $[\hat L,\hat y] = -i\hat x$. $\hat P = e^{i\pi\hat L}$ commutes with $\hat L$ and $\hat H$, while $\hat P\hat x = -\hat x\hat P$, $\hat P\hat y = -\hat y\hat P$. Compatible sets: $\{\hat\varphi,\hat x,\hat y\}$ (all functions of $\varphi$) and $\{\hat L,\hat H,\hat P\}$.</li>
          <li>$\hat x$ alone is degenerate, but the pair $(x,y)$ fixes $\varphi$ uniquely: $\{\hat x,\hat y\}$ is a complete set.</li>
          <li>$\psi_n(\varphi) = e^{in\varphi}/\sqrt{2\pi}$, $n\in\mathbb Z$ (periodicity). $\hat L\psi_n = n\psi_n$, $E_n = n^2/2I$ (doubly degenerate for $n\neq0$), $\hat P\psi_n = (-1)^n\psi_n$. $\langle m|\hat x|n\rangle = \tfrac d2(\delta_{m,n+1}+\delta_{m,n-1})$ – see the <i>Matrix elements</i> tab.</li>
          <li>$\hat A^\dagger = \hat x + i\hat y = d\,e^{i\hat\varphi}$, $[\hat L,\hat A^\dagger] = \hat A^\dagger$, so $\hat A^\dagger|n\rangle = d|n+1\rangle$ and $\hat A = d\,e^{-i\hat\varphi}$ lowers. $\hat A^\dagger\hat A = \hat A\hat A^\dagger = d^2$: no state is annihilated – the ladder is unbounded in both directions (unlike the harmonic oscillator).</li>
        </ol>
        Dynamics: $\psi(\varphi,t) = \sum_n a_n e^{-in^2t/2I}\psi_n(\varphi)$. At $t = 4\pi I$ all phases are $1$ (revival); at $t = 2\pi I$ they are $(-1)^n$, i.e. $\psi(\varphi+\pi)$.`,
      preset: { builder: 'gauss', phi0: Math.PI / 2, sigma: 0.5, n0: 3, nMax: 8, I: 0.5, view: 'dens', tab: 'matrix', matOp: 'x' },
    },
    {
      id: 'rotor-ladder', topic: 'ring', title: 'Rotor ladder operators, (x̂ − 3iŷ)|ψ₀⟩',
      source: 'Uloha3_2018.pdf', sheet: 'Úloha 3 (2018)', engine: 'ring',
      statement: R`Particle on a ring of radius $d$ with $\hat A^\dagger = \hat x + i\hat y$, $\hat A = \hat x - i\hat y$.
        <ol>
          <li>Compute $[\hat L,\hat A^\dagger]$ and $[\hat L,\hat A]$.</li>
          <li>If $\hat L|\psi\rangle = m\hbar|\psi\rangle$, show that $\hat A^\dagger|\psi\rangle$ and $\hat A|\psi\rangle$ are eigenvectors of $\hat L$ and find the eigenvalues.</li>
          <li>Find the norms of $\hat A^\dagger|\psi\rangle$ and $\hat A|\psi\rangle$.</li>
          <li>Write the matrix of $\hat y$ in the eigenbasis of $\hat L$.</li>
          <li>Let $\hat L|\psi_0\rangle = 0$ and $|\psi\rangle = (\hat x - 3i\hat y)|\psi_0\rangle$. With what probability does a measurement of $\hat L$ give a positive / negative value?</li>
        </ol>
        The lab starts from the builder $(\hat x + c\,i\hat y)|0\rangle$ with $c = -3$. Check the bar chart and the stats $P(\ell>0)$, $P(\ell<0)$, then measure $\hat L$ in the <i>Measurement</i> tab.`,
      solution: R`<ol>
          <li>$[\hat L,\hat A^\dagger] = \hbar\hat A^\dagger$, $[\hat L,\hat A] = -\hbar\hat A$ (from $[\hat L,\hat x] = i\hbar\hat y$, $[\hat L,\hat y] = -i\hbar\hat x$).</li>
          <li>$\hat L\hat A^\dagger|\psi\rangle = (\hat A^\dagger\hat L + \hbar\hat A^\dagger)|\psi\rangle = (m+1)\hbar\,\hat A^\dagger|\psi\rangle$; analogously $\hat A|\psi\rangle$ has eigenvalue $(m-1)\hbar$.</li>
          <li>$\|\hat A^\dagger\psi\|^2 = \langle\psi|\hat A\hat A^\dagger|\psi\rangle = \langle\psi|\hat x^2+\hat y^2 + i[\hat x,\hat y]|\psi\rangle = d^2$ since $[\hat x,\hat y]=0$; likewise $\|\hat A\psi\| = d$.</li>
          <li>$\langle m'|\hat y|m\rangle = \frac{d}{2i}\left(\delta_{m',m+1} - \delta_{m',m-1}\right)$.</li>
          <li>$\hat x - 3i\hat y = d\big[\tfrac{1+c}{2}e^{i\varphi} + \tfrac{1-c}{2}e^{-i\varphi}\big]_{c=-3} = d\,(-e^{i\varphi} + 2e^{-i\varphi})$, so $\psi \propto -e^{i\varphi} + 2e^{-i\varphi}$, i.e. $|\psi\rangle = (-|1\rangle + 2|{-1}\rangle)/\sqrt5$. Hence $P(\ell>0) = 1/5$, $P(\ell<0) = 4/5$, $P(\ell=0) = 0$.</li>
        </ol>`,
      preset: { builder: 'xiy', cRe: -3, cIm: 0, nMax: 4, I: 0.5, view: 'reim', tab: 'matrix', matOp: 'y' },
    },
    {
      id: 'dirichlet', topic: 'ring', title: 'Dirichlet kernel and the periodic δ',
      source: 'reseni_kvantovy_rotor.pdf', sheet: 'Appendix A – periodic δ', engine: 'ring',
      statement: R`The completeness of the rotor eigenbasis $\psi_n = e^{in\varphi}/\sqrt{2\pi}$ reads
        $\sum_n \psi_n(\varphi)\psi_n^*(\varphi') = \delta_{2\pi}(\varphi-\varphi')$ with the periodic delta
        $$\delta_{2\pi}(\varphi) = \sum_{k\in\mathbb Z}\delta(\varphi - 2\pi k) = \frac{1}{2\pi}\sum_{n\in\mathbb Z} e^{in\varphi}.$$
        Study the partial sums (Dirichlet kernel) $D_N(\varphi) = \frac1{2\pi}\sum_{|n|\le N} e^{in\varphi}$: their closed form, height, width and integral, and in what sense $D_N\to\delta_{2\pi}$.
        Use the slider / animation; <i>Load ψ ∝ D_N into the ring</i> puts the corresponding state (all $a_n$ equal) into the rotor and shows its time evolution.`,
      solution: R`Geometric series: $\sum_{n=-N}^{N} e^{in\varphi} = e^{-iN\varphi}\frac{e^{i(2N+1)\varphi}-1}{e^{i\varphi}-1} = \frac{\sin((N+\frac12)\varphi)}{\sin(\varphi/2)}$, so
        $D_N(\varphi) = \frac{\sin((N+\frac12)\varphi)}{2\pi\sin(\varphi/2)}$.
        <ul>
          <li>Peak $D_N(2\pi k) = (2N+1)/2\pi$, first zeros at $\pm 2\pi/(2N+1)$: the peaks grow and narrow as $N\to\infty$.</li>
          <li>$\int_{-\pi}^{\pi} D_N\,d\varphi = 1$ for every $N$ (only the $n=0$ term survives).</li>
          <li>For a smooth $2\pi$-periodic $f$: $\int_{-\pi}^{\pi} D_N f\,d\varphi = \sum_{|n|\le N}\hat f_n \to f(0)$ (Fourier series), i.e. $D_N\to\delta_{2\pi}$ in the sense of distributions – not pointwise (the side lobes keep oscillating).</li>
        </ul>`,
      preset: { showDirichlet: true, dirichletN: 10, builder: 'comb', K: 5, nMax: 8, I: 0.5, view: 'dens' },
    },
    {
      id: 'lab-ring', topic: 'ring', title: 'Ring sandbox', engine: 'ring', sandbox: true,
      statement: R`Free play with the particle on a ring. Choose the Hamiltonian: the free rotor $\hat L^2/2I$, a rotor in a potential
        $V(\varphi) = V_1\cos(\varphi-\varphi_1) + V_2\cos2\varphi$ (pendulum / hindered rotor, diagonalised numerically in the truncated $|n\rangle$ basis),
        or a ring threaded by a magnetic flux, $\hat H = (\hat L-\Phi)^2/2I$ (Aharonov–Bohm). Build any initial state – packets, cat states, a top-hat on an arc,
        eigenstates of $\hat H$ and their superpositions, random states, or custom coefficients – apply $\hat A^\dagger$, $\hat A$, $\hat P$, watch the
        quantum carpet over a chosen time window, read the level diagram with the weights $|\langle E|\psi\rangle|^2$, inspect matrix elements and measure $\hat L,\hat H,\hat P,\hat x,\hat y$.
        The preset is a quantum pendulum: $V = -6\cos(\varphi - 3\pi/2)$ has its minimum at the bottom of the ring, and a Gaussian packet released at
        $\varphi_0 = 3\pi/2 + 1$ swings back and forth (compare with the classical period $2\pi/\omega$ quoted under the level diagram) while it slowly dephases.`,
      preset: { ham: 'pot', V1: -6, phi1: 3 * Math.PI / 2, V2: 0, I: 0.5, nMax: 12, builder: 'gauss', phi0: 3 * Math.PI / 2 + 1, sigma: 0.45, n0: 0,
        tAuto: false, tWin: 12, view: 'dens', tab: 'measure', probePhi: 3 * Math.PI / 2 },
    },
  ]);
})();
