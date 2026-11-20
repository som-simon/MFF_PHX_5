/* Problems for the 'wave1d' engine (preset schema documented in engines/wave1d.js). */
QM.addProblems([
  {
    id: 'gauss-packet', topic: 'wave1d', title: 'Gaussian wave packets', source: 'Uloha3.pdf', sheet: 'Úloha 3 – Gaussian packets', engine: 'wave1d',
    statement: 'A free particle is in the state given in momentum representation by $\\tilde\\psi(p)=\\exp\\{-(p-p_0)^2/(4\\Delta^2)\\}$ (up to normalisation). ' +
      '(1) Find $\\langle p\\rangle$ and $\\Delta p$. (2) Find $\\psi(x)$ and $\\Delta x$; show the packet has minimum uncertainty. ' +
      '(3) What does multiplying $\\tilde\\psi(p)$ by $e^{-ipx_0}$ do? (4) How large must $p_0$ be (in units of $\\Delta$) so that the particle moves to the right with probability above 99 %, i.e. $P(p>0)>0.99$?',
    solution: '<ol>' +
      '<li>$|\\tilde\\psi(p)|^2\\propto e^{-(p-p_0)^2/(2\\Delta^2)}$ is a normal distribution: $\\langle p\\rangle=p_0$, $\\Delta p=\\Delta$.</li>' +
      '<li>$\\psi(x)=\\frac1{\\sqrt{2\\pi}}\\int\\tilde\\psi(p)e^{ipx}dp\\propto e^{ip_0x}e^{-\\Delta^2x^2}$, so $|\\psi|^2\\propto e^{-2\\Delta^2x^2}$, $\\Delta x=\\frac1{2\\Delta}$ and $\\Delta x\\,\\Delta p=\\frac12$ – the minimum allowed by the uncertainty relation.</li>' +
      '<li>$e^{-ipx_0}$ in $p$-representation is the translation operator $e^{-i\\hat px_0}$: $\\psi(x)\\to\\psi(x-x_0)$, hence $\\langle x\\rangle=x_0$ (the momentum distribution is unchanged).</li>' +
      '<li>$P(p>0)=\\int_0^\\infty|\\tilde\\psi|^2dp=\\tfrac12\\left[1+\\mathrm{erf}\\!\\left(\\frac{p_0}{\\sqrt2\\Delta}\\right)\\right]>0.99 \\Rightarrow p_0\\gtrsim 2.33\\,\\Delta$.</li></ol>' +
      'In the simulation: move $p_0$ with threshold $p^*=0$ until the tail reads 0.99. Play the evolution: the packet moves with velocity $p_0/m$ and spreads, $\\Delta x(t)=\\Delta x(0)\\sqrt{1+(2\\Delta^2t/m)^2}$, while $|\\tilde\\psi(p)|^2$ stays constant.',
    preset: { pot: 'free', state: 'gauss', p0: 2, Delta: 1, x0: 0, pStar: 0, xStar: 0, L: 60, xView: 14, pView: 7, tMax: 5 },
  },
  {
    id: 'delta-well-p', topic: 'wave1d', title: 'δ-well bound state in momentum space', source: 'Vzorove_reseni_2014.pdf', sheet: 'Test 2014, Úloha 5', engine: 'wave1d',
    statement: 'A particle of mass $m$ is bound in the attractive potential $V(x)=-\\lambda\\delta(x)$. Find the bound state, its momentum distribution, and the probability that a momentum measurement gives $p>p_0$, where $p_0=m\\lambda$.',
    solution: 'Away from $x=0$: $\\psi\\propto e^{-\\kappa|x|}$; the δ gives the jump $\\psi\'(0^+)-\\psi\'(0^-)=-2m\\lambda\\psi(0)$, so $\\kappa=m\\lambda$, ' +
      '$\\psi(x)=\\sqrt\\kappa\\,e^{-\\kappa|x|}$, $E=-\\frac{\\kappa^2}{2m}=-\\frac{m\\lambda^2}2$.<br>' +
      'Fourier transform: $\\tilde\\psi(p)=\\sqrt{\\frac{2}{\\pi}}\\frac{\\kappa^{3/2}}{\\kappa^2+p^2}$, so $\\rho(p)=\\frac{2}{\\pi p_0}\\left[1+(p/p_0)^2\\right]^{-2}$ with $p_0=\\kappa=m\\lambda$.<br>' +
      '$P(p>p_0)=\\frac2\\pi\\int_1^\\infty\\frac{du}{(1+u^2)^2}=\\frac2\\pi\\left[\\frac\\pi8-\\frac14\\right]=\\frac14-\\frac1{2\\pi}\\approx0.0908$.<br>' +
      'The heavy $p^{-4}$ tail comes from the cusp of $\\psi$ at $x=0$. The state is stationary: playing the evolution leaves both densities unchanged (up to grid error).',
    preset: { pot: 'delta', lambda: 1, m: 1, state: 'delta', pStar: 1, xStar: 1, L: 40, xView: 5, pView: 6, tMax: 5 },
  },
  {
    id: 'displaced-ho', topic: 'wave1d', title: 'Displaced oscillator ground state (coherent state)', source: 'Vzorove_reseni_2014.pdf', sheet: 'Test 2014, Úloha 1', engine: 'wave1d',
    statement: 'For the harmonic oscillator consider $\\hat U=e^{\\alpha(\\hat a-\\hat a^\\dagger)}$ with real $\\alpha>0$. Show that $\\hat U$ is unitary, find what it does to the position, and compute $\\langle x\\rangle$ in the state $\\hat U|0\\rangle$. How does this state evolve in time?',
    solution: '$\\hat A=\\alpha(\\hat a-\\hat a^\\dagger)$ is anti-Hermitian ($\\hat A^\\dagger=-\\hat A$), so $\\hat U^\\dagger=e^{-\\hat A}=\\hat U^{-1}$. ' +
      'With $\\hat a=(\\hat x/x_0+i\\hat p x_0)/\\sqrt2$, $x_0=1/\\sqrt{m\\omega}$: $\\hat a-\\hat a^\\dagger=i\\sqrt2\\,x_0\\hat p$, so $\\hat U=e^{i\\sqrt2\\alpha x_0\\hat p}$ is a translation by $-\\sqrt2x_0\\alpha$: $\\hat U^\\dagger\\hat x\\hat U=\\hat x-\\sqrt2x_0\\alpha$. ' +
      'Hence $\\langle x\\rangle=-\\sqrt2\\,x_0\\alpha$ in $\\hat U|0\\rangle$. Equivalently $\\hat U=D(-\\alpha)$ and $\\hat U|0\\rangle=|\\beta=-\\alpha\\rangle$ is a coherent state, $\\langle x\\rangle=\\sqrt2x_0\\,\\mathrm{Re}\\,\\beta$.<br>' +
      'Time evolution: $|\\beta\\rangle\\to e^{-i\\omega t/2}|\\beta e^{-i\\omega t}\\rangle$, so $\\langle x\\rangle(t)=-\\sqrt2x_0\\alpha\\cos\\omega t$ – the packet oscillates rigidly like a classical particle, $\\Delta x=x_0/\\sqrt2$ stays constant. ' +
      'Fock distribution: Poisson with $\\langle n\\rangle=|\\alpha|^2$, conserved in time.',
    preset: { pot: 'harmonic', omega: 1, m: 1, state: 'coherent', aRe: 1, aIm: 0, xStar: 0, pStar: 0, L: 40, xView: 5, pView: 5 },
  },
  {
    id: 'lab-wave1d', topic: 'wave1d', title: '1D wavefunction sandbox', engine: 'wave1d', sandbox: true,
    statement: 'Free play with a general 1D Hamiltonian $H=\\hat p^2/2m+V(x)$: choose the potential (free, harmonic, δ-well, finite well, infinite box, barrier, step, double well, linear/gravity with optional floor, periodic) and the initial state (Gaussian, cat state, HO eigenstate $|n\\rangle$, coherent state, box eigenstate, numerical eigenstates of $V$ and superpositions of two of them, exponential packet). ' +
      'The preset sends a Gaussian with mean energy $\\langle E\\rangle\\approx p_0^2/2m = V_b$ onto a rectangular barrier: watch the packet split into a transmitted and a reflected part, with $T(t)$, $R(t)$ approaching the plane-wave prediction $\\int_{p>0}|\\tilde\\psi(p)|^2T(p)\\,dp$. ' +
      'Other ideas: double well with state <i>superposition of two eigenstates</i> $n_1=0$, $n_2=1$ (tunnelling oscillation); box with a Gaussian (quantum carpet and revival at $t=4mw^2/\\pi$); linear potential with floor (quantum bouncer); cat state in the harmonic potential.',
    preset: { pot: 'barrier', Vb: 2, bw: 1, state: 'gauss', p0: 2, Delta: 0.3, x0: -12, L: 80, xView: 30, pView: 5, xStar: 0.5, pStar: 0, tMax: 12 },
  },
]);
