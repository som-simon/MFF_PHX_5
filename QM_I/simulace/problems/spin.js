/* Problems for the 'spin' engine (single spin ½). Preset schema: see engines/spin.js. */
(function () {
  'use strict';
  const PI = Math.PI;
  const T = String.raw; // keep TeX backslashes literal (never write "$" directly followed by "{")

  QM.addProblems([
    {
      id: 'pauli-decomp', topic: 'spin', title: 'Pauli matrices as a basis of 2×2 matrices',
      source: 'QMIcv_09.30_zadani.pdf', sheet: 'Cvičení 1, Ú5', engine: 'spin',
      statement: T`<p>Show that $(A,B)=\mathrm{Tr}(A^\dagger B)$ is an inner product on the space of complex $2\times2$ matrices.
        Show that $\{I,\sigma_x,\sigma_y,\sigma_z\}$ is an orthogonal basis (each of norm² $2$), so every matrix can be written as
        $M=a_0I+\vec a\cdot\vec\sigma$ with $a_0=\tfrac12\mathrm{Tr}\,M$, $a_i=\tfrac12\mathrm{Tr}(\sigma_iM)$.</p>
        <p><i>Simulation:</i> type any complex matrix (or load one); the tool returns $a_\mu$, checks Hermiticity and, for a Hermitian $M$,
        shows the eigenvalues $a_0\pm|\vec a|$ and the direction $\hat a=\vec a/|\vec a|$ on the Bloch sphere. The table shows all products $(\sigma_\mu,\sigma_\nu)$.</p>`,
      solution: T`<ol>
        <li>Inner product: $\mathrm{Tr}(A^\dagger B)$ is linear in $B$; $\mathrm{Tr}(A^\dagger B)^*=\mathrm{Tr}(B^\dagger A)$;
        $\mathrm{Tr}(A^\dagger A)=\sum_{ij}|A_{ij}|^2\ge0$ with equality only for $A=0$.</li>
        <li>From $\sigma_i\sigma_j=\delta_{ij}I+i\varepsilon_{ijk}\sigma_k$ and $\mathrm{Tr}\,\sigma_i=0$: $(\sigma_\mu,\sigma_\nu)=2\delta_{\mu\nu}$ (with $\sigma_0=I$).
        Four orthogonal non-zero elements of a 4-dimensional space form a basis.</li>
        <li>Coefficients: $a_\mu=(\sigma_\mu,M)/(\sigma_\mu,\sigma_\mu)=\tfrac12\mathrm{Tr}(\sigma_\mu M)$. $M$ is Hermitian ⇔ all $a_\mu\in\mathbb R$.</li>
        <li>$(M-a_0I)^2=(\vec a\cdot\vec\sigma)^2=(\vec a\cdot\vec a)\,I$ ⇒ eigenvalues $a_0\pm\sqrt{\vec a\cdot\vec a}$; for Hermitian $M$: $a_0\pm|\vec a|$
        with eigenvectors $|\hat a\pm\rangle$.</li>
        <li>Example $M=\begin{pmatrix}2&1-i\\1+i&0\end{pmatrix}$: $a_0=a_1=a_2=a_3=1$, eigenvalues $1\pm\sqrt3$.</li></ol>`,
      preset: { tools: ['pauli'], showSections: ['pauli'], pauli: { matrix: [[[2, 0], [1, -1]], [[1, 1], [0, 0]]] } },
    },
    {
      id: 'pauli-eigen', topic: 'spin', title: 'Eigenvectors of the Pauli matrices',
      source: 'QMIcv_09.30_zadani.pdf', sheet: 'Cvičení 1, Ú6', engine: 'spin',
      statement: T`<p>Find the eigenvalues and eigenvectors of $\sigma_x,\sigma_y,\sigma_z$ and the projectors onto the eigenvectors.
        Check $P_++P_-=I$. Expand the eigenvectors $|y\pm\rangle$ of $\sigma_y$ in the eigenbasis $\{|x\pm\rangle\}$ of $\sigma_x$.</p>
        <p><i>Simulation:</i> the state is $|y+\rangle$ and the axis is $\vec n=x$. The eigenbasis panel shows $|x\pm\rangle$, the projectors
        and the amplitudes $\langle x\pm|y+\rangle$. Switch the axis (buttons <i>n = x, y, z</i>) and the state (quick-state buttons).</p>`,
      solution: T`<ol>
        <li>$\sigma_i^2=I$, $\mathrm{Tr}\,\sigma_i=0$ ⇒ eigenvalues $\pm1$.</li>
        <li>$|x\pm\rangle=\tfrac1{\sqrt2}(1,\pm1)^T$, $|y\pm\rangle=\tfrac1{\sqrt2}(1,\pm i)^T$, $|z+\rangle=(1,0)^T$, $|z-\rangle=(0,1)^T$.</li>
        <li>$P_\pm=|i\pm\rangle\langle i\pm|=\tfrac12(I\pm\sigma_i)$, hence $P_++P_-=I$, $P_+P_-=0$.</li>
        <li>$\langle x+|y+\rangle=\tfrac{1+i}2$, $\langle x-|y+\rangle=\tfrac{1-i}2$, so $|y+\rangle=\tfrac{1+i}2|x+\rangle+\tfrac{1-i}2|x-\rangle$;
        similarly $|y-\rangle=\tfrac{1-i}2|x+\rangle+\tfrac{1+i}2|x-\rangle$. All $|\cdot|^2=\tfrac12$: measuring $\sigma_x$ on $|y\pm\rangle$ gives $\pm1$ with equal probability.</li></ol>`,
      preset: { state: { theta: PI / 2, phi: PI / 2 }, axis: { theta: PI / 2, phi: 0 }, tools: ['eigen'] },
    },
    {
      id: 'spin-axis', topic: 'spin', title: 'Spin ½ along a general axis n',
      source: 'QMIcv_10.07_zadani.pdf', sheet: 'Cvičení 2, Ú1–3', engine: 'spin',
      statement: T`<p>For $\vec n=(\sin\theta\cos\varphi,\sin\theta\sin\varphi,\cos\theta)$ write $\hat S_{\vec n}=\vec n\cdot\hat{\vec S}$ ($\hbar=1$) as a matrix,
        find its eigenvalues and eigenvectors $|n\pm\rangle$. Compute the probabilities of $S_{\vec n}=\pm\tfrac12$ in the states $|z+\rangle$ and $|x+\rangle$.</p>
        <p><i>Simulation:</i> move the axis sliders $\theta_n,\varphi_n$; the eigenvalues stay $\pm\tfrac12$. Compare the probabilities with the formulas;
        use <i>Measure</i> / <i>1000 shots</i> to sample.</p>`,
      solution: T`<ol>
        <li>$\hat S_{\vec n}=\tfrac12\begin{pmatrix}\cos\theta&e^{-i\varphi}\sin\theta\\e^{i\varphi}\sin\theta&-\cos\theta\end{pmatrix}$, $\hat S_{\vec n}^2=\tfrac14I$
        ⇒ eigenvalues $\pm\tfrac12$ for every $\vec n$.</li>
        <li>$|n+\rangle=(\cos\tfrac\theta2,\;e^{i\varphi}\sin\tfrac\theta2)^T$, $|n-\rangle=(-e^{-i\varphi}\sin\tfrac\theta2,\;\cos\tfrac\theta2)^T$.</li>
        <li>State $|z+\rangle$: $P(+)=\cos^2\tfrac\theta2$, $P(-)=\sin^2\tfrac\theta2$.</li>
        <li>State $|x+\rangle$: $\langle n+|x+\rangle=\tfrac1{\sqrt2}(\cos\tfrac\theta2+e^{-i\varphi}\sin\tfrac\theta2)$ ⇒ $P(\pm)=\tfrac12(1\pm\sin\theta\cos\varphi)$.</li>
        <li>In general, for a state with Bloch vector $\vec m$: $P(\pm)=\tfrac12(1\pm\vec n\cdot\vec m)$.</li></ol>`,
      preset: { state: { theta: 0, phi: 0 }, axis: { theta: PI / 3, phi: PI / 4 }, tools: ['eigen'] },
    },
    {
      id: 'spin-precession', topic: 'spin', title: 'Spin precession in a magnetic field & dephasing',
      source: 'QM-II-2015-du1.pdf', sheet: 'QM II 2015, DÚ 1', engine: 'spin',
      statement: T`<p>Spin-½ particles leave a Stern–Gerlach filter polarised along $\vec n=(\sin\theta,0,\cos\theta)$ (state $|n+\rangle$).
        They fly along $y$ through a region with a uniform magnetic field along $x$; the time of flight is $T=L/v$ and in the field
        $\hat H=-\tfrac{\hbar\omega}2\sigma_x$. Afterwards $S_z$ is measured.</p>
        <ol><li>Find $|\psi(T)\rangle$ and its Bloch vector.</li>
        <li>Find $P(S_z=\pm\tfrac\hbar2)$.</li>
        <li>Write the density matrix $\rho(T)$.</li>
        <li>The velocities are spread so that $T$ is uniformly distributed on $\langle T,T+\Delta T\rangle$ with $\omega\Delta T\gg1$. Find $\rho$; is it pure?</li>
        <li>What is $\rho$ after a non-selective measurement of $S_z$?</li></ol>
        <p><i>Simulation:</i> $\theta=\pi/3$, field along $x$ with the sign $\hat H=-(\omega/2)\sigma_x$, i.e. $\vec\Omega=-\omega\hat x$. The time slider is the flight time $T$.
        Tick <i>Average over flight times</i> and increase $\Delta T$ to see the ensemble dephase.</p>`,
      solution: T`<ol>
        <li>$|\psi(T)\rangle=e^{i\omega T\sigma_x/2}|n+\rangle=[\cos\tfrac{\omega T}2+i\sin\tfrac{\omega T}2\,\sigma_x]|n+\rangle$.
        The Bloch vector obeys $\dot{\vec P}=\vec\Omega\times\vec P$ with $\vec\Omega=-\omega\hat x$ (rotation about $-x$):
        $\vec P(T)=(\sin\theta,\ \cos\theta\sin\omega T,\ \cos\theta\cos\omega T)$.</li>
        <li>$P(S_z=\pm\tfrac12)=\tfrac12(1\pm\cos\theta\cos\omega T)$.</li>
        <li>$\rho(T)=\tfrac12(I+\vec P\cdot\vec\sigma)=\tfrac12\begin{pmatrix}1+\cos\theta\cos\omega T&\sin\theta-i\cos\theta\sin\omega T\\\sin\theta+i\cos\theta\sin\omega T&1-\cos\theta\cos\omega T\end{pmatrix}$.</li>
        <li>Averaging, $\overline{\cos\omega t},\overline{\sin\omega t}=O(1/\omega\Delta T)\to0$, so $\rho=\tfrac12(I+\sin\theta\,\sigma_x)$:
        $\mathrm{Tr}\rho^2=\tfrac12(1+\sin^2\theta)<1$ – mixed unless $\theta=\pi/2$ (only the component along the field survives).</li>
        <li>Non-selective $S_z$ measurement: $\rho\to P_+\rho P_++P_-\rho P_-$ removes the off-diagonal elements:
        $\rho=\tfrac12(I+\cos\theta\cos\omega T\,\sigma_z)$ for a sharp $T$, and $\rho=\tfrac12I$ for the ensemble of 4.
        (Button <i>Non-selective $S_z$</i>.)</li></ol>`,
      preset: {
        state: { theta: PI / 3, phi: 0 }, axis: { theta: 0, phi: 0 },
        field: { theta: PI / 2, phi: 0, omega: 1, sign: -1, fixed: true },
        dynamics: true, tMax: 4 * PI, ensemble: { on: false, dT: 30 },
      },
    },
    {
      id: 'lab-spin', topic: 'spin', title: 'Spin-½ sandbox', sandbox: true, engine: 'spin',
      statement: T`<p>Free single-qubit lab: prepare a pure or mixed state $\rho=\tfrac12(I+\vec P\cdot\vec\sigma)$, measure $S_x,S_y,S_z$ or $S_{\vec n}$,
        let it precess under $\hat H=\tfrac12\vec\Omega\cdot\vec\sigma$, average over an ensemble of evolution times, and decompose any $2\times2$ matrix into Pauli matrices.</p>`,
      preset: {
        state: { theta: PI / 3, phi: PI / 4 }, axis: { theta: PI / 4, phi: 0 },
        field: { theta: PI / 2, phi: 0, omega: 1 }, dynamics: true, ensemble: { on: false, dT: 20 }, tools: ['pauli', 'eigen'],
      },
    },
  ]);
})();
