/* Problems for the 'angular' engine (preset schema documented in engines/angular.js). */
(function () {
  'use strict';
  QM.addProblems([
    {
      id: 'spherical-well-L2', topic: 'angular', engine: 'angular', source: 'Vzorove_reseni_2014.pdf', sheet: 'Test 2014 – Úloha 4',
      title: 'Spherical well: measuring L²',
      statement: `A particle in an infinite spherical well of radius $a$ is in the state $\\psi(\\vec r)=(x+y)^2\\sin(\\pi r/a)$ (unnormalised).
        Which values of $\\hat L^2$ can be measured, with which probabilities? What is the state after the smallest value was measured?
        <p class="note">Write $(x+y)^2=r^2(n_x+n_y)^2$: only the angular part $(n_x+n_y)^2$ matters for $\\hat L^2$; the radial factor $r^2\\sin(\\pi r/a)$ is untouched. Edit the expression or change $c$, $k$ in $(x+c\\,y)^k$ to explore.</p>`,
      solution: `$(n_x+n_y)^2=\\sin^2\\theta(1+\\sin2\\varphi)$ contains only $l=0$ and $l=2$:
        $P(l=0)=\\tfrac59$ ($L^2=0$), $P(l=2)=\\tfrac49$ ($L^2=6\\hbar^2$). Within $l=2$: $m=0$ carries $\\tfrac19$ of the total, $m=\\pm2$ each $\\tfrac16$, $m=\\pm1$ zero.
        After measuring $L^2=0$ the angular part becomes constant: $\\psi\\propto r^2\\sin(\\pi r/a)$ (radial part unchanged).`,
      preset: { expr: '(x + c*y)^k', c: 1, k: 2, projL: 0, radial: '$R(r)=r^2\\sin(\\pi r/a)$ – unaffected by $\\hat L^2$, $\\hat L_z$.' },
    },
    {
      id: 'lab-angular', topic: 'angular', engine: 'angular', sandbox: true, title: 'Y_lm decomposition sandbox',
      statement: 'Type any polynomial in the direction cosines $x,y,z$ (e.g. <code>z^3</code>, <code>(x+i*y)^2</code>, <code>x*y*z</code>) and see its $Y_{lm}$ content, $\\langle L^2\\rangle$, $\\langle L_z\\rangle$ and the effect of an $L^2$ measurement.',
      preset: { expr: 'x*y*z + z^2', c: 1, k: 1, projL: 2 },
    },
  ]);
})();
