# Processed problem sources

Tracker of which files in `QM_I/Problems/` have been turned into simulations.
Keep in sync with `SOURCES` in `problems/registry.js` (shown in the app under *Sources & progress*).

**Adding new problems:** list `QM_I/Problems/`, diff against the table below, then for each new problem
add a preset to an existing engine (`problems/<engine>.js`). Only write a new engine when no existing
family fits. Similar systems must reuse the same engine and controls.

| File | Content | Status | Problem ids (engine) |
|---|---|---|---|
| 3-Uloha_3_entanglement_Qbitu.pdf | entangling two qubits | done | entangle-qubits (finite) |
| DU1_2020.pdf | 4 quantum dots on a ring | done | dots-ring4 (finite) |
| QM-II-2015-du1.pdf | spin precession in B field, dephasing | done | spin-precession (spin) |
| QMIcv_09.30_zadani.pdf | linear algebra, Pauli matrices | partial – Ú1–4, Ú7 are proofs | pauli-decomp, pauli-eigen (spin) |
| QMIcv_10.07_zadani.pdf | spin ½ along general axis | done | spin-axis (spin) |
| reseni_kvantovy_rotor.pdf | quantum rotor (solution) | done | rotor-basics, dirichlet (ring) |
| U1.pdf | spin ⊗ cavity (Nobel 2012) | done | nobel-cavity (finite) |
| Uloha1_2018.pdf | three-dot chain | done | dots-chain3 (finite) |
| Uloha2.pdf | two spins, A = Sx² − Sy² | done | two-spins-A (finite) |
| Uloha2_2016.pdf | three-dot chain with spin | done | dots-chain3-spin (finite) |
| Uloha3.pdf | Gaussian wave packets | done | gauss-packet (wave1d) |
| Uloha3_2016.pdf | infinite chain, velocity operator | done | chain-infinite (finite) |
| Uloha3_2018.pdf | rotor ladder operators | done | rotor-ladder (ring) |
| UlohaV.pdf | reduced density matrix, 2p j=½ | done | reduced-2p (finite) |
| Vzorove_reseni_2014.pdf | test 2014, 5 problems | done | displaced-ho, delta-well-p (wave1d); dots-k4, two-spins-e (finite); spherical-well-L2 (angular) |
| Dirichletovo jádro … .htm | old standalone sim | ported | dirichlet (ring) |
| Kvantové tečky – časový vývoj.htm | old standalone sim | replaced | dots-ring4 (finite) |
| Kvantový rotor — … .htm | old standalone sim | ported | rotor-basics (ring) |

## Structure
One page (**lab**) per topic; problems are presets selectable inside the lab (`#/lab/<topic>/<problemId>`).
Each lab is a general simulator – Hamiltonian, parameters, system size and initial state are all choosable;
the presets only load values. Free-play preset = the `sandbox: true` entry of the topic.

| Lab (topic) | Engine | Generalisation |
|---|---|---|
| Spin ½ & measurement (`spin`) | spin | any state/mixed state, axis, field, ensemble dephasing, Pauli tool |
| Spins, qubits & entanglement (`composite`) | finite: `spins` / `orbitalSpin` | 1–6 spins with fields + XX/YY/ZZ/DM couplings; orbital l ≤ 3 ⊗ spin |
| Quantum-dot networks (`dots`) | finite: `network` | N dots, chain/ring/complete/star/grid, γ, tilt, disorder, optional spin |
| Particle on a ring (`ring`) | ring | free rotor / potential / flux, many initial states, Dirichlet kernel |
| 1D wavefunctions (`wave1d`) | wave1d | many potentials and initial states, split-step evolution |
| Angular momentum (`angular`) | angular | any polynomial angular function → Y_lm |
