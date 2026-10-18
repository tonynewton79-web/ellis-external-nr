# Ellis–Newton Phase B: freely evolving Z4 constraints

**Why this experiment matters.** The successful independent Einstein Toolkit Phase-A test used `Z4c::set_Theta_zero = yes`. The upstream code explicitly states that this **reduces Z4c to the BSSN limit**. Hence Phase A establishes successful external dust-coupled evolution with fourth-order Hamiltonian/momentum convergence, but does not verify free Z4 propagation.

Phase B v6 changes this to `Z4c::set_Theta_zero = no`, sets `kappa1=0.02`, and requests full-time CarpetX norms for `Theta`, `ZtC`, `allC`, `HC`, and `MtC` at N=16,24,32, every 40 RK4 steps until t=0.18. The independent `scripts/phaseB_z4_audit.py` script tests required data columns, synchronized time grids, nonzero Theta, positive density, resolution trends, and injected nonconvergence counterexamples. Failure is explicitly a **refusal**, not a fabricated pass.

Phase A v5 is preserved manual-only. Phase B is still comoving **dust** coupled to evolving geometry and an Ellis scale: no angular photon distribution or measured dust-radiation exchange is yet included. There is no gravitational-memory theorem.
