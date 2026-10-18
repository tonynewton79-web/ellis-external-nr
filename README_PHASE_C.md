# Phase C: homogeneous Bianchi-I null radiation control

A genuine Einstein + dynamical-matter test with two **individually
conserved, antipodal massless beams** and pressureless comoving dust.

The new thorn EllisBeamX supplies Ttt, Txx and Ttx to SpacetimeX/Z4c at
every Runge-Kutta stage and evolves rho_d, Eplus, Eminus, and R. The
initial flat metric with uniform extrinsic curvature obeys 8*pi*rho_total
=3 H0^2. In exact homogeneous Bianchi I, density invariants are
rho_d R^3, Eplus R^3 sqrt(gxx), and Eminus R^3 sqrt(gxx). The
phaseC_audit.py checks positivity, conservation, balanced opposing
beams, nonzero anisotropic geometric response, spatial uniformity,
volume geometry and Z4c constraint sanity across N=16/24/32.

This is a homogeneous two-stream reduction, not inhomogeneous transport,
scattering/absorption or full Einstein-Boltzmann; no gravitational-memory
theorem follows. The earlier successful Phase-B workflow is manual-only.
