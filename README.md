# deflsim

Interactive simulator of how CRT deflection-yoke geometry determines the path of an
accelerated electron. Open `index.html` in a browser — no build, no server.

- `physics.js` — toroidal deflection-coil geometry (turns wound on a ferrite ring
  core, closed loops including radial crossovers), screened-toroid core model
  (bore field = exact Biot–Savart of the inner legs × (1+κ) image factor,
  κ = (μᵣ−1)/(μᵣ+1); outer legs shunted by the ferrite), relativistic Boris
  trajectory integration. Defaults are the WG6100 vertical-coil operating point.
  Loads under node too, for testing.
- `app.js` — controls, canvas rendering (3D geometry view, field heatmap + trajectory
  side view, deflection-vs-current sweep), tooltips, light/dark.
- `edef.c` — the original C attempt, kept for reference. Its main physics gaps, fixed
  here: open (non-closed) current paths with no end turns, the full accelerating
  voltage applied as a uniform E field along the whole tube instead of the beam
  entering the yoke at anode energy, and a 4×μ₀ fudge factor standing in for turns.
