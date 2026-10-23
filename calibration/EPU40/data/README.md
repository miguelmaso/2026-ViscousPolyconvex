# Experimental testing of 3D-printed elastomeric polyurethane EPU40

The CSV data files in this folder have been obtained from [[1]](#1).
Figures 4b, 5, 6 and 10a of the paper are vector graphics, so the curves were read directly
from the plotted polylines in the PDF (no pixel digitisation) and resampled at equally
spaced stretch increments; see `assets/extraction` and the check `assets/check_overlay.png`.
All specimens were pre-treated with three loading-unloading cycles (Mullins effect removed),
so the data correspond to the stable material. Room temperature (20 °C).

## Quasi-static loading

- id 1: Fig. 4b, quasi-static loading at strain rate 0.00032 /s, up to stretch 1.5, stretch increment 0.02.
- id 2: Fig. 10a, equilibrium stresses at the end of each step (30 min hold) of the multi-step
  relaxation test, stretch 1.5 to 5 in increments of 0.5, plus the origin. These are the data the
  authors use to fit the equilibrium (Carroll) model, Fig. 10b.

Fields:
- id
- stretch
- stress (nominal, Pa)

## One-cycle loading-unloading

Figs. 5a, 5b, 6a, 6b: maximum stretch 2, 3, 4, 5 at strain rates 0.384, 0.128 and 0.032 /s
(ids 10-21, sorted by maximum stretch, then by decreasing rate), stretch increment 0.04,
which is the native sampling of the published curves.
- id
- stretch
- stress (nominal, Pa)
- rate (/s)

Notes:
- Loading ends exactly at the maximum stretch. In several tests the published curve shows a
  vertical stress drop at the reversal point; the CSV keeps the peak as the last loading point
  and the unloading branch follows the curve after the drop.
- Unloading stops at the last published point (stress close to, but not exactly, zero).
- Some curves start with a small stress offset at stretch 1 (up to 0.07 MPa), kept as published.

## References

<a id="1">[1]</a> M. Hossain, R. Navaratne and D. Perić, 3D printed elastomeric polyurethane: Viscoelastic experimental characterizations and constitutive modelling with nonlinear viscosity functions, International Journal of Non-Linear Mechanics, 2020, vol 126, 103546, [10.1016/j.ijnonlinmec.2020.103546](https://doi.org/10.1016/j.ijnonlinmec.2020.103546)
