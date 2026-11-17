# Implemented scientific conventions

The release preserves the existing calculation used by the figure workflows.
These implementation details and manuscript differences must remain explicit.

## Waveform and Fisher parameters

IMRPhenomXHM modes: (2,2), (2,1), (3,3), (3,2), (4,4).
The Fisher order is **M, eta, Dl, tc, phi_c, x1, x2, alpha_bar_dCS_squared**.
Inclination is fixed at the value supplied to each calculation, including the
inclination scan; it is not marginalized. This is eight parameters, whereas the
supplied manuscript lists nine including inclination.

The coupling coordinate has dimension length^4 and fiducial value zero.
`sqrt(cov[7,7])**0.25 * 1e-3` is the reported length uncertainty in kilometres.
This nonlinear legacy forecast convention is not a posterior credible interval.
Fisher coupling derivatives retain amplitude and phase terms. The amplitude
derivative is nonzero at the zero fiducial coupling; this differs from the
manuscript's phase-only statement. The nonspinning constraint is undefined because
the leading derivative vanishes; it is retained as NaN rather than interpolated.

## Frequency and noise

`f_min=max(f_PN(T_obs), f_detector_low)`, `deltaF=f_min/50`.
Fisher ends at `min(f_ISCO, f_detector_high)`.
SNR uses pure GR through the full detector band, ending at `f_detector_high`.

Nominal limits: TianQin/LISA/Taiji 1e-4–1 Hz, BBO/DECIGO 1e-3–100 Hz, ground
detectors 4–4999 Hz. CE/aLIGO noise suppresses frequencies below 6 Hz; ET's table
spans 1–1e4 Hz. LISA omits Galactic confusion noise; Taiji includes it.
The three ground-detector tables contain ASD; interpolation is squared for PSD.

The inherited network calculation uses `1/S1+1/S2` on each independently
constructed detector grid and adds the integrals. Noise functions are not all
truncated to their nominal bands. Possible overlap/double contributions need
scientific assessment; the cleanup does not silently change that convention.

Mass scans retain absolute eta/spin finite-difference steps for joint detectors,
versus relative steps for individual detectors. Inversion is unregularized;
singular matrices are not replaced by pseudoinverses. Nonfinite/nonpositive
coupling variance produces a failed constraint. Conditioning and step convergence
require care near equal mass or zero spin.

## Source conventions and release scope

Controlled figure scans treat the input mass and luminosity distance as independent
parameters. Masses are converted directly from solar masses to SI without an
implicit cosmological redshift factor. Inclinations are radians; user distances
are Mpc and internal distances are meters.

Population catalogs, random realizations, population-table calculations and
statistical-uncertainty summaries are not part of this figure-only source release.
They must not be inferred from these fixed-source plots. Earlier scientific audit
notes and raw results are retained in the local external archive.
