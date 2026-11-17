# Retained detector amplitude spectral density tables

These are the existing numerical assets actually read by `Sn/Sn.py`:

| File | Used as | Columns |
| --- | --- | --- |
| `ET_D.txt` | ET-D interpolation | frequency [Hz], ASD [Hz^-1/2] |
| `CE_D40.txt` | CE 40 km interpolation | frequency [Hz], ASD [Hz^-1/2] |
| `aLIGO_D.txt` | aLIGO interpolation | frequency [Hz], ASD [Hz^-1/2] |

Noise functions return the **square** of the interpolated ASD. Out-of-band
sentinel noise and the nominal integration limits are described in the scientific
conventions document. `CE_D.txt` and `CE_D20.txt` were unused and archived.

The working tree did not contain trustworthy source URLs, release identifiers or
redistribution-license records for these tables. No such attribution is guessed.
Before publishing/relicensing the assets, the authors must recover their original
upstream versions, citations and permissions. Existing repository possession is
not evidence of a particular redistribution grant. SHA-256 hashes are stored in
`Sn/data/SHA256SUMS`; they identify the local numerical inputs, not their license.
