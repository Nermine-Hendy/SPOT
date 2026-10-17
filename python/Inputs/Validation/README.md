# Optional observation inputs

These files are needed for `03_Validation_Native_SatSim.ipynb` and `04_Validation_Local_Comparator.ipynb`. Core simulations, GUI and tests use the other bundled inputs.

Place the files directly in this directory:

| Filename | Purpose | Required |
|---|---|---|
| `Capture_00938 03_05_42Z.fits` | Measured Frigate frame | Yes |
| `StarDB_500k.mat` | Validation star catalogue | Yes |
| `wcs.fits` | Astrometric WCS solution | Either WCS file |
| `new-image.fits` | Alternative solved image carrying WCS | Either WCS file |

The notebooks prefer `wcs.fits`. Paths are editable in setup cells. Missing files raise errors; synthetic data are not substituted.

The companion `SPOT_Validation_Inputs.zip`, when available to the project owner, contains these files under `Validation/`. Copy that folder's contents here. Its `CHECKSUMS_SHA256.txt` provides file hashes.

The two full-frame FITS files are approximately 117.6 MiB each, above GitHub's 100 MiB regular-Git limit. The local `.gitignore` excludes research data from source commits. If redistribution permission is confirmed, distribute them through a release asset or data repository and document its URL and terms.
