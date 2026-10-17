# Data sources and attribution

## Stellar catalogue

Both implementations include byte-identical copies of `Inputs/Stars_648000.csv`. The model configuration identifies the catalogue as a Gaia DR3 subset using Gaia G magnitudes and a 2016.0 reference epoch. Credit ESA/Gaia and the Gaia Data Processing and Analysis Consortium (DPAC).

The subset's extraction query and retrieval record are not included. The model does not propagate catalogue proper motion. Gaia G magnitudes and the instrument's effective band require an appropriate conversion for calibrated applications.

Follow the official [Gaia DR3 credit instructions](https://gea.esac.esa.int/archive/documentation/GDR3/Miscellaneous/sec_credit_and_citation_instructions/) and [Gaia data credits](https://www.cosmos.esa.int/web/gaia-users/credits) when publishing work using these data. The software license does not replace data-source attribution requirements.

## Sensor, target and trajectory inputs

Sensor/target pools and ephemeris/TLE fixtures define illustrative configurations. TLE and ephemeris files have explicit epochs and finite validity/coverage; they are not a live orbit-data feed.

## Optional observation-conditioned validation data

The Python validation notebooks require a Frigate observation, solved WCS and MAT star catalogue. Their required placement is documented in [python/Inputs/Validation/README.md](python/Inputs/Validation/README.md). These files are excluded from source distribution.

Redistribution terms for the observation files and MAT catalogue are not identified in the project materials. Maintainers should document the source, attribution and redistribution permission before publishing those files as release assets or in a data repository. Core simulations and tests do not require them.
