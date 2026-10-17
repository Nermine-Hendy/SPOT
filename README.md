# SPOT

**Optical Imaging Emulation Framework for Space Domain Awareness**

SPOT generates synthetic optical observations of resident space objects and co-registered truth products. It connects observer–target geometry and finite-exposure integration with target photometry, stellar and diffuse backgrounds, optics and detector response. A common pipeline supports ground-sky, space-sky, Earth-limb and Earth-background observations.

## Implementations

| Implementation | Entry points | Documentation |
|---|---|---|
| MATLAB | `SPOT_Setup`, `RUN_SPOT`, scripts in `Examples/` | [MATLAB README](matlab/README.md) |
| Python | `ssa_emulator`, Jupyter examples, Tk desktop GUI | [Python README](python/README.md) |

Each implementation includes its own inputs and can be run independently. Bitwise equality between languages is not assumed.

## Capabilities

- Trajectories from seeded states, Cartesian states, Keplerian elements, sampled ephemerides and TLEs.
- Ground and space pointing, finite-exposure trajectories, visibility and illumination.
- Surface-reflection photometry, body attitude/spin and optional atmospheric effects.
- Stellar, sky, Earth-disk and Earth-limb backgrounds.
- Gaussian optics, detector noise, saturation and quantization.
- Raw images, intermediate diagnostics, target histories, bounding boxes and binary/instance/semantic masks.
- Individual scenes, exposure sequences, CSV selection and dataset generation.

## Start a simulation

In MATLAB, set the Current Folder to `matlab/`:

```matlab
SPOT_Setup;
run(fullfile('Examples','RUN_Example_GroundSky.m'));
```

For Python, create an environment and install `python/requirements.txt` as described in the [Python README](python/README.md). Start Jupyter from `python/` and open `01_Examples.ipynb`. The GUI is available through `02_GUI.ipynb` or `gui.py`.

Full presets use the bundled stellar catalogue. Preview configurations use smaller grids, synthetic stars and simplified settings. Results are generated under each implementation's `Outputs/` directory.

## Validation and reproducibility

MATLAB includes `SPOT.check`; Python includes a pytest suite. Deterministic cross-language comparisons are described in [Verification](python/Verification/README.md).

Optional Python notebooks compare SPOT with native SatSim optical/noise routines or a clearly labelled local Gaussian/Poisson comparator. Both use observation-conditioned detector maps and require separate observation files. Their scope is documented in the [Python README](python/README.md#observation-conditioned-validation).

Record configurations, seeds, dependencies and the repository commit when reporting results. Random realizations need not match across languages. Assess sampling convergence and calibration for the application.

## Model scope

Targets are unresolved, catalogue directions are static, Earth visibility is spherical and PSFs are Gaussian. Diffuse radiances and default effective-band/material parameters are illustrative. Quantitative sensor predictions require application-specific inputs and calibration.

## Citation

Related work: Nermine Hendy, Ke Wang, Kandeepan Sithamparanathan, and Akram Al-Hourani, *SPOT: Optical Imaging Emulation Framework for Space Domain Awareness*.

Software citation metadata is provided in [CITATION.cff](CITATION.cff). Cite the exact release or commit used. Bibliographic identifiers should be taken from the published manuscript when available.

## License and data

Software is distributed under the [MIT License](LICENSE). Catalogue attribution and optional observation-data information are in [DATA_CREDITS.md](DATA_CREDITS.md). Third-party data retain their own terms.

The source tree excludes generated outputs and optional observation FITS/MAT files. Maintainer instructions for publication and tagged releases are in [PUBLISHING.md](PUBLISHING.md).
