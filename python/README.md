# SPOT — Python

SPOT generates synthetic optical images for space domain awareness, together with geometry, radiometric and detector-stage diagnostics, target histories and annotation masks. Ground-sky, space-sky, Earth-limb and Earth-background regimes use a common finite-exposure image-formation model. This implementation provides Python functions, Jupyter examples and a Tk desktop interface.

## Requirements and installation

Use Python 3.11 or newer. Dependencies are in `requirements.txt`; Tk support and a desktop display are required for the GUI. The optional native SatSim comparison requires SatSim and compatible TensorFlow dependencies in the selected environment.

Open a terminal in this `python` directory. On Windows:

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python -m jupyterlab
```

On macOS/Linux:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m jupyterlab
```

Select that environment's kernel. Keep `ssa_emulator/`, `Inputs/`, notebooks and `gui.py` together. The package is imported from this directory; no editable installation is required. Linux installations may require Tk through the operating system's package manager.

## Quick start

Open `01_Examples.ipynb`, run setup, then one scene. The notebook defaults to `QUALITY='full'`. Running all cells also executes a ten-sample dataset batch.

For a compact programmatic run:

```python
from ssa_emulator import default_config, run

c = default_config('ground', 'preview')
out = run(c)
image_dn = out.Exposures[0].Image_DN
print(image_dn.shape)
```

`full` uses catalogue stars and 1024 × 1024 sampling. `preview` uses a 256 × 256 grid, synthetic stars, approximate frame conversion and disabled Poisson noise. It is a different acquisition configuration.

## Entry points

| File | Purpose | Additional requirement |
|---|---|---|
| `01_Examples.ipynb` | Seven workflows, extra trajectory presets and convergence analysis | Bundled inputs |
| `02_GUI.ipynb` | Launch the desktop GUI from Jupyter | Tk and a display |
| `gui.py` | Launch the same GUI from a terminal | Tk and a display |
| `03_Validation_Native_SatSim.ipynb` | Native SatSim PSF, FFT and photon-noise comparison | Observation files, SatSim/TensorFlow |
| `04_Validation_Local_Comparator.ipynb` | Local Gaussian/Poisson comparison or import of external SatSim results | Observation files; external arrays/provenance in external mode |

Launch the GUI with `.venv\Scripts\python gui.py` on Windows or `.venv/bin/python gui.py` on macOS/Linux. Select preset and quality, choose **Load example**, edit controls, then **Run simulation**. **Full configuration** exposes model parameters; **Edit target** exposes trajectory, material and attitude. **Export results** saves the completed acquisition. The GUI dataset preset runs one sample; use the notebook for batches.

## Examples

| Case | Purpose | Exposure / cadence |
|---|---|---|
| `ground` | Ground sky, four seeded targets | 1.5 / 1.5 s |
| `space_sky` | Space observer, 115° off nadir | 1.5 / 1.5 s |
| `space_earth` | Earth background, 30° off nadir | 1.5 / 1.5 s |
| `space_limb` | Sunlit limb, 70° off nadir | 1 / 0.8 s |
| `multi_exposure` | Three ground exposures | 1.5 / 1 s |
| `dataset` | Ten `current_ground` samples, 20 pool targets each | 1 / 0.8 s |
| `csv_pools` | `current_ground`, selected sensor, five targets from row 10 | 1 / 1 s |

Cadence is the interval between exposure starts. Exposures overlap when cadence is shorter than duration. The Earth-limb preset uses a 2026-05-01 epoch, a 200,000-star cap and explicit gain of 1e6 electrons per model unit. The other three principal regimes use effective-band photon conversion and a cap of 648,000, subject to magnitude filtering.

`current_*` profiles use linear seed propagation, midpoint camera attitude, explicit radiometric gain and a 200,000-star cap. They are distinct configurations. `default_config('list')` lists additional mixed-morphology, tracking and trajectory presets. Sampled ephemerides have limited coverage; the optional notebook demonstration uses 0.2 s. Extend coverage for longer acquisitions.

## Configure a scene

```python
from ssa_emulator.examples import example_config, run_example

c = example_config('ground', 'full')
c.Time.Exposure_s = 1.0
c.ExampleTarget.AngleToLOS_deg = [90, 0.2, 90, 0.2]
c.ExampleTarget.Direction_deg = [-35, 0, 35, -20]
c.Output.SaveResults = False
c.Output.Plot = True
c.Plot.Save = False
out = run_example(c)
```

Overrides include `AngleToLOS_deg`, `RelativeSpeed_mps`, `Direction_deg`, `InitialAngles_deg` and `Range_m`. The angle to LOS refers to initial velocity relative to the observer: 0° receding, 90° transverse, 180° approaching. Changing only the angle preserves speed magnitude. Empty values preserve settings; scalars broadcast and vectors follow selected-target order.

Overrides apply after CSV selection and before seed-state resolution. Seed overrides require seeded targets; edit trajectory inputs directly for other types. CSV modes are `ordered`, `random` and `ids`, without replacement. Build a fresh configuration when editing seeds; `out.Configuration` contains resolved runtime states.

All targets use exposure-integrated rendering. Compact/trail labels depend on brightness-weighted motion and pixel/PSF scales. Target class does not select a renderer. Python and MATLAB use the same model in scope, but random streams, frame libraries and numerical backends do not guarantee identical pixels. Cross-language agreement requires fixture comparisons.

## Saving and outputs

`run(c)` computes in memory. `run(c, export=True)` honours `c.Output.SaveResults`. `run_example(c)` also handles display and figure saving. The notebook exposes independent `SAVE_RESULTS`, `SAVE_PLOTS` and `SHOW_PLOTS` controls; saving defaults to false.

| Product | Contents |
|---|---|
| `exposure_####.npz` | Numerical DN, rate maps, stages, target histories and masks |
| Target/selection CSVs | Target properties, annotations and selected input rows |
| `metadata.json` | Configuration and provenance |
| PNG quicklooks | Display-scaled detector image and masks |
| Saved figures | Camera truth, detector, RA/DEC, 24-stage diagnostics, dataset products and photometry |

Directories under `Outputs/` use unique timestamps. An explicitly named existing export directory is rejected. `batch(...)` returns a sample index without retaining every image; `dataset_configurations(...)` constructs batches.

Use NPZ for numerical analysis. PNG detector intensities are display-scaled. NPZ preserves uint32 instance IDs; instance PNGs are exported only when IDs fit uint16. Arrays use zero-based `[row, column]` indexing. Geometry uses one-based `(x,y)` pixel centres and inclusive `[xmin,ymin,width,height]` boxes. Semantic labels are 0 background, 1 compact and 2 trail. Masks represent noiseless target support after optics.

## Observation-conditioned validation

Both validation notebooks require the files described in [Inputs/Validation/README.md](Inputs/Validation/README.md). Core examples, GUI and tests do not need these files. Simulator examples use the CSV catalogue; observation comparison uses a separate MAT catalogue and solved WCS.

The measured 6422 × 9600 frame is reduced to 1605 × 2400 using 4 × 4 block means, discarding incomplete edge blocks. The working ROI is 250 × 205 pixels. Positions, six trail endpoint pairs, empirical background, amplitude scales, PSF and gain are conditioned on the observation. Results assess detector-plane consistency and agreement with that observation, not independent validation of orbit geometry, target BRDF or absolute radiometry.

`03_Validation_Native_SatSim.ipynb` imports native SatSim PSF, FFT and photon-noise functions. It does not execute the `gen_multi` catalogue/orbit scene generator. Install SatSim using its [official instructions](https://github.com/ssc-ai/satsim) in a compatible environment; base requirements do not install it. The notebook records package/source information and checks deterministic consistency before calculating metrics.

`04_Validation_Local_Comparator.ipynb` defaults to a Gaussian/Poisson comparator labelled **SatSim-like**, which does not execute SatSim. Its `external_satsim` mode reads registered full-working-frame FITS/NPY/NPZ arrays and per-image JSON provenance. Required fields are `simulator`, `version`, `noise_on`, `downsample`, `units` and `scene_config`. It applies no automatic registration, resizing or intensity fitting to those arrays.

WCS inversion may warn for some catalogue entries. Inspect warnings and projected positions before interpreting counts or correspondence. Each run saves arrays, hashes and configurations.

## Checks and model scope

From this directory, using the configured environment:

```bash
python -m pytest -q
```

Tests cover geometry, photometry, detector behaviour, exports, overrides, Earth-limb behaviour and GUI controller/worker paths. They do not require optional observation files or a native window. See [Verification/README.md](Verification/README.md) for cross-language checks.

Earth visibility is spherical, catalogue directions are static, targets are unresolved and PSFs are Gaussian. Diffuse radiances and effective-band/material defaults are illustrative. Gaia G magnitudes require appropriate instrument-band conversion. Long exposures, rapid motion and narrow glints require convergence assessment; `convergence(c)` evaluates tighter sampling settings.

## License and citation

Software uses the [MIT License](LICENSE). See the [repository overview](../README.md) for citation information and [data credits](../DATA_CREDITS.md) for attribution and optional observation-data requirements.
