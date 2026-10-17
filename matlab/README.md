# SPOT — MATLAB

SPOT is an optical imaging emulation framework for space domain awareness. It generates detector images and co-registered truth by combining observer–target geometry, finite-exposure integration, target photometry, stellar and diffuse backgrounds, optics and detector response. Ground-sky, space-sky, Earth-limb and Earth-background observations share the same image-formation pipeline.

## Requirements

MATLAB R2022a or newer is the intended baseline. Dependencies depend on the configuration:

| Feature | Dependency |
|---|---|
| IAU2006 frame conversion in `full` presets | Aerospace Toolbox (`dcmeci2ecef`) |
| Poisson shot-noise sampling | Statistics and Machine Learning Toolbox (`poissrnd`) |
| TLE trajectory inputs | `satelliteScenario` and SGP4/SDP4 support in the installed MATLAB/toolbox combination |

`preview` uses approximate frame conversion, synthetic stars and disabled Poisson sampling. It changes the acquisition model as well as resolution. Full examples require the bundled CSV catalogue. Input-loading failures are reported explicitly.

## Quick start

Set MATLAB's Current Folder to this `matlab` directory:

```matlab
SPOT_Setup;
run(fullfile('Examples','RUN_Example_GroundSky.m'));
```

Each example contains editable exposure, motion, pointing, display and saving controls. For an in-memory run without figures or exports:

```matlab
SPOT_Setup;
C = SPOT_Config('ground','full');
C.Output.SaveResults = false;
C.Output.Plot = false;
C.Plot.Save = false;
Out = RUN_SPOT(C);
imageDN = Out.Exposures{1}.Image_DN;
```

Use `SPOT_Config('ground','preview')` for exploration. `SPOT_Config('list')` lists presets. Keep one installation active on the MATLAB path; `SPOT_Setup` checks which configuration and runner take precedence.

## Folder structure

| Path | Purpose |
|---|---|
| `+SPOT/` | Numerical model, labels, exports, batches and checks |
| `Configuration/` | Shared defaults and example overrides |
| `Examples/` | Scenes, exposure sequences, datasets, CSV selection and glints |
| `Visualization/` | Camera, detector, RA/DEC, stage, dataset and photometry plots |
| `Inputs/` | Catalogue, sensor/target pools, ephemeris and TLE fixtures |
| `Outputs/` | Runtime results and figures |
| `ZPaperPlots.m` | Publication figures from an existing acquisition |
| `Diagnosis.m` | Peak diffuse/specular irradiance summary from `Out` |

`RUN_SPOT.m` is the public launcher. The `+SPOT` name is required by MATLAB's package system: add the `matlab` root to the path, not `+SPOT` itself.

## Examples

| Script in `Examples/` | Configuration | Exposure / cadence |
|---|---|---|
| `RUN_Example_GroundSky.m` | `ground`, four seeded targets | 1.5 / 1.5 s |
| `RUN_Example_SpaceSky.m` | `space_sky`, 115° off nadir | 1.5 / 1.5 s |
| `RUN_Example_EarthBackground.m` | `space_earth`, 30° off nadir | 1.5 / 1.5 s |
| `RUN_Example_EarthLimb.m` | `space_limb`, sunlit limb, 70° off nadir | 1 / 0.8 s |
| `RUN_Example_MultiExposure.m` | Three ground exposures | 1.5 / 1 s |
| `RUN_Example_Dataset.m` | Ten `current_ground` samples, 20 pool targets each | 1 / 0.8 s |
| `RUN_Example_CSVPools.m` | `current_ground`, selected sensor, five targets from row 10 | 1 / 1 s |

Cadence is measured between exposure starts. Intervals overlap when cadence is shorter than exposure duration. Full presets use 1024 × 1024 sampling. The Earth-limb preset uses a 2026-05-01 epoch, a 200,000-star cap and explicit gain of 1e6 electrons per model unit. The other three principal regimes use effective-band photon conversion and a catalogue cap of 648,000, subject to magnitude filtering.

`current_*` profiles use linear seed propagation, midpoint camera attitude, explicit radiometric gain and a 200,000-star cap. They are distinct configurations. Additional presets cover tracking, mixed apparent morphologies, Cartesian/Keplerian states, geodetic/ECEF states, sampled ephemerides and TLEs. Extend the finite ephemeris coverage before requesting longer acquisitions.

## Motion, pointing and input pools

Example controls use `E`:

```matlab
E.Exposure_s = 1;
E.Target.AngleToLOS_deg = [90 0.2 90 0.2];
E.Target.Direction_deg = [-35 0 35 -20];
C = SPOT_ApplyExample(C,E);
Out = RUN_SPOT(C);
```

`AngleToLOS_deg` is the angle between initial target velocity relative to the observer and the observer-to-target line of sight: 0° receding, 90° transverse, 180° approaching. Changing only the angle preserves speed magnitude. Additional controls include `RelativeSpeed_mps`, `Range_m` and `InitialAngles_deg`.

Empty overrides preserve configured values. Scalars apply to every selected target; vectors follow selected-target order. `InitialAngles_deg` accepts one two-angle row or one row per target: ground seeds use azimuth/elevation offsets; space seeds use off-boresight/around-boresight angles. Seed overrides require seeded targets; edit trajectory inputs directly for other trajectory types.

Every target uses exposure-integrated rendering. Compact/trail labels follow brightness-weighted image motion and pixel/PSF scales; target class does not select a renderer. Apparent morphology also depends on pointing, range and exposure.

CSV selection precedes explicit overrides. Selection modes are `ordered`, `random` and `ids`, without replacement. A selected row does not guarantee visibility or detectability. CSV `area_m2` specifies total six-face cube area; `Surface.Area_m2` specifies individual face areas.

## Display, saving and results

```matlab
E.ShowPlots = true;
E.SavePlots = false;
E.SaveResults = false;
E.PlotProducts = {'camera','detector','radec','stages','dataset'};
E.PlotPhotometry = false;
E.PlotExposures = 1;  % also 'all' or an index vector
```

Display, figure saving and numerical export are independent. When `SaveResults` is enabled, `C.Output.SaveMAT`, `SaveCSV` and `SavePNG` select formats. `C.Plot.Formats` and `C.Plot.DPI` control figures. Results are written beneath `Outputs/<run timestamp>/<sample ID>/`.

| Product | Contents |
|---|---|
| `acquisition.mat` | Raw DN, intermediate stages, target histories, labels, configuration and provenance |
| `metadata.json` | Configuration, annotations and export scales |
| Target/selection CSVs | Target metadata and selected input rows |
| Image PNG | Percentile-scaled display preview |
| Binary/instance/semantic PNGs | Target support and labels |
| Valid-pixel PNG | Detector support |
| `dataset_index.csv` | Batch IDs, split, saved directories and exposure counts |

Use MAT arrays for numerical analysis. PNG intensities are display-scaled. MAT retains uint32 instance IDs; instance PNGs use uint16. Arrays use row/column indexing; pixel centres and inclusive `[xmin ymin width height]` boxes are one-based. Semantic labels are 0 background, 1 compact and 2 trail. Masks represent noiseless target support after optics, not detections from noisy images.

To replot, load `Out` from `acquisition.mat`, then call `SPOT_Plots(Out,1,{'camera','detector','radec'})`. `ZPaperPlots` provides publication layouts and export controls at the top of its script; `Diagnosis` prints peak irradiances.

For the glint demonstration, define `C` and `E` using an example's settings section, then run `Examples/RUN_SPOT_AllTargetGlints.m` instead of that example's final apply-and-run block. It configures eligible target attitudes and executes the acquisition. Material settings are illustrative.

`SPOT.batch(configurations)` processes a cell array and returns its index without retaining every acquisition. `SPOT.run(C,false)` runs without export or visualization. Start with a fresh configuration when editing seeds; `Out.Configuration` contains resolved runtime states.

## Checks and model scope

Run `SPOT.check` for native regression checks and `SPOT.convergence(C)` to examine tighter sampling. Cross-language fixture comparison is documented in [Python verification](../python/Verification/README.md).

Earth visibility is spherical, catalogue directions are static, targets are unresolved and PSFs are Gaussian. Diffuse radiances and default effective-band/material parameters are illustrative. `band_photon` derives conversion from aperture, effective wavelength and quantum efficiency; `legacy_gain` uses an explicit coefficient. Gaia G magnitudes are not automatically equivalent to the instrument band. Rapid motion, long exposures and narrow glints require sampling-convergence assessment and application-specific calibration.

## License and citation

Software uses the [MIT License](LICENSE). See the [repository overview](../README.md) for citation information and [data credits](../DATA_CREDITS.md) for attribution and optional observation-data requirements.
