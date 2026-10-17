# Cross-language verification

`Export_MATLAB_Parity.m` runs MATLAB and exports deterministic fixtures. `compare_matlab.py` uses their configurations to execute Python and compare image, scene, optics, geometry and annotation fields.

## Export MATLAB fixtures

Set MATLAB's Current Folder to the repository root:

```matlab
root = pwd;
addpath(fullfile(root,'python','Verification'));
Export_MATLAB_Parity(fullfile(root,'matlab'), ...
    fullfile(root,'python','Verification','matlab_fixtures'));
```

The default frame model is `iau2006`, requiring its MATLAB toolbox. An explicit third argument `'gmst_precession'` selects the approximate alternative. Fixtures use a 128 × 128 grid with scintillation and stochastic detector sampling disabled.

## Compare in Python

From `python/`, using the configured environment:

```bash
python Verification/compare_matlab.py Verification/matlab_fixtures
```

Numerical defaults are `rtol=1e-5` and `atol=1e-9`; masks require exact equality. Flags can override numerical tolerances. Failed comparisons exit nonzero; missing fixtures raise an error. Results go to `Verification/matlab_parity_results.json`, excluded from source commits.

Actual fixture comparisons establish agreement only for the tested configurations and tolerances. Python tests and source inspection alone do not establish MATLAB/Python numerical parity. Random outputs are not expected to match sample by sample.
