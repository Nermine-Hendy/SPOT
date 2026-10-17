"""Public examples and overrides, corresponding to the MATLAB example launchers."""
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from .config import default_config, struct

EXAMPLES = ('ground', 'space_earth', 'space_limb', 'space_sky',
            'multi_exposure', 'dataset', 'csv_pools')
TARGET_FIELDS = ('AngleToLOS_deg', 'RelativeSpeed_mps', 'Direction_deg',
                 'InitialAngles_deg', 'Range_m')
POINTING_FIELDS = ('AzEl_deg', 'OffNadir_deg', 'SpaceAzimuth_deg', 'Roll_deg')


def example_config(name='ground', quality='full'):
    """Create a fresh configuration. No simulation, output or state resolution."""
    if name not in EXAMPLES:
        return default_config(name, quality)
    c = default_config('current_ground' if name in ('dataset', 'csv_pools') else name, quality)
    c.ExampleTarget = struct({key: [] for key in TARGET_FIELDS})
    c.ExamplePointing = struct({key: [] for key in POINTING_FIELDS})
    c.Plot.Exposures = 'all' if name == 'multi_exposure' else 1
    if name in ('dataset', 'csv_pools'):
        c.Input.UseTargetCSV = True
        if name == 'dataset':
            c.Input.NumTargets = 20
            c.Input.TargetSelection.Mode = 'random'
            c.Output.Plot = c.Plot.Save = False
        else:
            c.Time.Cadence_s = 1
            c.Input.UseSensorCSV = True
            c.Input.NumTargets = 5
            c.Input.TargetSelection.StartRow = 10
            c.Input.TargetSelection.Seed = 100
            c.Output.SaveResults = False
    return c


def apply_example_overrides(c):
    """Apply AFTER CSV resolution, BEFORE converting seed geometry into states."""
    pointing = c.get('ExamplePointing', {})
    if set(pointing) - set(POINTING_FIELDS):
        raise ValueError('Unsupported example pointing override')
    for name, value in pointing.items():
        if value is not None and np.size(value):
            c.Pointing[name] = deepcopy(value)
    target = c.get('ExampleTarget', {})
    if set(target) - set(TARGET_FIELDS):
        raise ValueError('Unsupported example target override')
    active = {key: value for key, value in target.items()
              if value is not None and np.size(value)}
    if not active:
        return c
    n = len(c.Targets)
    if not n:
        raise ValueError('Target overrides require at least one target')
    values = {}
    for name, value in active.items():
        a = np.asarray(value, dtype=float)
        if not np.isfinite(a).all():
            raise ValueError(f'{name} must be finite')
        if name == 'InitialAngles_deg':
            if a.shape == (2,):
                a = a[None, :]
            if a.shape not in ((1, 2), (n, 2)):
                raise ValueError('InitialAngles_deg needs one [angle1, angle2] row or N rows')
            a = np.broadcast_to(a, (n, 2))
        else:
            if a.size not in (1, n) or a.ndim > 1:
                raise ValueError(f'{name} must be scalar or have one value per selected target')
            a = np.broadcast_to(a, (n,))
            if name == 'AngleToLOS_deg' and np.any((a < 0) | (a > 180)):
                raise ValueError('AngleToLOS_deg must be in [0, 180]')
            if name == 'RelativeSpeed_mps' and np.any(a < 0):
                raise ValueError('RelativeSpeed_mps must be nonnegative')
            if name == 'Range_m' and np.any(a <= 0):
                raise ValueError('Range_m must be positive')
        values[name] = a
    for j, target in enumerate(c.Targets):
        if target.Trajectory.Type != 'seeded_state':
            raise ValueError('Motion/angle overrides require fresh seeded targets; leave them empty for orbital, ephemeris or TLE inputs')
        s = target.Seed
        if 'AngleToLOS_deg' in values or 'RelativeSpeed_mps' in values:
            speed = np.hypot(s.TransverseSpeed_mps, s.RadialSpeed_mps)
            angle = np.rad2deg(np.arctan2(s.TransverseSpeed_mps, s.RadialSpeed_mps))
            speed = values['RelativeSpeed_mps'][j] if 'RelativeSpeed_mps' in values else speed
            angle = values['AngleToLOS_deg'][j] if 'AngleToLOS_deg' in values else angle
            s.TransverseSpeed_mps = float(speed * np.sin(np.deg2rad(angle)))
            s.RadialSpeed_mps = float(speed * np.cos(np.deg2rad(angle)))
        for name, field in [('Direction_deg', 'Direction_deg'), ('Range_m', 'Range_m'),
                            ('InitialAngles_deg', 'Angles_deg')]:
            if name in values:
                s[field] = values[name][j].tolist()
    return c


def run_example(c, progress=None):
    """Independent result saving, figure saving and display, as in RUN_SPOT."""
    from .pipeline import run
    from .plotting import plot_product
    import matplotlib.pyplot as plt

    out = run(c, export=True, progress=progress)
    c = out.Configuration
    if not (c.Output.Plot or c.Plot.Save):
        return out
    indices = range(len(out.Exposures)) if c.Plot.Exposures == 'all' else np.atleast_1d(c.Plot.Exposures) - 1
    products = list(c.Plot.Products) + (['photometry'] if c.Plot.Photometry else [])
    folder = None
    if c.Plot.Save:
        base = Path(out.ExportDirectory) if out.get('ExportDirectory') else Path(c.Output.Directory) / (datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f') + '_' + c.Dataset.SampleID)
        folder = base / 'Figures'
        folder.mkdir(parents=True, exist_ok=False)
        out.FigureDirectory = str(folder)
    for k in indices:
        if int(k) != k or not 0 <= k < len(out.Exposures):
            raise ValueError('Plot.Exposures uses one-based indices or "all"')
        for kind in products:
            fig = plot_product(out, kind, int(k))
            if folder:
                for extension in c.Plot.Formats:
                    if extension not in ('png', 'pdf', 'svg'):
                        raise ValueError('Python figure formats: png, pdf, svg')
                    fig.savefig(folder / f'exposure_{int(k)+1:04d}_{kind}.{extension}', dpi=c.Plot.DPI)
            if c.Output.Plot:
                plt.show()
            plt.close(fig)
    return out


def batch(configurations, progress=None):
    """Process acquisitions sequentially and return their index without retaining images."""
    rows = []
    for c in configurations:
        out = run_example(c, progress=progress)
        rows.append(dict(SampleID=c.Dataset.SampleID, Split=c.Dataset.Split,
                         Directory=out.get('ExportDirectory', ''), NumExposures=len(out.Exposures)))
        del out
    return pd.DataFrame(rows)


def dataset_configurations(c, count=10, first_selection_seed=100):
    for n in range(1, count + 1):
        item = deepcopy(c)
        item.Input.TargetSelection.Seed = first_selection_seed + n
        item.RandomSeed = n
        item.Dataset.SampleID = f'sample_{n:04d}'
        item.Dataset.Split = 'train'
        yield item
