"""Regression tests for the sunlit Earth-limb preset."""
from copy import deepcopy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pytest

from ssa_emulator.config import default_config, struct
from ssa_emulator.examples import example_config
from ssa_emulator.inputs import validate
from ssa_emulator.geometry import sample
from ssa_emulator.background import ray_scene
from ssa_emulator.plotting import plot_product


@pytest.mark.parametrize('factory', [default_config, example_config])
def test_public_limb_reference_acquisition(factory):
    c = factory('space_limb', 'full')
    assert c.Time.StartUTC == '2026-05-01T12:00:00Z'
    assert (c.Time.Exposure_s, c.Time.Cadence_s) == (1, .8)
    assert c.Observer.Trajectory.Elements == [6878137, .001, 97.5, 40, 0, 0]
    assert c.Stars.MaxStars == 200000 and c.Stars.Source == 'csv'
    assert c.Camera.Nx == c.Camera.Ny == 1024
    assert c.Frames.Model == 'iau2006' and c.Detector.PoissonEnabled
    assert c.Radiometry.DetectorCalibration == 'legacy_gain'
    assert c.Detector.RadiometricGain_e_per_unit == 1e6
    assert c.Seeding.Propagation == 'two_body'
    assert c.Integration.AttitudeMode == 'dynamic'
    assert all(t.Photometry.Model == 'surface_brdf' for t in c.Targets)
    angles = [t.Seed.Angles_deg for t in c.Targets]
    assert angles == [[8, 20], [10, 90], [14, 160], [12, 230]]
    speed = [np.hypot(t.Seed.TransverseSpeed_mps, t.Seed.RadialSpeed_mps) for t in c.Targets]
    np.testing.assert_allclose(speed, [7500, 7500, 6000, 8000])
    assert c.Plot.DetectorColorbar


def test_limb_brightness_is_physical_and_extends_to_both_sides():
    c = default_config('space_limb', 'full')
    c.Camera.Nx = c.Camera.Ny = 128
    c.Targets = []
    c = validate(c)
    d, rate = ray_scene(c, sample(c, [.5], .5), 0)
    assert d.GeometryRegime == 'space_limb'
    assert .42 < d.EarthFraction < .46
    assert d.EarthRadiance.max() > 1000
    # Both sides of the geometric boundary contain the atmospheric rim.
    near = abs(d.LimbDistance_deg) < .4
    assert np.median(d.EarthRadiance[near & d.EarthMask]) > 800
    assert np.median(d.EarthRadiance[near & ~d.EarthMask]) > 800
    assert np.max(d.EarthRadiance[d.LimbDistance_deg > 8]) == 0
    np.testing.assert_allclose(rate, d.DiffuseRadiance*d.SolidAngle_sr)
    # A night-side epoch and orbit orientation provide a dark-limb comparison.
    night = deepcopy(c)
    night.Time.StartUTC = night.Observer.Trajectory.EpochUTC = '2026-01-01T12:00:00Z'
    night.Observer.Trajectory.Elements = [6878137, .001, 97.5, 20, 0, 0]
    dark, _ = ray_scene(night, sample(night, [.5], .5), 0)
    assert dark.EarthRadiance.max() < 2


def test_limb_plot_uses_real_dn_and_does_not_modify_data():
    c = default_config('space_limb')
    c.Camera.Nx = c.Camera.Ny = 4
    image = np.arange(16, dtype=float).reshape(4, 4)*30+100
    e = struct(dict(Image_DN=image, Targets=[], Index=1,
                    Truth=dict(GeometryRegime='space_limb'),
                    Labels=dict(Binary=np.zeros((4, 4)), Instance=np.zeros((4, 4)),
                                Semantic=np.zeros((4, 4))), TargetRate=np.zeros((4, 4))))
    out = struct(dict(Configuration=c, Exposures=[e]))
    original = image.copy()
    fig = plot_product(out, 'detector')
    assert len(fig.axes) == 2
    np.testing.assert_array_equal(fig.axes[0].images[0].get_array(), image)
    np.testing.assert_allclose(fig.axes[0].images[0].get_clim(), np.percentile(image, [1, 99.7]))
    assert fig.axes[1].get_ylabel() == 'Digital number (DN)'
    np.testing.assert_array_equal(image, original)
    plt.close(fig)
    # Adding the standalone DN bar must not alter the six-panel product layout.
    fig = plot_product(out, 'dataset')
    assert len(fig.axes) == 6
    plt.close(fig)


def test_legacy_limb_and_user_overrides_remain_available():
    legacy = default_config('current_space_limb')
    assert legacy.Seeding.Propagation == 'linear'
    assert legacy.Integration.AttitudeMode == 'midpoint'
    c = example_config('space_limb', 'preview')
    c.Time.Exposure_s = .7
    c.ExamplePointing.Roll_deg = 12
    c.ExampleTarget.AngleToLOS_deg = [90, 1, 90, 90]
    c = validate(c)
    assert c.Time.Exposure_s == .7 and c.Pointing.Roll_deg == 12
    assert c.Camera.Nx == c.Camera.Ny == 256
