from copy import deepcopy
from pathlib import Path
import json
import numpy as np
import pytest
from ssa_emulator.config import default_config
from ssa_emulator.examples import example_config, apply_example_overrides, run_example
from ssa_emulator.inputs import resolve_pools
from ssa_emulator.imaging import labels
from ssa_emulator.config import struct
from ssa_emulator.frigate_conditioning import configure_spot_for_detector_maps, apply_spot_to_conditioned_maps
from ssa_emulator.validation import direct_pair_metrics, target_integrated_flux, load_external_satsim


@pytest.mark.parametrize('name',['ground','space_sky','space_earth'])
def test_final_matlab_defaults(name):
    c=default_config(name)
    assert (c.Time.Exposure_s,c.Time.Cadence_s)==(1.5,1.5)
    assert c.Camera.Nx==c.Camera.Ny==1024
    assert c.Stars.Source=='csv' and c.Stars.MaxStars==648000
    assert c.Frames.Model=='iau2006' and c.Detector.PoissonEnabled
    assert c.Radiometry.DetectorCalibration=='band_photon'


def test_example_schedules():
    c=example_config('multi_exposure');assert (c.Time.Exposure_s,c.Time.Cadence_s,c.Time.NumExposures)==(1.5,1,3)
    c=example_config('dataset');assert (c.Time.Exposure_s,c.Time.Cadence_s,c.Input.NumTargets)==(1,.8,20)
    c=example_config('csv_pools');assert (c.Time.Exposure_s,c.Time.Cadence_s,c.Input.NumTargets,c.Input.TargetSelection.StartRow)==(1,1,5,10)


def test_empty_overrides_preserve_inputs():
    c=example_config('ground');before=deepcopy(c)
    apply_example_overrides(c)
    assert c.Targets==before.Targets and c.Pointing==before.Pointing


def test_csv_override_precedence_and_speed():
    c=example_config('csv_pools');resolve_pools(c)
    speed=np.array([np.hypot(t.Seed.TransverseSpeed_mps,t.Seed.RadialSpeed_mps) for t in c.Targets])
    c.ExamplePointing.Roll_deg=17;c.ExampleTarget.AngleToLOS_deg=[90,.2,90,.2,90]
    apply_example_overrides(c)
    assert c.Input.SelectedTargetRows==[10,11,12,13,14] and c.Pointing.Roll_deg==17
    np.testing.assert_allclose([t.Seed.TransverseSpeed_mps for t in c.Targets],speed*np.sin(np.deg2rad([90,.2,90,.2,90])),atol=1e-8)
    c.ExampleTarget.AngleToLOS_deg=[90,0]
    with pytest.raises(ValueError):apply_example_overrides(c)


def test_nonseed_override_rejected():
    c=example_config('orbit_tracking');c.ExampleTarget=struct(dict(AngleToLOS_deg=90))
    with pytest.raises(ValueError):apply_example_overrides(c)


@pytest.mark.parametrize('save_results,save_plots',[(False,False),(True,False),(False,True),(True,True)])
def test_save_switches(tmp_path,save_results,save_plots):
    c=example_config('ground','preview');c.Camera.Nx=16;c.Camera.Ny=12;c.Targets=[];c.Stars.Enabled=False
    c.Time.Exposure_s=c.Time.Cadence_s=.02;c.Time.NumExposures=2;c.Atmosphere.Enabled=False
    c.Output.Directory=str(tmp_path/'runs');c.Output.Plot=False;c.Output.SaveResults=save_results
    c.Plot.Save=save_plots;c.Plot.Products=['detector'];c.Plot.Exposures='all';c.Plot.Formats=['png'];c.Plot.DPI=50
    out=run_example(c)
    assert bool(out.get('ExportDirectory'))==save_results
    assert bool(out.get('FigureDirectory'))==save_plots
    assert len(list(tmp_path.rglob('metadata.json')))==int(save_results)
    assert len(list(tmp_path.rglob('*_detector.png')))==2*int(save_plots)
    if not save_results and not save_plots:assert not (tmp_path/'runs').exists()


def test_matlab_mask_metadata_and_tie_order():
    c=default_config();c.Camera.Nx=c.Camera.Ny=4;a=np.zeros((4,4));a[1,2]=4
    records=[struct(dict(ID=8,Appearance='object')),struct(dict(ID=2,Appearance='trail'))]
    m=labels(c,[a,a],records)
    assert m.Instance[1,2]==8 and m.Semantic[1,2]==1
    assert m.ObjectUnion[1,2] and m.TrailUnion[1,2]
    assert records[0].MaskThreshold==.01 and records[0].AbsoluteMaskThreshold==.04
    assert records[0].BoundingBox_xywh==[3,2,1,1]


def test_conditioned_adapter_uses_production_stages(monkeypatch):
    from ssa_emulator import imaging
    calls=[];old_optics=imaging.optics;old_detector=imaging.detector
    def optics(*a):calls.append('optics');return old_optics(*a)
    def detector(*a):calls.append('detector');return old_detector(*a)
    monkeypatch.setattr(imaging,'optics',optics);monkeypatch.setattr(imaging,'detector',detector)
    c=configure_spot_for_detector_maps(default_config(),(16,16),exposure_s=.5,deterministic=True)
    point=np.zeros((16,16));point[8,8]=30;bg=np.full((16,16),120.)
    a=apply_spot_to_conditioned_maps(c,point,bg,.5)
    assert calls==['optics','optics','detector']
    assert np.min(a['Image'])>=120 and a['Image'][8,8]>120
    c.Time.Exposure_s=1.5
    b=apply_spot_to_conditioned_maps(c,point,bg,1.5)
    np.testing.assert_allclose(a['Image'],b['Image'],atol=1e-5)


def test_validation_metrics_and_positive_flux():
    a=np.full((16,16),100.);b=a+2;metrics=direct_pair_metrics(a,b)
    assert metrics['RMSE_DN']==metrics['MAE_DN']==2
    assert metrics['PSNR_dB']==pytest.approx(20*np.log10(102/2))
    support=np.zeros_like(a,bool);support[4:6,4:6]=True
    a[4,4]=110;a[4,5]=90
    assert target_integrated_flux(a,support)==(10.,100.)


def test_external_satsim_provenance_and_shape(tmp_path):
    p=tmp_path/'s.npy';np.save(p,np.ones((16,16)))
    m=tmp_path/'provenance.json';m.write_text(json.dumps(dict(simulator='SatSim',version='test',noise_on=False,downsample=4,units='DN',scene_config={'test_fixture':True})))
    a,meta=load_external_satsim(p,m,(16,16),False)
    assert a.shape==(16,16) and len(meta['image_sha256'])==64
    with pytest.raises(ValueError):load_external_satsim(p,m,(32,32),False)
    with pytest.raises(ValueError):load_external_satsim(p,m,(16,16),True)
    with pytest.raises(ValueError):load_external_satsim(p,None,(16,16),False)
