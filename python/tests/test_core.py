from copy import deepcopy
from pathlib import Path
import json
import numpy as np
import pytest
from ssa_emulator import default_config,run,convergence
from ssa_emulator.config import ROOT,struct,save_config,load_config
from ssa_emulator.geometry import geodetic,state,project,rx,ry,rz,rotation,iso
from ssa_emulator.photometry import surface_reflection,target_attitude,atmosphere
from ssa_emulator.imaging import detector_grid,deposit,optics,detector,labels
from ssa_emulator.sampling import time_grid
from ssa_emulator.inputs import validate,resolve_pools
from ssa_emulator.export import export_run

def small(name='ground',size=32):
 c=default_config(name,'preview');c.Camera.Nx=c.Camera.Ny=size;c.Stars.Enabled=False;c.Atmosphere.Scintillation.Enabled=False;c.Detector.ReadEnabled=False;return c

def test_projection_roundtrip():
 c=small();A=rz(33)@ry(72);d=detector_grid(c,A);p=project(c,A,d.LOS_GCRF);y,x=np.indices(d.Active.shape)
 np.testing.assert_allclose(p.x.reshape(d.Active.shape,order='F'),x+1,atol=1e-12);np.testing.assert_allclose(p.y.reshape(d.Active.shape,order='F'),y+1,atol=1e-12);assert p.InFOV.all()
 assert abs(d.SolidAngle_sr.sum()-d.FOVSolidAngle_sr)<1e-14

def test_geodetic_roundtrip():
 c=small();v=np.array([[-37.8,0,89.99],[145,-170,20],[50,5e5,100]])
 np.testing.assert_allclose(geodetic(c,geodetic(c,v),True),v,atol=1e-5)

def test_orbit_conservation():
 c=small('space_sky');tr=c.Observer.Trajectory;period=2*np.pi*np.sqrt(tr.Elements[0]**3/c.Constants.Mu_m3_s2);s=state(c,tr,[0,period/3,period])
 np.testing.assert_allclose(s.Position_m[:,0],s.Position_m[:,-1],atol=1e-7)
 energy=(s.Velocity_mps**2).sum(axis=0)/2-c.Constants.Mu_m3_s2/np.linalg.norm(s.Position_m,axis=0);np.testing.assert_allclose(energy,energy[0],rtol=1e-12)

def test_two_body_matches_kepler():
 c=small('space_sky');tr=deepcopy(c.Observer.Trajectory);s=state(c,tr,[0]);tr.Type='state';tr.Motion='two_body';tr.Position_m=s.Position_m[:,0];tr.Velocity_mps=s.Velocity_mps[:,0]
 a=state(c,c.Observer.Trajectory,[-20,0,20]);b=state(c,tr,[-20,0,20]);np.testing.assert_allclose(a.Position_m,b.Position_m,atol=.01,rtol=0)

def test_surface_equations():
 S=struct(dict(NormalBody=[[0],[0],[1]],Area_m2=[2],DiffuseReflectance=[.4],SpecularReflectance=[0],SpecularExponent=[40]));z=np.array([[0.],[0.],[1.]]);R=np.eye(3)[None]
 a=surface_reflection(S,R,z,z,[500],[1e6]);np.testing.assert_allclose(a.Irradiance_W_m2,500*2*.4/np.pi/1e12,rtol=1e-14)
 b=surface_reflection(S,R,z,z,[500],[2e6]);np.testing.assert_allclose(b.Irradiance_W_m2,a.Irradiance_W_m2/4)
 assert surface_reflection(S,R,-z,z,[500],[1e6]).Irradiance_W_m2[0]==0
 S.DiffuseReflectance=[0];S.SpecularReflectance=[.2];s=ry(30)@z;v=ry(-30)@z
 a=surface_reflection(S,R,s,v,[500],[1e6]);b=surface_reflection(S,R,s,s,[500],[1e6]);assert b.Irradiance_W_m2[0]<a.Irradiance_W_m2[0]*1e-8

def test_reflectance_energy_validation():
 c=small();c.Targets[0].Surface.DiffuseReflectance=[1]*6
 with pytest.raises(ValueError):validate(c)

def test_attitude_spin():
 c=small();a=c.Targets[0].Attitude;a.SpinRate_deg_s=90;R=target_attitude(c,a,[0,1]);np.testing.assert_allclose(R[1]@[1,0,0],[0,1,0],atol=1e-12)

def test_quadrature():
 t,w=time_grid(0,.23,.07);assert t[-1]==.23;assert abs(w.sum()-1)<1e-14;assert abs(w@t-.115)<1e-14

def test_deposition_edges_and_mask():
 c=small();active=np.ones((32,32),bool);active[10,10]=False;a=deposit(c,[1,32,11.2],[1,32,11.2],[2,3,4],active);assert abs(a.sum()-9)<1e-12;assert a[10,10]==0

def test_optics_linearity():
 c=small(size=64);c.Time.Exposure_s=1;c.Optics.Transmission=1;c.Optics.VignettingEnabled=False;a=np.zeros((64,64));a[32,32]=1;b=np.roll(a,1,axis=1);active=np.ones_like(a,bool)
 x=optics(c,a,active)[0];y=optics(c,b,active)[0];z=optics(c,a+b,active)[0];np.testing.assert_allclose(x+y,z,atol=1e-15);assert abs(x.sum()-1)<1e-12

def test_detector_order():
 c=small();c.Detector.PoissonEnabled=False;c.Detector.DarkEnabled=False;c.Detector.ReadEnabled=False;c.Detector.RadiometricGain_e_per_unit=1;c.Detector.FullWell_e=5;c.Detector.Bias_DN=10;c.Detector.ADCEnabled=True;c.Detector.ADCBits=4
 a=detector(c,np.ones((32,32))*7,1,np.ones((32,32),bool));assert np.all(a.AfterFullWell_e==5);assert np.all(a.Image_DN==15);assert a.FullWellSaturated.all()

def test_mask_overlap():
 c=small();a=np.zeros((32,32));a[10,10]=2;b=a/2;r=[struct(dict(ID=5,Appearance='trail')),struct(dict(ID=3,Appearance='object'))];L=labels(c,[a,b],r);assert L.Instance[10,10]==5;assert L.Semantic[10,10]==2

def test_pool_selection():
 c=small('current_ground');c.Input.UseTargetCSV=True;c.Input.NumTargets=5;c.Input.TargetSelection.StartRow=10;a=resolve_pools(deepcopy(c));assert a.Input.SelectedTargetRows==[10,11,12,13,14]
 c.Input.TargetSelection.Mode='random';c.Input.TargetSelection.Seed=100;a=resolve_pools(deepcopy(c));b=resolve_pools(deepcopy(c));assert a.Input.SelectedTargetRows==b.Input.SelectedTargetRows;assert len(set(a.Input.SelectedTargetRows))==5
 c.Input.TargetSelection.Mode='ids';c.Input.NumTargets=2;c.Input.TargetSelection.IDs=['SEED0003','SEED0001'];a=resolve_pools(c);assert a.Input.SelectedTargetRows==[3,1]

def test_scintillation_grid_invariance():
 c=small();c.Atmosphere.Scintillation.Enabled=True
 _,a=atmosphere(c,np.array([0,.1,.2]),np.full(3,45),np.ones(3,bool),2);_,b=atmosphere(c,np.array([0,.05,.1,.15,.2]),np.full(5,45),np.ones(5,bool),2);np.testing.assert_array_equal(a,b[::2])

@pytest.mark.parametrize('case,regime',[('ground','ground_sky'),('space_sky','space_sky'),('space_limb','space_limb'),('space_earth','space_earth')])
def test_regimes(case,regime):
 o=run(small(case));assert o.Exposures[0].Truth.GeometryRegime==regime;assert np.isfinite(o.Exposures[0].Image_DN).all()

def test_figure8_appearances():
 c=small('figure8',256);o=run(c);assert [r.Appearance for r in o.Exposures[0].Targets]==['trail','object','trail','unobserved'];assert [r.HasAnnotation for r in o.Exposures[0].Targets]==[True,True,True,False]
 c.Time.Exposure_s=.001;o=run(c);assert o.Exposures[0].Targets[0].Appearance=='object'

def test_zero_targets_and_exports(tmp_path):
 c=small();c.Targets=[];o=run(c);assert not o.Exposures[0].Labels.Binary.any();p=export_run(o,tmp_path/'run');z=np.load(p/'exposure_0001.npz');assert 'Image_DN' in z;assert json.loads((p/'metadata.json').read_text())['Provenance']['Reference']=='MATLAB ZPublic-SPOT 7.2.0 (schema 7.0.0)'

@pytest.mark.parametrize('name',['trajectory_geodetic','trajectory_ecef','trajectory_state','trajectory_ephemeris','trajectory_ephemeris_ecef','trajectory_ephemeris_csv','trajectory_keplerian','trajectory_tle','trajectory_tle_deep_space'])
def test_trajectory_inputs(name):
 c=small(name);c.Time.Exposure_s=.2;o=run(c);assert np.isfinite(o.Exposures[0].Image_DN).all()

def test_config_roundtrip(tmp_path):
 c=default_config('figure8','preview');p=tmp_path/'c.json';save_config(c,p);d=load_config(p);assert d.Targets[0].Surface==c.Targets[0].Surface

def test_convergence():
 c=small('figure8',64);r=convergence(c);assert max(r[0]['TargetFluxRelativeChange'])<.01

def test_iau_rotation():
 c=small();c.Frames.Model='iau2006';R=rotation(c,[0,.2]);np.testing.assert_allclose(R[0].T@R[0],np.eye(3),atol=1e-12)

def test_gui_control_merge_preserves_advanced_edits():
 from gui import SPOTApp
 class Variable:
  def __init__(self,s):self.s=s
  def get(self):return self.s
 c=small();c.Time.Exposure_s=.75;c.Time.Cadence_s=1
 app=object.__new__(SPOTApp);app.get_config=lambda:c
 app.variables={('Time','Exposure_s'):Variable('.2'),('Time','Cadence_s'):Variable('2')}
 app.form_values={('Time','Exposure_s'):'.2',('Time','Cadence_s'):'1'}
 app.set_config=lambda value:None
 result=app.apply_controls()
 assert result.Time.Exposure_s==.75
 assert result.Time.Cadence_s==2
