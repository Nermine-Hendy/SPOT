"""Compare deterministic MATLAB fixtures exported by Export_MATLAB_Parity.m.

Run from the python directory:
    python Verification/compare_matlab.py Verification/matlab_fixtures
"""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ssa_emulator.config import default_config,struct
from ssa_emulator.pipeline import run


def config_from_matlab(data):
    c=default_config(data['Name'],'full')
    def merge(a,b):
        for key,value in b.items():
            if isinstance(value,dict) and isinstance(a.get(key),dict):merge(a[key],value)
            else:a[key]=value
    merge(c,data);c=struct(c);c.RootDir=str(ROOT)
    if isinstance(c.Targets,dict):c.Targets=[c.Targets]
    if c.Visibility.RangeLimits_m[1] is None:c.Visibility.RangeLimits_m[1]=float('inf')
    c.Stars.File=str(ROOT/'Inputs/Stars_648000.csv')
    for key in ('SensorCSV','TargetCSV'):
        old=str(c.Input[key]).replace('\\','/');c.Input[key]=str(ROOT/'Inputs'/old.split('/')[-1])
    c.Output.SaveResults=c.Output.Plot=c.Plot.Save=False
    return c


def compare(folder,rtol=1e-5,atol=1e-9):
    files=sorted(Path(folder).glob('*.mat'))
    if not files:raise FileNotFoundError('No genuine MATLAB fixtures found. Run Export_MATLAB_Parity.m first.')
    rows=[];cached={}
    for file in files:
        fixture=loadmat(file,simplify_cells=True);raw=str(fixture['config_json'])
        if raw not in cached:cached[raw]=run(config_from_matlab(json.loads(raw)))
        E=cached[raw].Exposures[int(fixture['exposure_index'])-1]
        actual={k:E[k] for k in ('Image_DN','SceneRate','TargetRate','StarRate','SkyRate','EarthRate')}
        actual.update(OpticalSignal=E.OpticalStages.PSF,ExpectedSignal_e=E.DetectorStages.ExpectedSignal_e)
        actual.update({k:E.Labels[k] for k in ('Binary','Instance','Semantic')})
        actual.update({k:E.Truth[k] for k in ('Xi_deg','Eta_deg','RA_deg','DEC_deg','EarthMask')})
        for key,expected in fixture['F'].items():
            a=np.asarray(actual[key]);b=np.asarray(expected)
            exact=key in ('Binary','Instance','Semantic','EarthMask')
            equal=np.array_equal(a,b) if exact else np.allclose(a,b,rtol=rtol,atol=atol,equal_nan=True)
            delta=np.abs(a.astype(float)-b.astype(float))
            rows.append(dict(fixture=file.name,field=key,passed=bool(equal),exact=exact,
                             max_absolute_error=float(np.nanmax(delta)),
                             relative_l1=float(np.nansum(delta)/max(np.nansum(abs(b)),1e-300))))
    return dict(passed=all(r['passed'] for r in rows),rtol=rtol,atol=atol,comparisons=rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('folder',type=Path)
    parser.add_argument('--rtol',type=float,default=1e-5);parser.add_argument('--atol',type=float,default=1e-9)
    args=parser.parse_args();report=compare(args.folder,args.rtol,args.atol)
    path=ROOT/'Verification/matlab_parity_results.json';path.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('PASS' if report['passed'] else 'FAIL',path)
    sys.exit(0 if report['passed'] else 1)
