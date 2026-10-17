from pathlib import Path
from datetime import datetime,timezone
import json
import numpy as np
import pandas as pd
from .config import plain

def json_safe(x):
    x=plain(x)
    if isinstance(x,dict):return {k:json_safe(v) for k,v in x.items()}
    if isinstance(x,list):return [json_safe(v) for v in x]
    if isinstance(x,float) and not np.isfinite(x):return None
    return x

def array_products(E):
    a={k:np.asarray(E[k]) for k in ['Image_DN','SceneRate','TargetRate','StarRate','SkyRate','EarthRate','DiffuseRate','TotalTargetOpticalSignal','TargetTime_s','TargetWeights','BackgroundTime_s','BackgroundWeights']}
    for parent in ['Labels','DetectorStages','OpticalStages','Truth']:
        for k,v in E[parent].items():
            if isinstance(v,np.ndarray):a[parent+'__'+k]=v
    for k,v in E.Truth.AngularDiagnosticGrid.items():a['Truth__AngularDiagnosticGrid__'+k]=v
    for j,T in enumerate(E.Targets):
        for k,v in T.Samples.items():
            if isinstance(v,np.ndarray):a[f'Target{j+1}__'+k]=v
        if 'Reflection' in T.Samples:
            for k,v in T.Samples.Reflection.items():a[f'Target{j+1}__Reflection__'+k]=v
    for j,m in enumerate(E.TargetOpticalSignal):a[f'Target{j+1}__OpticalSignal']=m
    for j,m in enumerate(E.Labels.PerTarget):a[f'Target{j+1}__Mask']=m
    return a

def export_run(out,folder=None):
    c=out.Configuration
    if not c.Output.get('SaveResults',True):
        return None
    if folder is None:folder=Path(c.Output.Directory)/(datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f')+'_'+c.Dataset.SampleID)
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=False)
    meta=dict(SchemaVersion=out.SchemaVersion,Configuration=c,Provenance=out.Provenance,CoordinateConvention='One-based x/y pixel centres and inclusive xywh boxes; arrays row-major Ny by Nx',Exposures=[])
    for k,E in enumerate(out.Exposures,1):
        if c.Output.SaveNPZ:np.savez_compressed(folder/f'exposure_{k:04d}.npz',**array_products(E))
        records=[{key:v for key,v in T.items() if key not in ['Samples','Morphology']} for T in E.Targets]
        if c.Output.SaveCSV:pd.DataFrame(records).to_csv(folder/f'exposure_{k:04d}_targets.csv',index=False)
        if c.Output.SavePNG:
            from PIL import Image
            lo,hi=np.percentile(E.Image_DN,c.Output.PNGPercentiles);im=np.clip((E.Image_DN-lo)/max(hi-lo,1e-12),0,1);Image.fromarray(np.uint8(im*255)).save(folder/f'exposure_{k:04d}_quicklook.png')
            Image.fromarray(E.Labels.Semantic).save(folder/f'exposure_{k:04d}_semantic.png');Image.fromarray(np.uint8(E.Labels.Binary)*255).save(folder/f'exposure_{k:04d}_binary.png')
            Image.fromarray(np.uint8(E.DetectorStages.ValidPixel)*255).save(folder/f'exposure_{k:04d}_valid_pixels.png')
            if int(E.Labels.Instance.max())<=65535:
                Image.fromarray(E.Labels.Instance.astype(np.uint16)).save(folder/f'exposure_{k:04d}_instance.png')
        meta['Exposures'].append(dict(Index=k,StartUTC=E.StartUTC,Regime=E.Truth.GeometryRegime,Sampling=E.Sampling,Targets=records))
    (folder/'metadata.json').write_text(json.dumps(json_safe(meta),indent=2,allow_nan=False))
    for field,name in [('SelectedTargets','selected_targets'),('SelectedSensor','selected_sensor')]:
        if field in c.Input:
            rows=c.Input[field];rows=[rows] if isinstance(rows,dict) else rows;pd.DataFrame(rows).to_csv(folder/(name+'.csv'),index=False)
    return folder
