% Dataset generation: current_ground, 20 pool targets per sample.
% Run this file from its folder or use run(fullfile(root,'Examples',filename)).
clear;clc;close all hidden
addpath(fileparts(fileparts(mfilename('fullpath'))));
SPOT_Setup;
C=SPOT_Config('current_ground','full'); % FULL: CSV catalogue, 1024 px, frames and noise.
%% EDIT HERE: exposure, target motion, display and saving
E=struct;
E.Exposure_s=C.Time.Exposure_s; % Default: 1 s. Exposure duration.
E.Cadence_s=C.Time.Cadence_s; % Default: 0.8 s between exposure starts.
E.NumExposures=C.Time.NumExposures; % Default: 1.

% [] means keep the configured target setting (or the selected CSV value).
% AngleToLOS: 0 = receding along LOS, 90 = transverse, 180 = approaching.
% Small angles reduce initial transverse motion; exposure/pointing/PSF also matter.
% Scalars apply to every target; vectors must match the selected target count.
E.Target.AngleToLOS_deg=[]; % Example override: 90
E.Target.RelativeSpeed_mps=[]; % [] preserves configured speed magnitude.
E.Target.Direction_deg=[];     % Motion direction from +xi toward +eta [deg].
E.Target.InitialAngles_deg=[]; % Initial source position; one [angle1 angle2] row per target.
E.Target.Range_m=[];           % Initial observer-to-target range [m].
% Empty pointing overrides retain the selected sensor row / base observer.
E.Pointing.AzEl_deg=[];
E.Pointing.OffNadir_deg=[];
E.Pointing.SpaceAzimuth_deg=[];
E.Pointing.Roll_deg=[];

E.ShowPlots=false;   % Open figures on screen.
E.SavePlots=false;   % Save figure PNG/FIG files independently of results.
E.SaveResults=true; % Save acquisition MAT, metadata JSON, target CSV and image/mask PNG.
E.PlotProducts={'camera','detector','radec','stages','dataset'};
E.PlotPhotometry=false; % Optional diffuse/specular and received-irradiance histories.
E.PlotExposures=1; % First exposure; use 'all' if NumExposures is increased.

NumSamples=10;
C.Input.UseSensorCSV=false; % Keep the configured ground observer; true selects a sensor row.
C.Input.UseTargetCSV=true;
C.Input.NumTargets=20;
C.Input.TargetSelection.Mode='random';
FirstSelectionSeed=100; % Sample n uses seed 100+n.
% Default pool is TargetPool_Seeded.csv. For Keplerian targets, leave E.Target fields empty.
% C.Input.TargetCSV=fullfile(C.RootDir,'Inputs','TargetPool_Keplerian.csv');

%% APPLY AND RUN (normally leave this section unchanged)
C=SPOT_ApplyExample(C,E);
configurations=cell(1,NumSamples);
for n=1:NumSamples
    D=C; D.Input.TargetSelection.Seed=FirstSelectionSeed+n; D.RandomSeed=n;
    D.Dataset.SampleID=sprintf('sample_%04d',n); D.Dataset.Split='train';
    configurations{n}=D;
end
index=SPOT.batch(configurations);
if E.SaveResults
    writetable(index,fullfile(C.Output.Directory,'dataset_index.csv'));
end
% index remains in memory even when saving is disabled.
