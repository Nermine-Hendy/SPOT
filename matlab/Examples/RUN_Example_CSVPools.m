% Ordered CSV selection: five targets beginning at row 10.
% Run this file from its folder or use run(fullfile(root,'Examples',filename)).
clear;clc;close all hidden
addpath(fileparts(fileparts(mfilename('fullpath'))));
SPOT_Setup;
C=SPOT_Config('current_ground','full'); % FULL: CSV catalogue, 1024 px, frames and noise.
%% EDIT HERE: exposure, target motion, display and saving
E=struct;
E.Exposure_s=C.Time.Exposure_s; % Default: 1 s. Exposure duration.
E.Cadence_s=1; % Default CSV launcher sets start-to-start cadence to 1 s.
E.NumExposures=C.Time.NumExposures; % Default: 1.

% [] means keep the configured target setting (or the selected CSV value).
% AngleToLOS: 0 = receding along LOS, 90 = transverse, 180 = approaching.
% Small angles reduce initial transverse motion; exposure/pointing/PSF also matter.
% Scalars apply to every target; vectors must match the selected target count.
E.Target.AngleToLOS_deg=[]; % Example override: [90 0.2 90 0.2 90]
E.Target.RelativeSpeed_mps=[]; % [] preserves configured speed magnitude.
E.Target.Direction_deg=[];     % Motion direction from +xi toward +eta [deg].
E.Target.InitialAngles_deg=[]; % Initial source position; one [angle1 angle2] row per target.
E.Target.Range_m=[];           % Initial observer-to-target range [m].
% Empty pointing overrides retain the selected sensor row / base observer.
E.Pointing.AzEl_deg=[];
E.Pointing.OffNadir_deg=[];
E.Pointing.SpaceAzimuth_deg=[];
E.Pointing.Roll_deg=[];

E.ShowPlots=true;   % Open figures on screen.
E.SavePlots=true;   % Save figure PNG/FIG files independently of results.
E.SaveResults=false; % Save acquisition MAT, metadata JSON, target CSV and image/mask PNG.
E.PlotProducts={'camera','detector','radec','stages','dataset'};
E.PlotPhotometry=false; % Optional diffuse/specular and received-irradiance histories.
E.PlotExposures=1; % First exposure; use 'all' if NumExposures is increased.

C.Input.UseSensorCSV=true;
C.Input.UseTargetCSV=true;
C.Input.NumTargets=5;
C.Input.TargetSelection.Mode='ordered'; % 'random' and 'ids' are also supported.
C.Input.TargetSelection.StartRow=10;
C.Input.TargetSelection.Seed=100;
% For random selection: C.Input.TargetSelection.Mode='random';
% For explicit IDs, also set NumTargets to the number of IDs:
% C.Input.NumTargets=2; C.Input.TargetSelection.Mode='ids';
% C.Input.TargetSelection.IDs=["SEED0003" "SEED0001"];
% Default pool is TargetPool_Seeded.csv. For Keplerian targets, leave E.Target fields empty.
% C.Input.TargetCSV=fullfile(C.RootDir,'Inputs','TargetPool_Keplerian.csv');

%% APPLY AND RUN (normally leave this section unchanged)
C=SPOT_ApplyExample(C,E);
Out=RUN_SPOT(C);
