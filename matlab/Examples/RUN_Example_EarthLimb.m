% Sunlit Earth-limb observation with finite-exposure integration.
% Run this file from its folder or use run(fullfile(root,'Examples',filename)).
clear;clc;close all hidden
addpath(fileparts(fileparts(mfilename('fullpath'))));
SPOT_Setup;
C=SPOT_Config('space_limb','full'); % FULL: CSV catalogue, 1024 px, frames and noise.
%% EDIT HERE: exposure, target motion, display and saving
E=struct;
E.Exposure_s=C.Time.Exposure_s; % Reference: 1 s. Exposure duration.
E.Cadence_s=C.Time.Cadence_s; % Reference: 0.8 s between exposure starts.
E.NumExposures=C.Time.NumExposures; % Default: 1.

% [] means keep the configured target setting (or the selected CSV value).
% AngleToLOS: 0 = receding along LOS, 90 = transverse, 180 = approaching.
% Small angles reduce initial transverse motion; exposure/pointing/PSF also matter.
% Scalars apply to every target; vectors must match the selected target count.
E.Target.AngleToLOS_deg=[]; % Example override: [90 0.2 90 0.2]
E.Target.RelativeSpeed_mps=[]; % [] preserves configured speed magnitude.
E.Target.Direction_deg=[];     % Motion direction from +xi toward +eta [deg].
E.Target.InitialAngles_deg=[]; % Initial source position; one [angle1 angle2] row per target.
E.Target.Range_m=[];           % Initial observer-to-target range [m].
% InitialAngles_deg above means [off-boresight theta, around-boresight phi] for space seeds.
E.Pointing.OffNadir_deg=C.Pointing.OffNadir_deg; % Default 70 deg.
E.Pointing.SpaceAzimuth_deg=C.Pointing.SpaceAzimuth_deg; % Default 0 deg.
E.Pointing.Roll_deg=C.Pointing.Roll_deg; % Default 0 deg.

E.ShowPlots=true;   % Open figures on screen.
E.SavePlots=false;   % Save figure PNG/FIG files independently of results.
E.SaveResults=false; % Save acquisition MAT, metadata JSON, target CSV and image/mask PNG.
E.PlotProducts={'camera','detector','radec','stages','dataset'};
E.PlotPhotometry=false; % Optional diffuse/specular and received-irradiance histories.
E.PlotExposures=1; % First exposure; use 'all' if NumExposures is increased.

%% APPLY AND RUN (normally leave this section unchanged)
C=SPOT_ApplyExample(C,E);
Out=RUN_SPOT(C);
