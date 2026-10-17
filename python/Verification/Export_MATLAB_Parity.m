function Export_MATLAB_Parity(matlabRoot, outputDirectory, frameModel)
% Export deterministic fixtures from the MATLAB SPOT implementation.
% Example:
% Export_MATLAB_Parity('C:\path\to\SPOT\matlab', ...
%                     'C:\path\to\SPOT\python\Verification\matlab_fixtures')
% Requires the MATLAB simulator and its toolbox dependencies.
if nargin<2, outputDirectory=fullfile(fileparts(mfilename('fullpath')),'matlab_fixtures'); end
if nargin<3, frameModel='iau2006'; end
assert(isfolder(fullfile(matlabRoot,'+SPOT')),'Expected the MATLAB folder containing +SPOT.');
addpath(matlabRoot); SPOT_Setup;
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
for item={'ground','space_sky','space_limb','space_earth','figure8','multi_exposure','current_ground'}
    C=SPOT_Config(item{1},'full');
    % Deliberate parity fixture changes; public notebook presets remain full.
    C.Camera.Nx=128; C.Camera.Ny=128;
    C.Frames.Model=frameModel;
    C.Atmosphere.Scintillation.Enabled=false;
    C.Detector.PoissonEnabled=false; C.Detector.ReadEnabled=false;
    C.Output.SaveResults=false; C.Output.Plot=false; C.Plot.Save=false;
    config_json=jsonencode(SPOT.serializable(C,true));
    Out=SPOT.run(C,false);
    for k=1:numel(Out.Exposures)
        E=Out.Exposures{k};
        F=struct('Image_DN',E.Image_DN,'SceneRate',E.SceneRate,'TargetRate',E.TargetRate, ...
            'StarRate',E.StarRate,'SkyRate',E.SkyRate,'EarthRate',E.EarthRate, ...
            'OpticalSignal',E.OpticalStages.PSF,'ExpectedSignal_e',E.DetectorStages.ExpectedSignal_e, ...
            'Binary',E.Labels.Binary,'Instance',E.Labels.Instance,'Semantic',E.Labels.Semantic, ...
            'Xi_deg',E.Truth.Xi_deg,'Eta_deg',E.Truth.Eta_deg,'RA_deg',E.Truth.RA_deg, ...
            'DEC_deg',E.Truth.DEC_deg,'EarthMask',E.Truth.EarthMask);
        exposure_index=k;
        filename=fullfile(outputDirectory,sprintf('%s_%02d.mat',item{1},k));
        save(filename,'F','config_json','exposure_index','-v7');
        fprintf('Exported %s\n',filename);
    end
end
end
