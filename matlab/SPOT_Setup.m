function root=SPOT_Setup
%SPOT_SETUP Add the SPOT root, configuration and visualization interfaces.
% The +SPOT package and Visualization/private resolve automatically.
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'Configuration'), ...
    fullfile(root,'Configuration','internal'),fullfile(root,'Visualization'),'-begin');
expected=fullfile(root,'Configuration','SPOT_Config.m');
assert(strcmpi(which('SPOT_Config'),expected),'SPOT:Path', ...
    'Another SPOT_Config is taking precedence. Change Current Folder to this SPOT root or its Examples folder.');
expected=fullfile(root,'+SPOT','run.m');
assert(strcmpi(which('SPOT.run'),expected),'SPOT:Path', ...
    'Another SPOT package is taking precedence. Remove the old installation from the MATLAB path.');
if nargout==0, clear root; end
end
