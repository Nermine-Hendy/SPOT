function Out=RUN_SPOT(config,quality)
%RUN_SPOT Simulate with independent controls for displaying and saving.
% Out=RUN_SPOT('ground');                 % full-quality settings
% C=SPOT_Config('space_limb','full'); Out=RUN_SPOT(C);
% C.Output.SaveResults controls MAT/CSV/PNG/JSON result export.
% C.Output.Plot controls figure display; C.Plot.Save controls figure export.
SPOT_Setup;
if nargin==0
    C=SPOT_Config;
elseif isstruct(config)
    assert(nargin==1,'SPOT:Run','Pass quality to SPOT_Config before supplying a structure.');
    C=config;
elseif nargin<2
    C=SPOT_Config(config);
else
    C=SPOT_Config(config,quality);
end
if ~isfield(C.Output,'SaveResults'), C.Output.SaveResults=true; end
if ~isfield(C.Plot,'Exposures'), C.Plot.Exposures=1; end
if ~isfield(C.Plot,'Photometry'), C.Plot.Photometry=false; end
for value={C.Output.SaveResults,C.Output.Plot,C.Plot.Save,C.Plot.Photometry}
    validateattributes(value{1},{'logical','numeric'},{'scalar','binary'});
end
indices=C.Plot.Exposures;
if ischar(indices) || isstring(indices)
    assert(strcmp(indices,'all'),'SPOT:Plot','PlotExposures must be an index vector or ''all''.');
    indices=1:C.Time.NumExposures;
else
    validateattributes(indices,{'numeric'},{'vector','integer','>=',1,'<=',C.Time.NumExposures});
    indices=unique(indices,'stable');
end
if C.Output.Plot, C.Plot.Visible='on'; else, C.Plot.Visible='off'; end
% Reserve a fresh sample path only when files will actually be written.
if C.Output.SaveResults || C.Plot.Save
    base=C.Output.Directory; repeat=1;
    while isfolder(fullfile(C.Output.Directory,C.Dataset.SampleID))
        C.Output.Directory=fullfile(base,sprintf('repeat_%03d',repeat));
        repeat=repeat+1;
    end
end
Out=SPOT.run(C,C.Output.SaveResults);
Out.PlotFiles={}; Out.FigureHandles=gobjects(0);
if C.Output.Plot || C.Plot.Save
    for k=indices
        [figs,files]=SPOT_Plots(Out,k);
        Out.PlotFiles=[Out.PlotFiles files];
        if C.Plot.Photometry
            [f,files]=SPOT_PlotPhotometry(Out,k);
            figs=[figs f]; Out.PlotFiles=[Out.PlotFiles files];
        end
        if C.Output.Plot
            Out.FigureHandles=[Out.FigureHandles reshape(figs,1,[])];
        else
            close(figs);
        end
    end
end
if C.Output.SaveResults || C.Plot.Save
    fprintf('SPOT output: %s\n',fullfile(C.Output.Directory,C.Dataset.SampleID));
else
    fprintf('SPOT acquisition returned in Out; saving is disabled.\n');
end
if nargout==0, assignin('base','Out',Out); clear Out; end
end
