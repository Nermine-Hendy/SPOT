function [figs,files]=SPOT_Plots(Out,k,products,stem)
%SPOT_PLOTS Camera, detector, RA/DEC, stage and dataset visualizations.
% The plots read the SAME truth arrays and masks used by the dataset exporter.
% No artificial line lengths, forced object sizes, fixed-radius masks or
% independent plotting configuration are used. Colored marks are annotations.
if nargin<2, k=1; end
C=Out.Configuration; E=Out.Exposures{k};
if nargin<3
    enabled=[C.Plot.Camera C.Plot.Detector C.Plot.RADEC C.Plot.Stages C.Plot.Dataset];
    names={'camera','detector','radec','stages','dataset'}; products=names(enabled);
end
if ischar(products) || isstring(products), products=cellstr(products); end
if nargin<4, stem=sprintf('%s_exp%02d',C.Name,k); end
if isfield(C.Plot,'Exposures'), validateattributes(k,{'numeric'},{'scalar','integer','>=',1,'<=',numel(Out.Exposures)}); end
folder=fullfile(C.Output.Directory,C.Dataset.SampleID,'Figures');
figs=gobjects(1,numel(products)); files={};
for j=1:numel(products)
    kind=products{j}; sz=C.Plot.SingleSize_px;
    if strcmp(kind,'stages'), sz=C.Plot.StageSize_px; end
    if strcmp(kind,'dataset'), sz=C.Plot.DatasetSize_px; end
    f=figure('Name',['SPOT ' kind],'Color','w','Visible',C.Plot.Visible, ...
        'Position',[60 60 sz],'Tag','SPOTVisualization'); figs(j)=f;
    switch kind
        case 'camera', SPOT_PlotCamera(axes('Parent',f),C,E);
        case 'detector'
            showDN=isfield(C.Plot,'DetectorColorbar') && C.Plot.DetectorColorbar;
            SPOT_PlotDetector(axes('Parent',f),C,E,true,showDN);
        case 'radec', SPOT_PlotRADEC(axes('Parent',f),C,E);
        case 'stages', SPOT_PlotStages(f,C,E);
        case 'dataset', SPOT_PlotDataset(f,C,E);
        otherwise, error('SPOT:Plot','Unknown product %s.',kind);
    end
    drawnow;
    tag=[stem '_' kind]; if isscalar(products), tag=stem; end
    files=[files saveFigure(f,C,folder,tag)]; %#ok<AGROW>
    if strcmp(kind,'detector') && C.Plot.Save && C.Plot.ExportRawDetector
        % Annotation-free preview accompanies the annotated scientific figure.
        rawFile=fullfile(folder,[tag '_unannotated.png']);
        imwrite(displayImage(E.Image_DN,C),rawFile); files{end+1}=rawFile;
    end
end
end
