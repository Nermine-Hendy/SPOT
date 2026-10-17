function [f,files]=SPOT_PlotPhotometry(Out,k)
%SPOT_PLOTPHOTOMETRY Recorded diffuse/specular irradiance and received brightness.
% Supplemental diagnostic: uses actual per-target acquisition histories.
if nargin<2, k=1; end
C=Out.Configuration; E=Out.Exposures{k};
f=figure('Color','w','Visible',C.Plot.Visible,'Name','SPOT target photometry', ...
    'Position',[50 50 1200 800],'Tag','SPOTVisualization');
tl=tiledlayout(f,2,1,'TileSpacing','compact'); a=nexttile(tl); hold(a,'on'); b=nexttile(tl); hold(b,'on');
for j=1:numel(E.Targets)
    T=E.Targets(j); H=T.Samples; t=seconds(H.TimeUTC-E.StartUTC); color=C.Plot.Colors(1+mod(j-1,size(C.Plot.Colors,1)),:);
    if isfield(H,'Reflection')
        plot(a,t,H.Reflection.DiffuseIrradiance_W_m2,'Color',color,'LineWidth',1.5,'DisplayName',[T.Name ' diffuse']);
        plot(a,t,H.Reflection.SpecularIrradiance_W_m2,'--','Color',color,'LineWidth',1.5,'DisplayName',[T.Name ' specular']);
    end
    plot(b,t,H.ReceivedIrradiance_W_m2.*H.Visible,'Color',color,'LineWidth',1.5,'DisplayName',T.Name);
end
xlabel(a,'Time within exposure (s)'); ylabel(a,'Unattenuated irradiance (W m^{-2})'); grid(a,'on'); legend(a,'show','Location','best');
xlabel(b,'Time within exposure (s)'); ylabel(b,'Visible received irradiance (W m^{-2})'); grid(b,'on'); legend(b,'show','Location','best');
title(a,'Surface reflection before eclipse and atmosphere'); title(b,'After eclipse, atmosphere and visibility');
folder=fullfile(C.Output.Directory,C.Dataset.SampleID,'Figures');
files=saveFigure(f,C,folder,sprintf('Photometry_exp%02d',k));
end
