function SPOT_PlotDetector(ax,C,E,overlay,showColorbar)
if nargin<5, showColorbar=false; end
[shown,limits]=displayImage(E.Image_DN,C);
if showColorbar
    imagesc(ax,E.Image_DN); clim(ax,limits);
    cb=colorbar(ax); cb.Label.String='Digital number (DN)';
    title(ax,sprintf('Detector image plane | %s | exposure %d',E.Truth.GeometryRegime,E.Index),'Interpreter','none');
else
    imagesc(ax,shown); clim(ax,[0 1]);
end
colormap(ax,gray(256));
axis(ax,'image'); set(ax,'YDir','reverse'); hold(ax,'on');
if overlay, plotTruth(ax,C,E,'detector',false,0); end
xlabel(ax,'x pixel'); ylabel(ax,'y pixel'); styleAxes(ax,C);
end
