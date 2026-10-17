function plotMask(ax,a,map,C)
imagesc(ax,a); axis(ax,'image'); set(ax,'YDir','reverse'); colormap(ax,map);
xlabel(ax,'x pixel'); ylabel(ax,'y pixel'); styleAxes(ax,C);
end
