function SPOT_PlotCamera(ax,C,E)
D=E.Truth; support=.07*double(D.Active)+.13*double(D.EarthMask);
imagesc(ax,D.Xi_deg(1,:),D.Eta_deg(:,1),support); axis(ax,'image');
set(ax,'YDir','normal'); colormap(ax,gray(256)); clim(ax,[0 1]); hold(ax,'on');
S=E.StarTruth;
scatter(ax,S.Xi_deg,S.Eta_deg,C.Plot.StarPointSize,[.8 .8 .8],'filled','MarkerFaceAlpha',.35);
if strcmp(C.Camera.Shape,'circular')
    th=linspace(0,2*pi,361); q=C.Camera.FOV_deg(1)/2;
    plot(ax,atand(tand(q)*cos(th)),atand(tand(q)*sin(th)),'g-','LineWidth',2);
else
    x=C.Camera.FOV_deg(1)/2; y=C.Camera.FOV_deg(2)/2;
    plot(ax,[-x x x -x -x],[-y -y y y -y],'g-','LineWidth',2);
end
plot(ax,0,0,'wx','MarkerSize',12,'LineWidth',3);
plotTruth(ax,C,E,'camera',false,0);
xlabel(ax,'Cross-boresight camera angle, \xi [deg]');
ylabel(ax,'Along-boresight camera angle, \eta [deg]'); styleAxes(ax,C);
end
