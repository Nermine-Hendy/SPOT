function SPOT_PlotRADEC(ax,C,E)
D=E.Truth; [~,k]=min(abs(seconds(E.TargetTimeUTC-E.MidUTC)));
b=E.CameraToGCRF(:,3,k); ra0=mod(atan2d(b(2),b(1)),360); dec0=asind(b(3));
s=C.Plot.RADECStride; ix=unique([1:s:C.Camera.Nx C.Camera.Nx]); iy=unique([1:s:C.Camera.Ny C.Camera.Ny]);
ra=wrapAngle(D.RA_deg(iy,ix)-ra0); dec=D.DEC_deg(iy,ix);
c=.08+.13*double(D.EarthMask(iy,ix)); c(~D.Active(iy,ix))=NaN;
seam=[false(size(ra,1),1) abs(diff(ra,1,2))>180] | [false(1,size(ra,2));abs(diff(ra,1,1))>180];
ra(seam)=NaN;
surface(ax,ra,dec,zeros(size(ra)),c,'EdgeColor','none','FaceColor','flat');
view(ax,2); colormap(ax,gray(256)); clim(ax,[0 1]); set(ax,'Color','k'); hold(ax,'on');
S=E.StarTruth; scatter(ax,wrapAngle(S.RA_deg-ra0),S.DEC_deg,C.Plot.StarPointSize,[.8 .8 .8],'filled','MarkerFaceAlpha',.35);
plot(ax,0,dec0,'wx','MarkerSize',12,'LineWidth',3); plotTruth(ax,C,E,'radec',false,ra0);
if C.Plot.RADECFullScale, xlim(ax,[-180 180]); ylim(ax,[-90 90]); end
grid(ax,'on'); set(ax,'GridColor',[.65 .65 .65]);
xlabel(ax,sprintf('Right ascension, \\alpha - %.2f^{\\circ} [deg]',ra0));
ylabel(ax,'Declination, \delta [deg]'); styleAxes(ax,C);
end
