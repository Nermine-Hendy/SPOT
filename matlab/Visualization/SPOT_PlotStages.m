function SPOT_PlotStages(f,C,E)
D=E.Truth; O=E.OpticalStages; Q=E.DetectorStages;
blocked=D.EarthMask;
if strcmp(D.GeometryRegime,'ground_sky'), blocked=D.Active & D.Elevation_deg<C.Visibility.MinElevation_deg; end
sky=E.SkyRate./max(D.SolidAngle_sr,realmin); earth=E.EarthRate./max(D.SolidAngle_sr,realmin);
bias=0; if C.Detector.BiasEnabled, bias=C.Detector.Bias_DN; end
data={D.Active,blocked,D.SkyMask,D.EarthMask,E.StarRate,sky, ...
    earth,E.TargetRate,E.StarRate+E.TargetRate,sky+earth,E.SceneRate,O.Vignetting, ...
    O.Vignetted,O.PSF,Q.ExpectedSignal_e,Q.ExpectedDark_e,Q.ReadNoise_e,Q.Analog_DN, ...
    E.Image_DN-bias-Q.ExpectedDark_e/C.Detector.ConversionGain_e_per_DN,E.Image_DN, ...
    Q.SignalElectrons,Q.BeforeFullWell_e,Q.FullWellSaturated|Q.ADCSaturated,zeros(size(E.Image_DN))};
titles={'In-FoV mask','Earth / horizon blocked mask','Visible sky mask','Earth-in-FoV mask','Stars flux only','Sky background radiance', ...
    'Earth radiance only','Targets flux only','Point-source flux input','Diffuse radiance input','Radiometric input to optics','Vignetting map', ...
    'After vignetting','After optics / PSF','Electron mean image','Dark current mean','Read noise realization','Pre-ADC signal', ...
    'Bias/dark-subtracted image','Final detector image','After shot noise','After dark current','Saturation mask','Detector truth map'};
tl=tiledlayout(f,4,6,'TileSpacing','compact','Padding','compact');
for j=1:24
    ax=nexttile(tl); a=double(data{j});
    if ismember(j,[5:11 13 14]) && any(a>0,'all'), a=log10(1+max(0,a)/max(max(a,[],'all')*1e-5,realmin)); end
    imagesc(ax,a); axis(ax,'image'); set(ax,'YDir','reverse'); colormap(ax,gray(256));
    title(ax,titles{j},'FontSize',max(7,C.Plot.FontSize-2)); xlabel(ax,'x pixel'); ylabel(ax,'y pixel');
    styleAxes(ax,C); set(ax,'FontSize',max(6,C.Plot.FontSize-3));
    if j==24, clim(ax,[0 1]); plotTruth(ax,C,E,'detector',false,0); end
end
end
