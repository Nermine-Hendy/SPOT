function SPOT_PlotDataset(f,C,E)
tl=tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); SPOT_PlotDetector(ax,C,E,true); title(ax,'Detector image');
ax=nexttile(tl); plotMask(ax,E.Labels.Binary,[0 0 0;1 1 1],C); clim(ax,[0 1]); title(ax,'Binary target mask');
% Display categorical IDs with compact color indices; stored IDs stay intact.
im=zeros(size(E.Labels.Instance)); colors=zeros(numel(E.Targets)+1,3);
for j=1:numel(E.Targets)
    im(E.Labels.Instance==E.Targets(j).ID)=j;
    colors(j+1,:)=C.Plot.Colors(1+mod(j-1,size(C.Plot.Colors,1)),:);
end
if size(colors,1)==1, colors=[colors;1 1 1]; end
ax=nexttile(tl); plotMask(ax,im,colors,C); clim(ax,[0 max(1,numel(E.Targets))]); title(ax,'Instance mask');
ax=nexttile(tl); plotMask(ax,E.Labels.Semantic,[0 0 0;.55 .55 .55;1 1 1],C); clim(ax,[0 2]); title(ax,'Semantic mask: compact / trail');
ax=nexttile(tl); SPOT_PlotDetector(ax,C,E,false); plotTruth(ax,C,E,'detector',true,0); title(ax,'Detector, target truth and bounding boxes');
ax=nexttile(tl); axis(ax,'off'); R=E.Targets;
lines={sprintf('Exposure: %d   duration: %.4g s',E.Index,C.Time.Exposure_s), ...
    ['Scene: ' E.Truth.GeometryRegime],sprintf('Image: %d x %d px',C.Camera.Nx,C.Camera.Ny), ...
    sprintf('Annotated targets: %d / %d',nnz([R.HasAnnotation]),numel(R)), ...
    sprintf('Compact / trail pixels: %d / %d',nnz(E.Labels.Semantic==1),nnz(E.Labels.Semantic==2)), ...
    ' ', 'ID  appearance  smear / threshold [px]'};
for j=1:numel(R)
    m=R(j).Morphology;
    lines{end+1}=sprintf('%d   %-10s  %.3g / %.3g',R(j).ID,R(j).Appearance,m.EquivalentSmear_pix,m.Threshold_pix); %#ok<AGROW>
end
lines{end+1}=' '; lines{end+1}='One exposure-integrated renderer for all targets.';
text(ax,0,1,strjoin(lines,newline),'Units','normalized','VerticalAlignment','top', ...
    'FontName','FixedWidth','FontSize',max(8,C.Plot.FontSize-1),'Interpreter','none');
title(ax,'Label and metadata summary');
end
