function plotTruth(ax,C,E,mode,boxes,ra0)
hold(ax,'on');
for j=1:numel(E.Targets)
    R=E.Targets(j); H=R.Samples; ok=H.Visible & H.RadiometricRate>0;
    if ~any(ok), continue; end
    switch mode
        case 'camera', x=H.Xi_deg; y=H.Eta_deg;
        case 'radec', x=wrapAngle(H.RA_deg-ra0); y=H.DEC_deg;
        otherwise, x=H.x; y=H.y;
    end
    color=C.Plot.Colors(1+mod(j-1,size(C.Plot.Colors,1)),:);
    x(~ok)=NaN; y(~ok)=NaN; label=find(ok,1,'last');
    lx=x(label); ly=y(label);
    if strcmp(mode,'radec')
        seam=[false abs(diff(x))>180]; x(seam)=NaN; y(seam)=NaN;
    end
    plot(ax,x,y,'-','Color',color,'LineWidth',C.Plot.TruthLineWidth);
    if strcmp(R.Appearance,'object')
        % A small plotting symbol marks a compact integrated trace, not its size.
        plot(ax,lx,ly,'o','Color',color,'MarkerFaceColor',color,'MarkerSize',C.Plot.MarkerSize);
    end
    text(ax,lx,ly,sprintf(' T%d',R.ID),'Color',color,'FontWeight','bold', ...
        'FontSize',C.Plot.TargetFontSize,'Clipping','on');
    if boxes && R.HasAnnotation
        b=R.BoundingBox_xywh;
        rectangle(ax,'Position',[b(1:2)-.5 b(3:4)],'EdgeColor',color,'LineWidth',1.2);
    end
end
end
