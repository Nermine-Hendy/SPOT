function files=saveFigure(fig,C,folder,stem)
%SAVEFIGURE Save one scientific figure using the unified plotting settings.
files={};
if ~C.Plot.Save, return; end
if ~isfolder(folder), mkdir(folder); end
for j=1:numel(C.Plot.Formats)
    ext=char(C.Plot.Formats{j}); path=fullfile(folder,[stem '.' ext]);
    switch ext
        case 'fig', savefig(fig,path);
        case 'png', exportgraphics(fig,path,'Resolution',C.Plot.DPI,'BackgroundColor','white');
        case 'pdf', exportgraphics(fig,path,'ContentType','vector','BackgroundColor','white');
        otherwise, error('SPOT:Plot','Supported formats: png, fig, pdf.');
    end
    files{end+1}=path; %#ok<AGROW>
end
end
