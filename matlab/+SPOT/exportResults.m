function folder=exportResults(Out,C)
%SPOT_EXPORT Raw MAT, display PNG, exact masks, CSV records, JSON metadata.
% PNG images are quicklooks; physical DN and stage arrays remain in acquisition.mat.
if nargin<2, C=Out.Configuration; end
id=char(C.Dataset.SampleID);
assert(~isempty(regexp(id,'^[A-Za-z0-9_-]+$','once')),'SPOT:Export','SampleID must contain only letters, digits, underscores and hyphens.');
folder=fullfile(C.Output.Directory,id);
if isfolder(folder) && isfile(fullfile(folder,'metadata.json'))
    error('SPOT:ExportExists','Sample already exists: %s. Choose a new SampleID or deliberately remove the old exported sample.',folder);
end
if ~isfolder(folder), mkdir(folder); end
if C.Output.SaveMAT, save(fullfile(folder,'acquisition.mat'),'Out','-v7.3'); end
if C.Output.SaveCSV && isfield(C.Input,'SelectedTargets')
    writetable(struct2table(C.Input.SelectedTargets),fullfile(folder,'selected_targets.csv'));
end
if C.Output.SaveCSV && isfield(C.Input,'SelectedSensor')
    writetable(struct2table(C.Input.SelectedSensor),fullfile(folder,'selected_sensor.csv'));
end
M=struct('SchemaVersion',C.SchemaVersion,'Configuration',Out.Configuration,'Provenance',Out.Provenance, ...
    'ImageUnits','DN; PNG is a separately scaled display preview', ...
    'Coordinates','one-based detector pixel centres; inclusive xywh boxes', ...
    'SemanticClasses',struct('Background',0,'Object',1,'Trail',2),'Exposures',{{}});
for ie=1:numel(Out.Exposures)
    E=Out.Exposures{ie}; prefix=sprintf('exposure_%04d',ie);
    records=E.Targets;
    if C.Output.SaveCSV, writetable(targetTable(records),fullfile(folder,[prefix '_targets.csv'])); end
    annotations=records([records.HasAnnotation]);
    if ~isempty(annotations), annotations=rmfield(annotations,'Samples'); end
    m=struct('Index',ie,'StartUTC',E.StartUTC,'MidUTC',E.MidUTC,'EndUTC',E.EndUTC, ...
        'GeometryRegime',E.Truth.GeometryRegime,'ConservativeRegime',E.Truth.ConservativeRegime, ...
        'EarthFraction',E.Truth.EarthFraction,'Sampling',E.Sampling,'Annotations',annotations, ...
        'Stars',E.StarMetadata,'PNGScale_DN',[]);
    if C.Output.SavePNG
        [im,lim]=quicklook(E.Image_DN,C.Output.PNGPercentiles); m.PNGScale_DN=lim;
        imwrite(im,fullfile(folder,[prefix '_image.png']));
        imwrite(uint8(E.Labels.Binary)*255,fullfile(folder,[prefix '_binary.png']));
        imwrite(E.Labels.Semantic,fullfile(folder,[prefix '_semantic.png']));
        assert(max(E.Labels.Instance,[],'all')<=65535,'SPOT:Export','PNG instance IDs exceed uint16. Disable PNG and use exact uint32 MAT masks.');
        imwrite(uint16(E.Labels.Instance),fullfile(folder,[prefix '_instance.png']));
        imwrite(uint8(E.DetectorStages.ValidPixel)*255,fullfile(folder,[prefix '_valid_pixels.png']));
    end
    M.Exposures{ie}=m;
end
fid=fopen(fullfile(folder,'metadata.json'),'w'); assert(fid>=0,'SPOT:Export','Cannot create metadata.json.');
cleanup=onCleanup(@() fclose(fid)); fprintf(fid,'%s',jsonencode(SPOT.serializable(M,true),'PrettyPrint',true));
clear cleanup
end
function T=targetTable(R)
% Stable scalar-column schema, including scenes with zero targets.
n=numel(R); B=zeros(n,4);
if n>0, B=vertcat(R.BoundingBox_xywh); end
smear=zeros(n,1); psf=smear; threshold=smear; path=smear;
for j=1:n, m=R(j).Morphology; smear(j)=m.EquivalentSmear_pix; psf(j)=m.PSFFWHM_pix; threshold(j)=m.Threshold_pix; path(j)=m.PathLength_pix; end
column=@(field) reshape([R.(field)],n,1);
T=table(column('ID'),reshape(string({R.Name}),n,1),reshape(string({R.Class}),n,1), ...
    reshape(string({R.Appearance}),n,1), ...
    logical(column('AnyVisible')),column('VisibleFraction'),column('IlluminatedFraction'), ...
    column('MeanVisibleRate'),column('MeanReceivedIrradiance_W_m2'),column('MeanApparentMagnitude'),column('IntegratedOpticalSignal'),logical(column('HasAnnotation')), ...
    B(:,1),B(:,2),B(:,3),B(:,4),smear,psf,threshold,path, ...
    'VariableNames',{'target_id','name','object_class','appearance', ...
    'any_visible','visible_fraction','illuminated_fraction','mean_visible_rate','mean_received_irradiance_W_m2','mean_apparent_magnitude','integrated_optical_signal', ...
    'has_annotation','xmin','ymin','width','height','equivalent_smear_pix','psf_fwhm_pix','compact_threshold_pix','visible_path_length_pix'});
end
function [im,lim]=quicklook(a,p)
x=sort(a(isfinite(a))); assert(~isempty(x),'SPOT:Export','Image has no finite pixels.');
q=1+(numel(x)-1)*p/100; lim=interp1(1:numel(x),x,q,'linear');
if lim(2)<=lim(1), lim(2)=lim(1)+1; end
im=uint8(round(255*max(0,min(1,(a-lim(1))/(lim(2)-lim(1))))));
end
