function [L,records]=labels(C,maps,records)
%LABELS Eqs.139-153: target-only post-optics masks and deterministic overlap.
shape=[C.Camera.Ny C.Camera.Nx]; L.Binary=false(shape); L.Instance=zeros(shape,'uint32');
L.Semantic=zeros(shape,'uint8'); L.ObjectUnion=false(shape); L.TrailUnion=false(shape);
L.PerTarget=cell(1,numel(maps)); best=zeros(shape); ids=[C.Targets.ID];
order=1:numel(ids); % smallest target array index wins exact ties (paper Eq.149).
for j=order
    m=maps{j}; peak=max(m,[],'all'); kind=records(j).Appearance;
    if strcmp(kind,'object'), threshold=C.Masks.RelativeThreshold; cls=uint8(1);
    else, threshold=C.Masks.RelativeThreshold; cls=uint8(2); end
    mask=peak>0 & m>=threshold*peak; L.PerTarget{j}=mask;
    win=mask & (L.Instance==0 | m>best); L.Instance(win)=uint32(ids(j)); L.Semantic(win)=cls; best(win)=m(win);
    L.Binary=L.Binary|mask;
    if cls==1, L.ObjectUnion=L.ObjectUnion|mask; else, L.TrailUnion=L.TrailUnion|mask; end
    [y,x]=find(mask); records(j).HasAnnotation=~isempty(x); records(j).PeakOpticalSignal=peak;
    records(j).MaskPixelCount=numel(x); records(j).MaskThreshold=threshold;
    if isempty(x)
        records(j).BoundingBox_xywh=[NaN NaN 0 0]; records(j).BoxCenter_xy=[NaN NaN];
    else
        box=[min(x) min(y) max(x)-min(x)+1 max(y)-min(y)+1];
        records(j).BoundingBox_xywh=box; records(j).BoxCenter_xy=box(1:2)+(box(3:4)-1)/2;
    end
end
L.Schema=struct('Background',0,'Object',1,'Trail',2,'InstanceMeaning','Target ID', ...
    'PixelCoordinates','one-based pixel centres; boxes [xmin ymin inclusive_width inclusive_height]', ...
    'OverlapRule','largest target-only optical signal among candidate masks; ties smallest target array index');
end
