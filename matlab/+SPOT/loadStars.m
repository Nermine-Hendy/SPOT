function S=loadStars(C)
%LOADSTARS Load one catalogue or reproducible full-sphere synthetic field.
S=struct('ID',strings(0,1),'Direction_GCRF',zeros(3,0),'Magnitude',zeros(0,1), ...
    'Metadata',struct('Enabled',C.Stars.Enabled,'Source',C.Stars.Source,'Count',0));
if ~C.Stars.Enabled, return; end
if strcmp(C.Stars.Source,'csv')
    assert(isfile(C.Stars.File),'SPOT:Catalogue','Catalogue file is missing: %s.',C.Stars.File);
    opts=detectImportOptions(C.Stars.File,'VariableNamingRule','preserve');
    if ismember('source_id',opts.VariableNames), opts=setvartype(opts,'source_id','string'); end
    A=readtable(C.Stars.File,opts); names=A.Properties.VariableNames;
    assert(all(ismember({'ra','dec','magnitude'},names)),'SPOT:Catalogue','Catalogue needs ra, dec, magnitude (degrees, magnitude).');
    ra=double(A.ra); dec=double(A.dec); mag=double(A.magnitude);
    if ismember('source_id',names), id=string(A.source_id); else, id="row_"+string((1:height(A))'); end
    S.Metadata.Catalogue=C.Stars.Catalogue; S.Metadata.Band=C.Stars.Band;
    S.Metadata.EpochUTC=C.Stars.EpochUTC;
    S.Metadata.Astrometry='Static GCRF/ICRS-aligned directions at catalogue epoch; proper motion, parallax and aberration not applied.';
else
    rs=RandStream('mt19937ar','Seed',mod(C.RandomSeed+900001,2^32)); n=C.Stars.SyntheticCount;
    ra=360*rand(rs,n,1); dec=asind(2*rand(rs,n,1)-1);
    limits=C.Stars.SyntheticMagnitudeRange; mag=limits(1)+diff(limits)*rand(rs,n,1);
    id="synthetic_"+string((1:n)'); S.Metadata.Catalogue='Uniform sphere, uniform magnitude distribution';
end
assert(all(isfinite([ra;dec;mag])) && all(abs(dec)<=90),'SPOT:Catalogue','Catalogue contains invalid astrometry or magnitudes.');
use=find(mag<=C.Stars.MaxMagnitude); S.Metadata.CountBeforeCap=numel(use);
[~,order]=sort(mag(use),'ascend'); use=use(order(1:min(numel(order),C.Stars.MaxStars)));
S.ID=id(use); S.Magnitude=mag(use);
S.Direction_GCRF=[(cosd(dec(use)).*cosd(ra(use)))';(cosd(dec(use)).*sind(ra(use)))';sind(dec(use))'];
S.Metadata.Count=numel(use); S.Metadata.Truncated=numel(use)<S.Metadata.CountBeforeCap;
end
