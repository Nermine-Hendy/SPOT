function b=serializable(a,forJSON)
%SERIALIZABLE Drop transient backend handles; optionally encode UTC ISO times.
if nargin<2, forJSON=true; end
if isstruct(a)
    b=a;
    if isfield(b,'Runtime'), b=rmfield(b,'Runtime'); end
    f=fieldnames(b);
    for j=1:numel(b)
        for k=1:numel(f), b(j).(f{k})=SPOT.serializable(b(j).(f{k}),forJSON); end
    end
elseif iscell(a)
    b=cell(size(a)); for k=1:numel(a), b{k}=SPOT.serializable(a{k},forJSON); end
elseif isdatetime(a) && forJSON
    a.TimeZone='UTC'; a.Format="yyyy-MM-dd'T'HH:mm:ss.SSSSSS'Z'"; b=string(a);
else, b=a;
end
end
