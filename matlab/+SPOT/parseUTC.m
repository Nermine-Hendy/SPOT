function t=parseUTC(value)
%PARSEUTC Parse ISO calendar timestamps without locale or format inference.
% Accept Z or +/-HH:MM / +/-HHMM, with optional 1-9 fractional-second digits.
% A timezone-aware datetime is accepted directly; unzoned values are rejected.
if isdatetime(value)
    assert(~isempty(value.TimeZone) && all(~isnat(value(:))),'SPOT:UTC','Datetime input must have a timezone and contain no NaT.');
    t=value(:); t.TimeZone='UTC'; return;
end
s=string(value); s=strtrim(s(:)); t=NaT(numel(s),1,'TimeZone','UTC');
pattern=['^(?<year>\d{4})-(?<month>\d{2})-(?<day>\d{2})[Tt ]' ...
    '(?<hour>\d{2}):(?<minute>\d{2}):(?<second>\d{2}(?:\.\d{1,9})?)' ...
    '(?<zone>[Zz]|[+-]\d{2}:?\d{2})$'];
for k=1:numel(s)
    assert(~ismissing(s(k)),'SPOT:UTC','Missing timestamp at row %d.',k);
    p=regexp(char(s(k)),pattern,'names','once');
    assert(~isempty(p),'SPOT:UTC','Timestamp row %d must use yyyy-MM-ddTHH:mm:ss[.fraction]Z or an explicit numeric UTC offset.',k);
    yyyy=str2double(p.year); mm=str2double(p.month); dd=str2double(p.day);
    hh=str2double(p.hour); mi=str2double(p.minute); ss=str2double(p.second);
    assert(yyyy>=1 && mm>=1 && mm<=12 && dd>=1 && dd<=31 && hh<24 && mi<60 && ss<60, ...
        'SPOT:UTC','Invalid calendar/time fields at row %d; leap-second text is not supported by this UTC model.',k);
    whole=floor(ss);
    try
        local=datetime(yyyy,mm,dd,hh,mi,whole,'TimeZone','UTC');
    catch
        error('SPOT:UTC','Invalid calendar date at row %d.',k);
    end
    assert(year(local)==yyyy && month(local)==mm && day(local)==dd,'SPOT:UTC','Invalid calendar date at row %d.',k);
    offset=0;
    if ~strcmpi(p.zone,'Z')
        zone=erase(p.zone,':'); oh=str2double(zone(2:3)); om=str2double(zone(4:5));
        assert(oh<24 && om<60,'SPOT:UTC','Invalid UTC offset at row %d.',k);
        direction=1; if zone(1)=='-', direction=-1; end
        offset=direction*(60*oh+om);
    end
    t(k)=local+seconds(ss-whole)-minutes(offset);
end
end
