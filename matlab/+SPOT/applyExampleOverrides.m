function C=applyExampleOverrides(C)
%APPLYEXAMPLEOVERRIDES Apply explicit user edits after sensor/target CSV rows.
% Empty fields preserve configured settings. No appearance is forced.
if isfield(C,'ExamplePointing')
    names=fieldnames(C.ExamplePointing);
    for k=1:numel(names)
        value=C.ExamplePointing.(names{k});
        if ~isempty(value), C.Pointing.(names{k})=value; end
    end
end
if ~isfield(C,'ExampleTarget'), return; end
T=C.ExampleTarget; names=fieldnames(T);
allowed={'AngleToLOS_deg','RelativeSpeed_mps','Direction_deg','InitialAngles_deg','Range_m'};
assert(all(ismember(names,allowed)),'SPOT:Example','Unsupported target override.');
active=names(~cellfun(@(name) isempty(T.(name)),names));
if isempty(active), return; end
n=numel(C.Targets);
assert(n>0,'SPOT:Example','Nonempty target overrides require at least one target.');
for k=1:numel(active)
    name=active{k}; value=T.(name);
    if strcmp(name,'InitialAngles_deg')
        assert(isnumeric(value) && all(isfinite(value(:))) && ...
            (isequal(size(value),[1 2]) || isequal(size(value),[n 2])), ...
            'SPOT:Example','InitialAngles_deg must be one [angle1 angle2] row or one row per selected target.');
        if size(value,1)==1, value=repmat(value,n,1); end
    else
        assert(isnumeric(value) && all(isfinite(value(:))) && ...
            (isscalar(value) || (isvector(value) && numel(value)==n)), ...
            'SPOT:Example','%s must be a scalar or have one entry per selected target.',name);
        if isscalar(value), value=repmat(value,1,n); end
        value=reshape(value,1,[]);
        if strcmp(name,'AngleToLOS_deg')
            assert(all(value>=0 & value<=180),'SPOT:Example','AngleToLOS_deg must be in [0,180].');
        elseif strcmp(name,'RelativeSpeed_mps')
            assert(all(value>=0),'SPOT:Example','RelativeSpeed_mps must be nonnegative.');
        elseif strcmp(name,'Range_m')
            assert(all(value>0),'SPOT:Example','Range_m must be positive.');
        end
    end
    T.(name)=value;
end
for j=1:n
    assert(strcmp(C.Targets(j).Trajectory.Type,'seeded_state'), ...
        'SPOT:Example','Motion/angle overrides require seeded targets. Target %d uses %s. Leave overrides empty for orbital/ephemeris/TLE inputs.', ...
        j,C.Targets(j).Trajectory.Type);
    d=C.Targets(j).Seed;
    if ismember('AngleToLOS_deg',active) || ismember('RelativeSpeed_mps',active)
        speed=hypot(d.TransverseSpeed_mps,d.RadialSpeed_mps);
        angle=atan2d(d.TransverseSpeed_mps,d.RadialSpeed_mps);
        if ismember('RelativeSpeed_mps',active), speed=T.RelativeSpeed_mps(j); end
        if ismember('AngleToLOS_deg',active), angle=T.AngleToLOS_deg(j); end
        d.TransverseSpeed_mps=speed*sind(angle);
        d.RadialSpeed_mps=speed*cosd(angle);
    end
    if ismember('Direction_deg',active), d.Direction_deg=T.Direction_deg(j); end
    if ismember('InitialAngles_deg',active), d.Angles_deg=T.InitialAngles_deg(j,:); end
    if ismember('Range_m',active), d.Range_m=T.Range_m(j); end
    C.Targets(j).Seed=d;
end
end
