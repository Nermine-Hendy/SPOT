function T = prepareTrajectory(C,T)
%PREPARETRAJECTORY Load an ephemeris or TLE once, preserving the input schema.
T.Type=lower(char(T.Type)); T.Frame=upper(char(T.Frame));
switch T.Type
    case 'ephemeris'
        if ~isempty(T.File)
            options=detectImportOptions(T.File,'VariableNamingRule','preserve');
            if ~ismember('time_utc',options.VariableNames), error('SPOT:Ephemeris','Missing time_utc column.'); end
            % Preserve original ISO text: automatic datetime inference can lose
            % fractional seconds or choose an incompatible display format.
            options=setvartype(options,'time_utc','string');
            A=readtable(T.File,options);
            T.TimeUTC=SPOT.parseUTC(A.time_utc);
            suffix=lower(T.Frame);
            if strcmp(suffix,'eci'), suffix='gcrf'; T.Frame='GCRF'; end
            keys={['x_' suffix '_m'],['y_' suffix '_m'],['z_' suffix '_m']};
            if ~all(ismember(keys,A.Properties.VariableNames)), error('SPOT:Ephemeris','Position columns must match the declared Frame.'); end
            T.Positions_m=double(A{:,keys}).';
            keys={['vx_' suffix '_mps'],['vy_' suffix '_mps'],['vz_' suffix '_mps']};
            if all(ismember(keys,A.Properties.VariableNames)), T.Velocities_mps=double(A{:,keys}).';
            elseif any(ismember(keys,A.Properties.VariableNames)), error('SPOT:Ephemeris','Supply all three velocity columns or none.');
            else, T.Velocities_mps=zeros(3,0);
            end
        end
        T.TimeUTC=T.TimeUTC(:);
        assert(isdatetime(T.TimeUTC) && ~isempty(T.TimeUTC.TimeZone) && all(~isnat(T.TimeUTC)),'SPOT:Ephemeris','Ephemeris TimeUTC must contain timezone-aware finite times.');
        T.TimeUTC.TimeZone='UTC';
        if numel(T.TimeUTC)<2 || any(seconds(diff(T.TimeUTC))<=0), error('SPOT:Ephemeris','Ephemeris needs at least two strictly increasing UTC times.'); end
        validateattributes(T.Positions_m,{'numeric'},{'size',[3 numel(T.TimeUTC)],'finite'});
        if ~isempty(T.Velocities_mps), validateattributes(T.Velocities_mps,{'numeric'},{'size',size(T.Positions_m),'finite'}); end
        if ~ismember(T.Frame,{'GCRF','ECEF'}), error('SPOT:Frame','Ephemeris frame must be GCRF or ECEF; TEME is not interchangeable with GCRF.'); end
    case 'tle'
        if exist('satelliteScenario','file')~=2, error('SPOT:Dependency','TLE input requires MATLAB satelliteScenario with SGP4/SDP4 support.'); end
        if ~isempty(T.TLELines), lines=string(T.TLELines(:)); else, lines=readlines(T.File); end
        lines=strtrim(lines); lines=lines(strlength(lines)>0);
        idx=find(startsWith(lines,'1 '));
        if isempty(idx) || T.TLEIndex>numel(idx), error('SPOT:TLE','Requested TLE index is absent.'); end
        p=idx(T.TLEIndex);
        if p==numel(lines) || ~startsWith(lines(p+1),'2 '), error('SPOT:TLE','TLE line 1 must be followed by line 2.'); end
        if any(strlength(lines(p:p+1))~=69), error('SPOT:TLE','TLE lines must each contain 69 columns, including checksum.'); end
        pair=char(lines(p:p+1));
        for j=1:2
            s=pair(j,:); d=s(1:68); digits=d>='0' & d<='9';
            chk=mod(sum(double(d(digits))-48)+sum(d=='-'),10);
            if chk~=str2double(s(69)), error('SPOT:TLE','Invalid checksum on TLE line %d.',j); end
        end
        if ~strcmp(pair(1,3:7),pair(2,3:7)), error('SPOT:TLE','TLE satellite numbers disagree.'); end
        yy=str2double(pair(1,19:20)); yy=yy+1900+100*(yy<57);
        epoch=datetime(yy,1,1,'TimeZone','UTC')+days(str2double(pair(1,21:32))-1);
        assert(~isnat(epoch),'SPOT:TLE','TLE epoch is invalid.');
        finish=C.Time.StartUTC+seconds((C.Time.NumExposures-1)*C.Time.Cadence_s+C.Time.Exposure_s);
        age=max(abs(days([C.Time.StartUTC finish]-epoch)));
        if age>T.MaxTLEAge_days, error('SPOT:TLEAge','TLE age %.2f days exceeds configured maximum %.2f.',age,T.MaxTLEAge_days); end
        tmp=[tempname '.tle']; cleanup=onCleanup(@() delete(tmp));
        fid=fopen(tmp,'w'); if fid<0, error('SPOT:TLE','Cannot create temporary TLE file.'); end
        fprintf(fid,'%s\n%s\n',pair(1,:),pair(2,:)); fclose(fid);
        sc=satelliteScenario(C.Time.StartUTC,finish,max(.001,min(C.Integration.Step_s,C.Time.Exposure_s)));
        sat=satellite(sc,tmp); % Genuine SGP4 or SDP4 selected from TLE period.
        if ~ismember(lower(string(sat.OrbitPropagator)),["sgp4","sdp4"]), error('SPOT:TLE','TLE must use SGP4/SDP4.'); end
        T.Runtime=struct('Scenario',sc,'Satellite',sat);
        T.TLELines=cellstr(lines(p:p+1)); T.EpochUTC=epoch; T.Frame='GCRF';
        clear cleanup
    case {'geodetic','ecef','state','keplerian','demo'}
        % Analytic/state inputs require no scenario objects.
    otherwise, error('SPOT:Trajectory','Unsupported trajectory type %s.',T.Type);
end
end
