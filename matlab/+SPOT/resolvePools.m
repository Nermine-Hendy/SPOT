function C=resolvePools(C)
%RESOLVEPOOLS Convert selected pool records into the SAME editable configuration.
% Pools are opt-in. Selection is recorded and never perturbs detector RNG.
if C.Input.UseSensorCSV
    T=readtable(C.Input.SensorCSV,'TextType','string'); uniqueIDs(T,'sensor_id');
    [q,index]=SPOT.selectPoolRows(T,C.Input.SensorSelection,1,'sensor_id');
    tr=C.Observer.Trajectory; tr.Type=char(q.position_type); tr.EpochUTC=C.Time.StartUTC;
    switch tr.Type
        case 'geodetic', tr.Frame='ECEF'; tr.Motion='fixed'; tr.Geodetic_deg_m=[q.lat_deg q.lon_deg q.alt_m];
        case 'ecef', tr.Frame='ECEF'; tr.Motion='fixed'; tr.Position_m=[q.x_ecef_m;q.y_ecef_m;q.z_ecef_m];
        case 'keplerian', tr.Frame='GCRF'; tr.Elements=[q.a_m q.e q.i_deg q.raan_deg q.argp_deg q.ta_deg];
        case 'ephemeris_csv'
            tr.Type='ephemeris'; tr.File=char(q.ephemeris_file);
            if ~isfile(tr.File), tr.File=fullfile(fileparts(C.Input.SensorCSV),tr.File); end
            opts=detectImportOptions(tr.File);
            if ismember('x_gcrf_m',opts.VariableNames), tr.Frame='GCRF';
            elseif ismember('x_ecef_m',opts.VariableNames), tr.Frame='ECEF';
            else, error('SPOT:Pool','Ephemeris requires canonical x_gcrf_m or x_ecef_m columns.'); end
        otherwise, error('SPOT:Pool','Unsupported sensor pool position_type: %s.',tr.Type);
    end
    C.Observer.Trajectory=tr; C.Observer.Name=char(q.name);
    C.Pointing.Mode=char(q.orientation_mode);
    if strcmp(C.Pointing.Mode,'boresight_eci'), C.Pointing.Mode='inertial'; end
    C.Pointing.AzEl_deg=[q.azimuth_deg q.elevation_deg]; C.Pointing.Roll_deg=q.roll_deg;
    C.Pointing.OffNadir_deg=q.offnadir_deg; C.Pointing.SpaceAzimuth_deg=q.space_azimuth_deg;
    C.Pointing.Boresight_GCRF=[q.bx_eci;q.by_eci;q.bz_eci];
    C.Camera.Shape=char(q.fov_shape); C.Camera.FOV_deg=[q.fov_azimuth_deg q.fov_elevation_deg];
    if strcmp(C.Camera.Shape,'circular'), C.Camera.FOV_deg=q.fov_diameter_deg*[1 1]; end
    C.Input.SelectedSensorRow=index; C.Input.SelectedSensor=table2struct(q);
end
if C.Input.UseTargetCSV
    T=readtable(C.Input.TargetCSV,'TextType','string'); uniqueIDs(T,'target_id');
    [rows,index]=SPOT.selectPoolRows(T,C.Input.TargetSelection,C.Input.NumTargets,'target_id');
    assert(~isempty(C.Targets) || isempty(index),'SPOT:Pool','Keep one target template when enabling a nonempty pool.');
    if ~isempty(index), template=C.Targets(1); C.Targets=repmat(template,1,numel(index)); else, C.Targets=C.Targets([]); end
    for j=1:numel(index)
        q=rows(j,:);
        if ismember('photometry_model',q.Properties.VariableNames)
            assert(strcmp(q.photometry_model,'surface_brdf'),'SPOT:Pool','Target pool photometry_model must be surface_brdf.');
        end
        tg=template; tg.ID=index(j); tg.Name=char(q.target_id); tg.Class=char(q.target_class);
        tr=tg.Trajectory; tr.EpochUTC=C.Time.StartUTC;
        switch char(q.source)
            case 'seeded'
                tr.Type='seeded_state'; d=tg.Seed;
                if strcmp(C.Pointing.Mode,'az_el_roll')
                    d.AngleMode='az_el_offset'; d.Angles_deg=[q.az_offset_deg q.el_offset_deg]; d.Range_m=q.range_m;
                else
                    d.AngleMode='theta_phi'; d.Angles_deg=[q.space_theta_deg q.space_phi_deg]; d.Range_m=q.space_range_m;
                end
                d.TransverseSpeed_mps=q.relative_speed_mps*sind(q.angle_to_los_deg);
                d.RadialSpeed_mps=q.relative_speed_mps*cosd(q.angle_to_los_deg);
                d.Direction_deg=q.position_angle_deg; tg.Seed=d;
            case 'keplerian'
                tr.Type='keplerian'; tr.Frame='GCRF'; tr.Elements=[q.a_m q.e q.i_deg q.raan_deg q.argp_deg q.ta_deg];
            otherwise, error('SPOT:Pool','Unsupported target pool source. Use direct Trajectory configuration for other inputs.');
        end
        tg.Trajectory=tr; tg.Photometry.Model='surface_brdf';
        % CSV area is TOTAL area of a convex six-face cube in this schema.
        tg.Surface.NormalBody=[eye(3) -eye(3)]; tg.Surface.Area_m2=q.area_m2/6*ones(1,6);
        tg.Surface.DiffuseReflectance=q.albedo*ones(1,6);
        tg.Surface.SpecularReflectance=value(q,'specular_reflectance',0.1)*ones(1,6);
        tg.Surface.SpecularExponent=value(q,'specular_exponent',40)*ones(1,6);
        tg.Attitude.Mode='spin'; tg.Attitude.EpochUTC=C.Time.StartUTC;
        tg.Attitude.SpinRate_deg_s=value(q,'spin_rate_deg_s',5);
        tg.Attitude.SpinAxisBody=[value(q,'spin_axis_x',0);value(q,'spin_axis_y',0);value(q,'spin_axis_z',1)];
        angles=[value(q,'body_roll_deg',0) value(q,'body_pitch_deg',0) value(q,'body_yaw_deg',0)];
        tg.Attitude.BodyToGCRF0=SPOT.eulerBodyRotation(angles);
        tg.Metadata.PoolID=char(q.target_id); tg.Metadata.PoolRow=index(j); tg.Metadata.PoolRecord=table2struct(q);
        tg.Metadata.AreaConvention='CSV area_m2 = total surface area of convex cube; per-face area = area_m2/6.';
        C.Targets(j)=tg;
    end
    C.Input.SelectedTargetRows=index; C.Input.SelectedTargets=table2struct(rows);
end
C.Input.Resolved=true;
end
function z=value(q,name,fallback)
z=fallback;
if ismember(name,q.Properties.VariableNames)
    z=q.(name); assert(isnumeric(z) && isscalar(z) && isfinite(z),'SPOT:Pool','Invalid %s.',name);
end
end
function uniqueIDs(T,name)
assert(ismember(name,T.Properties.VariableNames),'SPOT:Pool','Missing pool ID column.');
a=string(T.(name)); assert(all(~ismissing(a) & strlength(a)>0) && numel(unique(a))==height(T),'SPOT:Pool','Pool IDs must be nonempty and unique.');
end
