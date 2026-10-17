function C=SPOT_ConfigTrajectoryExample(C,name)
switch name
    case {'trajectory_geodetic','trajectory_ecef'}
        C.Observer.Trajectory.Type='geodetic'; 
        C.Observer.Trajectory.Frame='ECEF';
        C.Targets=C.Targets(1); 
        C.Targets(1).Trajectory=C.Observer.Trajectory;
        C.Targets(1).Trajectory.Geodetic_deg_m=[-37.8 144.96 800e3];
        if strcmp(name,'trajectory_ecef')
            C.Observer.Trajectory.Position_m=SPOT.geodetic(C,C.Observer.Trajectory.Geodetic_deg_m(:),'to_ecef');
            C.Targets(1).Trajectory.Position_m=SPOT.geodetic(C,C.Targets(1).Trajectory.Geodetic_deg_m(:),'to_ecef');
            C.Observer.Trajectory.Type='ecef'; C.Targets(1).Trajectory.Type='ecef';
        end
    case 'trajectory_state'
        O=SPOT.state(C,C.Observer.Trajectory,C.Time.StartUTC);
        C.Observer.Trajectory.Type='state'; C.Observer.Trajectory.Motion='two_body';
        C.Observer.Trajectory.Position_m=O.Position_m; C.Observer.Trajectory.Velocity_mps=O.Velocity_mps;
        for j=1:numel(C.Targets)
            S=SPOT.state(C,C.Targets(j).Trajectory,C.Time.StartUTC);
            C.Targets(j).Trajectory.Type='state'; C.Targets(j).Trajectory.Motion='two_body';
            C.Targets(j).Trajectory.Position_m=S.Position_m; C.Targets(j).Trajectory.Velocity_mps=S.Velocity_mps;
        end
    case {'trajectory_ephemeris','trajectory_ephemeris_ecef'}
        t=C.Time.StartUTC+seconds([-1 0 .1 .2 1]');
        C.Observer.Trajectory=SPOT_ConfigEphemeris(C,C.Observer.Trajectory,t,strcmp(name,'trajectory_ephemeris_ecef'));
        for j=1:numel(C.Targets), C.Targets(j).Trajectory=SPOT_ConfigEphemeris(C,C.Targets(j).Trajectory,t,strcmp(name,'trajectory_ephemeris_ecef')); end
    case {'trajectory_tle','trajectory_tle_deep_space'}
        prefix='leo'; if strcmp(name,'trajectory_tle_deep_space'), prefix='geo'; end
        C.Targets=C.Targets(1);
        C.Observer.Trajectory.Type='tle'; C.Observer.Trajectory.File=fullfile(C.RootDir,'Inputs','TLE',[prefix '_observer.tle']);
        C.Targets(1).Trajectory.Type='tle'; C.Targets(1).Trajectory.File=fullfile(C.RootDir,'Inputs','TLE',[prefix '_target.tle']);
    case 'trajectory_ephemeris_csv'
        C.Targets=C.Targets(1); C.Observer.Trajectory.Type='ephemeris'; C.Targets(1).Trajectory.Type='ephemeris';
        C.Observer.Trajectory.File=fullfile(C.RootDir,'Inputs','Ephemeris','observer_gcrf.csv');
        C.Targets(1).Trajectory.File=fullfile(C.RootDir,'Inputs','Ephemeris','target_gcrf.csv');
end
end

