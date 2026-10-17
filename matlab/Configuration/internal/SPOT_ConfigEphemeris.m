function T=SPOT_ConfigEphemeris(C,T,t,ecef)
S=SPOT.state(C,T,t); 
T.Type='ephemeris'; 
T.TimeUTC=t; 
T.Positions_m=S.Position_m; 
T.Velocities_mps=S.Velocity_mps;
T.File=''; 
T.Frame='GCRF';
if ecef
    T.Frame='ECEF'; T.Positions_m=S.ECEF_m;
    % Position-only ephemeris takes the documented piecewise-linear path.
    T.Velocities_mps=[]; T.Interpolation='linear';
end
end
