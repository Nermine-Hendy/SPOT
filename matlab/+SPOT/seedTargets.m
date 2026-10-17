function C=seedTargets(C)
%SEEDTARGETS Resolve seeds ONCE into physical GCRF Cartesian states.
% Velocity is relative to the observer at the epoch, not a desired pixel length.
use=find(arrayfun(@(a) strcmp(a.Trajectory.Type,'seeded_state'),C.Targets));
if isempty(use), return; end
assert(~strcmp(C.Pointing.Mode,'target_tracking'),'SPOT:Seed','Resolve seeds before selecting target tracking, or supply explicit states.');
t=C.Time.StartUTC; O=SPOT.state(C,C.Observer.Trajectory,t);
A=SPOT.pointing(C,O,t,[]); R=SPOT.rotation(C,t);
for j=use
    d=C.Targets(j).Seed;
    validateattributes(d.Range_m,{'numeric'},{'scalar','finite','positive'});
    validateattributes(d.TransverseSpeed_mps,{'numeric'},{'scalar','finite','nonnegative'});
    validateattributes(d.RadialSpeed_mps,{'numeric'},{'scalar','finite'});
    validateattributes(d.Direction_deg,{'numeric'},{'scalar','finite'});
    validateattributes(d.Angles_deg,{'numeric'},{'vector','numel',2,'finite'});
    switch d.AngleMode
        case 'az_el_offset'
            assert(strcmp(C.Pointing.Mode,'az_el_roll'),'SPOT:Seed','Az/el offsets require az_el_roll pointing.');
            E=R*SPOT.enu(O.Geodetic_deg_m(1),O.Geodetic_deg_m(2));
            ae=reshape(C.Pointing.AzEl_deg,1,2)+reshape(d.Angles_deg,1,2);
            p=SPOT.project(C,A,E*[cosd(ae(2))*sind(ae(1));cosd(ae(2))*cosd(ae(1));sind(ae(2))]);
            xy=[p.Xi_deg p.Eta_deg];
        case 'theta_phi'
            assert(d.Angles_deg(1)>=0 && d.Angles_deg(1)<89,'SPOT:Seed','Theta must be in [0,89) degrees.');
            xy=atand(tand(d.Angles_deg(1))*[cosd(d.Angles_deg(2)) sind(d.Angles_deg(2))]);
        case 'camera_angles', xy=d.Angles_deg;
        otherwise, error('SPOT:Seed','Unknown seed AngleMode.');
    end
    assert(all(abs(xy)<89),'SPOT:Seed','Seed must lie in the forward camera hemisphere.');
    q=[tand(xy(1));-tand(xy(2));1]; u=q/norm(q);
    v=[cosd(d.Direction_deg)/cosd(xy(1))^2;-sind(d.Direction_deg)/cosd(xy(2))^2;0];
    v=v-u*dot(u,v); v=v/norm(v);
    T=C.Targets(j).Trajectory; T.Type='state'; T.Frame='GCRF'; T.EpochUTC=t;
    T.Motion=C.Seeding.Propagation;
    T.Position_m=O.Position_m+d.Range_m*A*u;
    T.Velocity_mps=O.Velocity_mps+A*(d.TransverseSpeed_mps*v+d.RadialSpeed_mps*u);
    C.Targets(j).Trajectory=T;
    C.Targets(j).Metadata.Seed=d; C.Targets(j).Metadata.SeedEpochUTC=t;
    C.Targets(j).Metadata.SeedVelocityConvention='Observer-relative velocity at epoch; camera rotation is not subtracted.';
end
end
