function G=sample(C,t,midUTC)
%SAMPLE Shared time-indexed observer, target, pointing and visibility geometry.
t=t(:); n=numel(t); nt=numel(C.Targets); O=SPOT.state(C,C.Observer.Trajectory,t);
R=SPOT.rotation(C,t); up=zeros(3,n);
for k=1:n, E=SPOT.enu(O.Geodetic_deg_m(1,k),O.Geodetic_deg_m(2,k)); up(:,k)=R(:,:,k)*E(:,3); end
T=cell(1,nt); track=[];
for j=1:nt
    if ~strcmp(C.Targets(j).Trajectory.Type,'demo'), T{j}=SPOT.state(C,C.Targets(j).Trajectory,t); end
    if C.Targets(j).ID==C.Pointing.TrackTargetID && ~isempty(T{j}), track=T{j}.Position_m; end
end
A=SPOT.pointing(C,O,t,track);
if strcmp(C.Integration.AttitudeMode,'midpoint')
    Om=SPOT.state(C,C.Observer.Trajectory,midUTC); tm=[];
    if strcmp(C.Pointing.Mode,'target_tracking')
        j=find([C.Targets.ID]==C.Pointing.TrackTargetID,1); s=SPOT.state(C,C.Targets(j).Trajectory,midUTC); tm=s.Position_m;
    end
    A=repmat(SPOT.pointing(C,Om,midUTC,tm),1,1,n);
end
G.TimeUTC=t; G.Observer=O; G.CameraToGCRF=A; G.Zenith_GCRF=up;
G.IsGround=O.Geodetic_deg_m(3,:)<C.Visibility.PlatformAltitude_m;
G.SunPosition_m=SPOT.sun(C,t); G.Target=cell(1,nt);
for j=1:nt
    tr=C.Targets(j).Trajectory;
    if strcmp(tr.Type,'demo')
        r=zeros(3,n); dt=seconds(t-tr.Demo.ReferenceUTC)';
        for k=1:n
            d=tr.Demo; off=d.Offset_deg;
            switch d.AngleMode
                case 'az_el_offset'
                    E=R(:,:,k)*SPOT.enu(O.Geodetic_deg_m(1,k),O.Geodetic_deg_m(2,k));
                    az=C.Pointing.AzEl_deg(1)+off(1); el=C.Pointing.AzEl_deg(2)+off(2);
                    initial=E*[cosd(el)*sind(az);cosd(el)*cosd(az);sind(el)];
                    p=SPOT.project(C,A(:,:,k),initial); xy=[p.Xi_deg p.Eta_deg];
                case 'theta_phi'
                    xy=atand(tand(off(1))*[cosd(off(2)) sind(off(2))]);
                case 'camera_angles', xy=off;
            end
            if strcmp(d.MotionModel,'linear_angles')
                xy=xy+d.AngularRate_deg_s*dt(k);
                assert(all(abs(xy)<89.9),'SPOT:Demo','Linear camera-angle demo left its forward hemisphere. Use great_circle demo motion or a physical trajectory for long acquisitions.');
                local=[tand(xy(1));-tand(xy(2));1]; local=local/norm(local);
            else
                % Globally continuous artificial path. Initial position and
                % instantaneous angular rates agree with the CSV seed.
                raw=[tand(xy(1));-tand(xy(2));1]; local=raw/norm(raw);
                derivative=[d.AngularRate_deg_s(1)/cosd(xy(1))^2;-d.AngularRate_deg_s(2)/cosd(xy(2))^2;0]*pi/180;
                derivative=(derivative-local*dot(local,derivative))/norm(raw); speed=norm(derivative);
                if speed>0, local=local*cos(speed*dt(k))+derivative/speed*sin(speed*dt(k)); end
            end
            u=A(:,:,k)*local; u=u/norm(u);
            r(:,k)=O.Position_m(:,k)+d.Range_m*u;
        end
        velocity=nan(3,n);
        if n>1
            querySeconds=seconds(t-t(1));
            for component=1:3, velocity(component,:)=gradient(r(component,:),querySeconds); end
        end
        T{j}=struct('Position_m',r,'Velocity_mps',velocity,'Frame','GCRF');
    end
    r=T{j}.Position_m; dr=r-O.Position_m; range=vecnorm(dr); u=dr./max(range,realmin);
    H=struct('ID',C.Targets(j).ID,'Position_m',r,'Velocity_mps',T{j}.Velocity_mps, ...
        'Range_m',range,'LOS_GCRF',u,'x',zeros(1,n),'y',zeros(1,n), ...
        'Xi_deg',zeros(1,n),'Eta_deg',zeros(1,n),'InFOV',false(1,n));
    for k=1:n
        p=SPOT.project(C,A(:,:,k),u(:,k));
        H.x(k)=p.x; H.y(k)=p.y; H.Xi_deg(k)=p.Xi_deg; H.Eta_deg(k)=p.Eta_deg; H.InFOV(k)=p.InFOV;
    end
    H.Elevation_deg=asind(max(-1,min(1,sum(up.*u,1))));
    H.EarthOccluded=false(1,n);
    if C.Visibility.EarthOccultation
        s=-sum(O.Position_m.*u,1); closest=O.Position_m+max(0,min(range,s)).*u;
        H.EarthOccluded=vecnorm(closest)<C.Constants.EarthRadius_m & ~G.IsGround;
    end
    H.PlatformVisible=(~G.IsGround | H.Elevation_deg>=C.Visibility.MinElevation_deg) & ~H.EarthOccluded;
    H.RangeValid=range>=C.Visibility.RangeLimits_m(1) & range<=C.Visibility.RangeLimits_m(2);
    su=G.SunPosition_m-O.Position_m; su=su./vecnorm(su);
    H.SunExcluded=acosd(max(-1,min(1,sum(u.*su,1))))<C.Visibility.SunExclusion_deg;
    H.Visible=H.InFOV & H.PlatformVisible & H.RangeValid & ~H.SunExcluded;
    H.Illumination=SPOT.illumination(C,r,G.SunPosition_m);
    H.PhaseAngle_rad=acos(max(-1,min(1,sum((G.SunPosition_m-r)./vecnorm(G.SunPosition_m-r).*(-u),1))));
    G.Target{j}=H;
end
end
