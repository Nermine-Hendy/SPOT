function S = state(C,T,t)
%STATE Common observer/target trajectory evaluator. Output is GCRF, m, m/s.
t=t(:); n=numel(t); dt=seconds(t-T.EpochUTC).';
r=zeros(3,n); v=r; inputFrame=T.Frame;
switch lower(T.Type)
    case 'geodetic'
        g=T.Geodetic_deg_m(:); r0=SPOT.geodetic(C,g,'to_ecef');
        ve=SPOT.enu(g(1),g(2))*T.VelocityENU_mps(:); inputFrame='ECEF';
        if strcmp(T.Motion,'fixed'), ve(:)=0; end
        r=r0+ve*dt; v=repmat(ve,1,n);
    case {'ecef','state'}
        if strcmp(T.Type,'ecef'), inputFrame='ECEF'; end
        r0=T.Position_m(:); v0=T.Velocity_mps(:);
        switch lower(T.Motion)
            case 'fixed', r=repmat(r0,1,n);
            case 'linear', r=r0+v0*dt; v=repmat(v0,1,n);
            case 'two_body'
                if ~strcmp(inputFrame,'GCRF'), error('SPOT:Frame','Two-body Cartesian state must be GCRF.'); end
                mu=C.Constants.Mu_m3_s2; y0=[r0;v0];
                ode=@(~,y)[y(4:6);-mu*y(1:3)/norm(y(1:3))^3];
                opts=odeset('RelTol',1e-11,'AbsTol',1e-5);
                for j=1:n
                    if dt(j)==0, y=y0; else, sol=ode45(ode,[0 dt(j)],y0,opts); y=deval(sol,dt(j)); end
                    r(:,j)=y(1:3); v(:,j)=y(4:6);
                end
            otherwise, error('SPOT:Trajectory','State Motion must be fixed, linear or two_body.');
        end
    case 'keplerian'
        el=T.Elements; a=el(1); e=el(2); nu=deg2rad(el(6)); mu=C.Constants.Mu_m3_s2;
        E0=2*atan2(sqrt(1-e)*sin(nu/2),sqrt(1+e)*cos(nu/2));
        M=E0-e*sin(E0)+sqrt(mu/a^3)*dt; E=M;
        for j=1:40
            correction=(E-e*sin(E)-M)./(1-e*cos(E)); E=E-correction;
            if max(abs(correction),[],'all')<1e-13, break; end
        end
        if max(abs(E-e*sin(E)-M),[],'all')>1e-10, error('SPOT:Orbit','Kepler solver failed to converge.'); end
        q=sqrt(1-e^2); rp=[a*(cos(E)-e);a*q*sin(E);zeros(1,n)];
        edot=sqrt(mu/a^3)./(1-e*cos(E)); vp=[-a*sin(E).*edot;a*q*cos(E).*edot;zeros(1,n)];
        Rz=@(ang)[cosd(ang) -sind(ang) 0;sind(ang) cosd(ang) 0;0 0 1];
        Rx=@(ang)[1 0 0;0 cosd(ang) -sind(ang);0 sind(ang) cosd(ang)];
        Q=Rz(el(4))*Rx(el(3))*Rz(el(5)); r=Q*rp; v=Q*vp; inputFrame='GCRF';
    case 'ephemeris'
        x=seconds(T.TimeUTC-T.EpochUTC); q=dt(:);
        if strcmp(T.Extrapolation,'error') && (min(q)<x(1)-1e-7 || max(q)>x(end)+1e-7)
            error('SPOT:EphemerisRange','Query times lie outside the supplied ephemeris.');
        end
        q=min(max(q,x(1)),x(end)); % only tolerance-level excursions accepted.
        for j=1:n
            k=find(x<=q(j),1,'last'); k=min(k,numel(x)-1);
            h=x(k+1)-x(k); u=(q(j)-x(k))/h; p0=T.Positions_m(:,k); p1=T.Positions_m(:,k+1);
            if strcmp(T.Interpolation,'hermite') && ~isempty(T.Velocities_mps)
                v0=T.Velocities_mps(:,k); v1=T.Velocities_mps(:,k+1);
                r(:,j)=(2*u^3-3*u^2+1)*p0+(u^3-2*u^2+u)*h*v0+(-2*u^3+3*u^2)*p1+(u^3-u^2)*h*v1;
                v(:,j)=((6*u^2-6*u)*p0+(-6*u^2+6*u)*p1)/h+(3*u^2-4*u+1)*v0+(3*u^2-2*u)*v1;
            else
                r(:,j)=(1-u)*p0+u*p1; v(:,j)=(p1-p0)/h;
            end
        end
    case 'tle'
        for j=1:n
            [r(:,j),v(:,j)]=states(T.Runtime.Satellite,t(j),'CoordinateFrame','inertial');
        end
        inputFrame='GCRF';
    otherwise, error('SPOT:Trajectory','Demo states are evaluated through SPOT.sample, not the physical trajectory backend.');
end
R=SPOT.rotation(C,t); re=zeros(3,n);
if strcmp(inputFrame,'ECEF')
    re=r; Rp=SPOT.rotation(C,t+seconds(.5)); Rm=SPOT.rotation(C,t-seconds(.5));
    for j=1:n, v(:,j)=R(:,:,j)*v(:,j)+(Rp(:,:,j)-Rm(:,:,j))*r(:,j); r(:,j)=R(:,:,j)*r(:,j); end
elseif strcmp(inputFrame,'GCRF')
    for j=1:n, re(:,j)=R(:,:,j).'*r(:,j); end
else, error('SPOT:Frame','State frame must be GCRF or ECEF.');
end
if any(~isfinite(r(:))) || any(~isfinite(v(:))), error('SPOT:State','Nonfinite trajectory output.'); end
S=struct('TimeUTC',t,'Position_m',r,'Velocity_mps',v,'ECEF_m',re, ...
    'Geodetic_deg_m',SPOT.geodetic(C,re,'from_ecef'),'Frame','GCRF');
end
