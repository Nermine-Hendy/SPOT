function A=pointing(C,O,t,targetPosition)
%POINTING Camera-to-GCRF triad [column-right,row-down,boresight].
n=numel(t); A=zeros(3,3,n);
if strcmp(C.Pointing.Mode,'attitude')
    tt=seconds(C.Pointing.AttitudeUTC(:)-C.Pointing.AttitudeUTC(1));
    q=seconds(t(:)-C.Pointing.AttitudeUTC(1)); support(tt,q,'attitude');
    for k=1:n
        i=find(tt<=q(k),1,'last'); i=min(i,numel(tt)-1); f=(q(k)-tt(i))/(tt(i+1)-tt(i));
        B=C.Pointing.CameraToGCRF(:,:,i); D=B'*C.Pointing.CameraToGCRF(:,:,i+1);
        % Geodesic interpolation on SO(3), including 180-degree rotations.
        [V,L]=eig(D); [~,z]=min(abs(diag(L)-1)); axis=real(V(:,z)); axis=axis/norm(axis);
        ang=acos(max(-1,min(1,(trace(D)-1)/2)));
        if dot(axis,[D(3,2)-D(2,3);D(1,3)-D(3,1);D(2,1)-D(1,2)])<0, axis=-axis; end
        K=[0 -axis(3) axis(2);axis(3) 0 -axis(1);-axis(2) axis(1) 0];
        A(:,:,k)=B*(eye(3)+sin(f*ang)*K+(1-cos(f*ang))*K*K);
    end
else
    R=SPOT.rotation(C,t);
    for k=1:n
        E=R(:,:,k)*SPOT.enu(O.Geodetic_deg_m(1,k),O.Geodetic_deg_m(2,k));
        switch C.Pointing.Mode
            case 'az_el_roll'
                az=C.Pointing.AzEl_deg(1); el=C.Pointing.AzEl_deg(2);
                b=E*[cosd(el)*sind(az);cosd(el)*cosd(az);sind(el)]; a=E(:,3);
            case 'nadir_offset'
                nadir=-unit(O.Position_m(:,k)); normal=cross(O.Position_m(:,k),O.Velocity_mps(:,k));
                assert(norm(normal)>1e-3,'SPOT:Pointing','Nadir/orbital azimuth pointing needs a nonradial observer velocity.');
                normal=unit(normal); along=unit(cross(normal,-nadir));
                b=cosd(C.Pointing.OffNadir_deg)*nadir+sind(C.Pointing.OffNadir_deg)*(cosd(C.Pointing.SpaceAzimuth_deg)*along+sind(C.Pointing.SpaceAzimuth_deg)*normal); a=normal;
            case 'inertial'
                b=unit(C.Pointing.Boresight_GCRF); a=unit(C.Pointing.ReferenceUp_GCRF);
            case 'target_tracking'
                b=unit(targetPosition(:,k)-O.Position_m(:,k)); a=unit(C.Pointing.ReferenceUp_GCRF);
        end
        if norm(cross(a,b))<1e-8
            a=E(:,2); if norm(cross(a,b))<1e-8, [~,ii]=min(abs(b)); a=zeros(3,1); a(ii)=1; end
        end
        v=unit(-a+dot(a,b)*b); u=unit(cross(v,b)); roll=C.Pointing.Roll_deg;
        A(:,:,k)=[cosd(roll)*u+sind(roll)*v,-sind(roll)*u+cosd(roll)*v,b];
    end
end
if ~isempty(C.Pointing.JitterUTC)
    tt=seconds(C.Pointing.JitterUTC(:)-C.Pointing.JitterUTC(1)); q=seconds(t(:)-C.Pointing.JitterUTC(1)); support(tt,q,'jitter');
    j=interp1(tt,C.Pointing.Jitter_arcsec,q,'linear')*pi/(180*3600);
    for k=1:n
        a=j(k,1); b=j(k,2); c=j(k,3);
        Ry=[cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)];
        Rx=[1 0 0;0 cos(b) sin(b);0 -sin(b) cos(b)];
        Rz=[cos(c) -sin(c) 0;sin(c) cos(c) 0;0 0 1];
        A(:,:,k)=A(:,:,k)*Ry*Rx*Rz;
    end
end
end
function x=unit(x)
x=x(:); assert(norm(x)>0 && all(isfinite(x)),'SPOT:Pointing','Cannot normalize a zero or invalid vector.'); x=x/norm(x);
end
function support(x,q,label)
assert(all(q>=x(1)-1e-7 & q<=x(end)+1e-7),'SPOT:Time','Requested time lies outside %s samples.',label);
end
