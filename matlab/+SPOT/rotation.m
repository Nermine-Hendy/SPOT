function R = rotation(C,t)
%ROTATION ECEF -> GCRF, with explicitly selected Earth orientation accuracy.
t=t(:); assert(isdatetime(t) && ~isempty(t.TimeZone),'SPOT:Time','Frame conversion needs timezone-aware times.');
t.TimeZone='UTC'; n=numel(t); R=zeros(3,3,n);
switch lower(C.Frames.Model)
    case 'iau2006'
        if exist('dcmeci2ecef','file')~=2
            error('SPOT:Dependency','IAU2006 needs Aerospace Toolbox. For approximate work explicitly select Frames.Model=''gmst_precession''.');
        end
        % The backend requires one EOP row per UTC row, even for constant EOP.
        % Rebuild these arrays for every query: adaptive grids have different N.
        validateattributes(C.Frames.DeltaAT_s,{'numeric'},{'real','finite','scalar'});
        validateattributes(C.Frames.DeltaUT1_s,{'numeric'},{'real','finite','scalar'});
        validateattributes(C.Frames.PolarMotion_rad,{'numeric'},{'real','finite','vector','numel',2});
        deltaAT=repmat(double(C.Frames.DeltaAT_s),n,1);
        deltaUT1=repmat(double(C.Frames.DeltaUT1_s),n,1);
        polarMotion=repmat(reshape(double(C.Frames.PolarMotion_rad),1,2),n,1);
        D=dcmeci2ecef('IAU-2000/2006',datevec(t),deltaAT,deltaUT1,polarMotion);
        for k=1:n, R(:,:,k)=D(:,:,k).'; end
    case 'gmst_precession'
        for k=1:n
            jd=juliandate(t(k)); T=(jd-2451545)/36525;
            th=mod(280.46061837+360.98564736629*(jd-2451545)+ ...
                .000387933*T^2-T^3/38710000,360);
            Q=SPOT.precession(t(k)); % GCRF/J2000 -> mean equator/equinox of date.
            R(:,:,k)=Q.'*[cosd(th) -sind(th) 0;sind(th) cosd(th) 0;0 0 1];
        end
    otherwise, error('SPOT:Frame','Unknown Frames.Model %s.',C.Frames.Model);
end
end
