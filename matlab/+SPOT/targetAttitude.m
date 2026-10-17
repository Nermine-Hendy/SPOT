function R=targetAttitude(A,t)
%TARGETATTITUDE Body-to-GCRF active rotations, independent of camera pointing.
t=t(:); R=zeros(3,3,numel(t));
switch A.Mode
    case {'fixed','spin'}
        axis=A.SpinAxisBody(:)/norm(A.SpinAxisBody);
        K=[0 -axis(3) axis(2);axis(3) 0 -axis(1);-axis(2) axis(1) 0];
        for k=1:numel(t)
            angle=0;
            if strcmp(A.Mode,'spin'), angle=deg2rad(A.SpinRate_deg_s)*seconds(t(k)-A.EpochUTC); end
            R(:,:,k)=A.BodyToGCRF0*(eye(3)+sin(angle)*K+(1-cos(angle))*(K*K));
        end
    case 'samples'
        x=seconds(A.TimeUTC(:)-A.TimeUTC(1)); q=seconds(t-A.TimeUTC(1));
        assert(all(q>=x(1)-1e-7 & q<=x(end)+1e-7),'SPOT:TargetAttitude','Target attitude does not cover exposure.');
        for k=1:numel(t)
            z=min(max(q(k),x(1)),x(end)); i=find(x<=z,1,'last'); i=min(i,numel(x)-1);
            a=(z-x(i))/(x(i+1)-x(i)); R0=A.BodyToGCRF(:,:,i); R1=A.BodyToGCRF(:,:,i+1);
            % Principal SO(3) interpolation; rotations between knots must be <180 deg.
            L=real(logm(R0'*R1)); L=(L-L')/2; R(:,:,k)=R0*expm(a*L);
        end
    otherwise, error('SPOT:TargetAttitude','Unsupported target attitude mode.');
end
end
