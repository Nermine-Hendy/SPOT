function [G,w,requested]=integrateGeometry(C,startUTC)
%INTEGRATEGEOMETRY Spatial refinement plus intrinsic brightness error estimates.
% A spin/lobe-width step cap guards against missing narrow specular peaks.
% This is an error estimator, not a proof: use SPOT.convergence on final scenes.
duration=C.Time.Exposure_s; mid=startUTC+seconds(duration/2); extra=[];
series={C.Pointing.AttitudeUTC,C.Pointing.JitterUTC,C.Observer.Trajectory.TimeUTC};
stepLimit=C.Integration.Step_s;
for j=1:numel(C.Targets)
    T=C.Targets(j); series=[series {T.Trajectory.TimeUTC,T.Photometry.LightCurveUTC}]; %#ok<AGROW>
    if strcmp(T.Photometry.Model,'surface_brdf')
        A=T.Attitude; series=[series {A.TimeUTC}]; %#ok<AGROW>
        speed=0;
        if strcmp(A.Mode,'spin'), speed=abs(A.SpinRate_deg_s); end
        if strcmp(A.Mode,'samples')
            for k=1:numel(A.TimeUTC)-1
                angle=acosd(max(-1,min(1,(trace(A.BodyToGCRF(:,:,k)'*A.BodyToGCRF(:,:,k+1))-1)/2)));
                speed=max(speed,angle/seconds(A.TimeUTC(k+1)-A.TimeUTC(k)));
            end
        end
        if speed>0
            stepLimit=min(stepLimit,C.Integration.MaxBodyRotation_deg/speed);
            use=T.Surface.SpecularReflectance>0;
            if any(use)
                p=max(T.Surface.SpecularExponent(use)); halfWidth=acosd(2^(-1/p));
                stepLimit=min(stepLimit,halfWidth/(2*speed*C.Integration.SpecularSamplesPerHWHM));
            end
        end
    end
end
for j=1:numel(series)
    if ~isempty(series{j}), x=seconds(series{j}(:)-startUTC); extra=[extra;x(x>=0 & x<=duration)]; end %#ok<AGROW>
end
assert(ceil(duration/stepLimit)+1<=C.Integration.MaxSamples,'SPOT:Sampling','Spin/glint sampling exceeds MaxSamples; shorten exposure or increase budget.');
[requested,~]=SPOT.timeGrid(startUTC,duration,C.Integration.Step_s);
[t,w]=SPOT.timeGrid(startUTC,duration,stepLimit,extra);
for iteration=1:20
    assert(numel(t)<=C.Integration.MaxSamples,'SPOT:Sampling','Target sample budget exceeded.');
    G=SPOT.sample(C,t,mid); step=zeros(1,numel(t)-1); specularSplit=false(size(step));
    for j=1:numel(G.Target)
        H=G.Target{j}; u=zeros(3,numel(t));
        for k=1:numel(t), u(:,k)=G.CameraToGCRF(:,:,k)'*H.LOS_GCRF(:,k); end
        ang=acosd(max(-1,min(1,sum(u(:,1:end-1).*u(:,2:end),1))));
        scale=max(([C.Camera.Nx C.Camera.Ny]-1)./C.Camera.FOV_deg)/cosd(max(C.Camera.FOV_deg)/2)^2;
        d=ang*scale; near=u(3,:)>cosd(min(89,max(C.Camera.FOV_deg)));
        d(~(near(1:end-1)|near(2:end)))=0; step=max(step,d);
        T=C.Targets(j);
        if strcmp(T.Photometry.Model,'surface_brdf') && any(T.Surface.SpecularReflectance>0)
            su=G.SunPosition_m-H.Position_m; su=su./vecnorm(su);
            sunTurn=acosd(max(-1,min(1,sum(su(:,1:end-1).*su(:,2:end),1))));
            viewTurn=acosd(max(-1,min(1,sum(H.LOS_GCRF(:,1:end-1).*H.LOS_GCRF(:,2:end),1))));
            R=SPOT.targetAttitude(T.Attitude,t); bodyTurn=zeros(size(step));
            for k=1:numel(step)
                bodyTurn(k)=acosd(max(-1,min(1,(trace(R(:,:,k)'*R(:,:,k+1))-1)/2)));
            end
            power=max(T.Surface.SpecularExponent(T.Surface.SpecularReflectance>0));
            angularLimit=acosd(2^(-1/power))/C.Integration.SpecularSamplesPerHWHM;
            specularSplit=specularSplit | 2*bodyTurn+sunTurn+viewTurn>angularLimit;
        end
    end
    q=seconds(t-startUTC)'; qm=(q(1:end-1)+q(2:end))/2;
    split=specularSplit;
    if C.Integration.Adaptive, split=split | step>C.Integration.MaxStep_pix*1.001; end
    if C.Integration.PhotometryAdaptive && ~isempty(C.Targets)
        M=SPOT.sample(C,startUTC+seconds(qm(:)),mid);
        for j=1:numel(C.Targets)
            a=SPOT.intrinsicPhotometry(C,G,j); b=SPOT.intrinsicPhotometry(C,M,j);
            % Visibility gate contributes to integrated signal, including ingress.
            f=a.IntrinsicRate.*a.Visible; fm=b.IntrinsicRate.*b.Visible;
            coarse=(f(1:end-1)+f(2:end))/2;
            fine=(f(1:end-1)+2*fm+f(2:end))/4;
            % Global peak absolute scale prevents endless refinement at a step.
            tol=C.Integration.PhotometryAbsTol+C.Integration.PhotometryRelTol*max([f fm]);
            err=abs(fine-coarse).*diff(q)/duration;
            need=err>tol*diff(q)/duration;
            % A discontinuity may never meet pointwise error: allow tiny intervals.
            need(err<=tol & diff(q)<=duration*1e-5)=false;
            split=split | need;
        end
    end
    if ~any(split), return; end
    assert(numel(t)+nnz(split)<=C.Integration.MaxSamples,'SPOT:Sampling','Spatial/photometric refinement exceeds MaxSamples. Increase budget or relax tolerances explicitly.');
    [t,w]=SPOT.timeGrid(startUTC,duration,duration,[q qm(split)]);
end
error('SPOT:Sampling','Adaptive target sampling did not converge.');
end
