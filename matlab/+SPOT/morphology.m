function M=morphology(C,H,w)
%MORPHOLOGY Describe the integrated image; never choose a rendering algorithm.
% sqrt(12*lambda_max) is equivalent uniform straight-trail length from the
% flux-weighted projected motion covariance, including curved/returning paths.
% Changing the label threshold cannot change the image or target signal maps.
p=[H.x(:) H.y(:)]; a=w(:).*H.RadiometricRate(:).*double(H.Visible(:));
ok=all(isfinite(p),2) & isfinite(a) & a>0;
M=struct('Kind','unobserved','Centroid_xy',[NaN NaN],'Covariance_pix2',zeros(2), ...
    'EquivalentSmear_pix',0,'PathLength_pix',0,'EndToEnd_pix',0,'PSFFWHM_pix',0, ...
    'Threshold_pix',0,'DirectionDetector_deg',NaN);
if ~any(ok), return; end
q=p(ok,:); a=a(ok); a=a/sum(a); centre=sum(q.*a,1); d=q-centre;
V=d'*(d.*a); V=(V+V')/2; [U,L]=eig(V); [lambda,index]=max(diag(L));
smear=sqrt(12*max(0,lambda)); axisVector=U(:,index);
step=hypot(diff(p(:,1)),diff(p(:,2))); connected=ok(1:end-1)&ok(2:end);
path=sum(step(connected)); displacement=norm(q(end,:)-q(1,:));
radius=hypot(q(:,1)-(C.Camera.Nx+1)/2,q(:,2)-(C.Camera.Ny+1)/2);
if strcmp(C.Camera.Shape,'circular')
    [xx,yy]=meshgrid(1:C.Camera.Nx,1:C.Camera.Ny);
    xi=(xx-(C.Camera.Nx+1)/2)*C.Camera.FOV_deg(1)/(C.Camera.Nx-1);
    eta=-(yy-(C.Camera.Ny+1)/2)*C.Camera.FOV_deg(2)/(C.Camera.Ny-1);
    active=atand(hypot(tand(xi),tand(eta)))<=C.Camera.FOV_deg(1)/2+1e-10;
    rmax=max(hypot(xx(active)-(C.Camera.Nx+1)/2,yy(active)-(C.Camera.Ny+1)/2));
else
    rmax=hypot((C.Camera.Nx-1)/2,(C.Camera.Ny-1)/2);
end
sigma=zeros(size(radius));
if C.Optics.Enabled && C.Optics.PSFEnabled
    sigma(:)=C.Optics.SigmaCenter_pix;
    if strcmp(C.Optics.PSFModel,'field_dependent_gaussian')
        sigma=sigma+C.Optics.SigmaEdgeAdd_pix*min(1,radius/rmax);
    end
end
fwhm=2*sqrt(2*log(2))*sqrt(sum(a.*sigma.^2));
threshold=max(C.Morphology.MinSmear_pix,C.Morphology.PSFFWHMFactor*fwhm);
kind='object'; if smear>threshold, kind='trail'; end
M.Kind=kind; M.Centroid_xy=centre; M.Covariance_pix2=V;
M.EquivalentSmear_pix=smear; M.PathLength_pix=path; M.EndToEnd_pix=displacement;
M.PSFFWHM_pix=fwhm; M.Threshold_pix=threshold;
M.DirectionDetector_deg=mod(atan2d(axisVector(2),axisVector(1)),180);
end
