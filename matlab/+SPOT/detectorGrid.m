function D=detectorGrid(C,A)
%DETECTORGRID Equiangular centres, inverse projection and active solid angle.
[xi,eta]=meshgrid(linspace(-C.Camera.FOV_deg(1)/2,C.Camera.FOV_deg(1)/2,C.Camera.Nx), ...
    linspace(C.Camera.FOV_deg(2)/2,-C.Camera.FOV_deg(2)/2,C.Camera.Ny));
v=[tand(xi(:))';-tand(eta(:))';ones(1,numel(xi))]; v=v./vecnorm(v);
active=true(size(xi));
if strcmp(C.Camera.Shape,'circular'), active=reshape(acosd(v(3,:))<=C.Camera.FOV_deg(1)/2+1e-10,size(xi)); end
assert(any(active,'all'),'SPOT:Camera','The selected FoV contains no active detector pixel centres. Increase detector dimensions.');
f=deg2rad(C.Camera.FOV_deg);
if strcmp(C.Camera.Shape,'circular'), omega=2*pi*(1-cos(f(1)/2));
else, omega=4*asin(sin(f(1)/2)*sin(f(2)/2)); end
if strcmp(C.Camera.SolidAngle,'mean')
    w=double(active)*omega/nnz(active);
else
    % Local angular Jacobian; normalize sampled weights to exact total FoV.
    w=(1+tand(xi).^2).*(1+tand(eta).^2)./(1+tand(xi).^2+tand(eta).^2).^1.5;
    % Boundary pixel centres occupy half-cells at the stated FoV boundary.
    w([1 end],:)=w([1 end],:)/2; w(:,[1 end])=w(:,[1 end])/2;
    w=w.*active; w=w*(omega/sum(w,'all'));
end
los=A*v;
D=struct('Xi_deg',xi,'Eta_deg',eta,'Active',active,'SolidAngle_sr',w, ...
    'FOVSolidAngle_sr',omega,'LOS_GCRF',los, ...
    'RA_deg',reshape(mod(atan2d(los(2,:),los(1,:)),360),size(xi)), ...
    'DEC_deg',reshape(asind(max(-1,min(1,los(3,:)))),size(xi)));
% Independent coarser angular diagnostic grid (Table III); detector truth above
% always stays at the exact detector resolution.
step=C.Camera.DiagnosticGridStep_deg;
gx=unique([-C.Camera.FOV_deg(1)/2:step:C.Camera.FOV_deg(1)/2 C.Camera.FOV_deg(1)/2]);
gy=unique([-C.Camera.FOV_deg(2)/2:step:C.Camera.FOV_deg(2)/2 C.Camera.FOV_deg(2)/2]);
[xx,yy]=meshgrid(gx,fliplr(gy)); v=[tand(xx(:))';-tand(yy(:))';ones(1,numel(xx))]; v=v./vecnorm(v); u=A*v;
D.AngularDiagnosticGrid=struct('Xi_deg',xx,'Eta_deg',yy, ...
    'RA_deg',reshape(mod(atan2d(u(2,:),u(1,:)),360),size(xx)), ...
    'DEC_deg',reshape(asind(max(-1,min(1,u(3,:)))),size(xx)));
end
