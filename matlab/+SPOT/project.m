function P=project(C,A,los)
%PROJECT Exact forward mapping paired with detectorGrid, Eqs.27-40.
d=A'*los; xi=atan2d(d(1,:),d(3,:)); eta=atan2d(-d(2,:),d(3,:));
P.Xi_deg=xi; P.Eta_deg=eta;
P.x=(C.Camera.Nx+1)/2+(C.Camera.Nx-1)*xi/C.Camera.FOV_deg(1);
P.y=(C.Camera.Ny+1)/2-(C.Camera.Ny-1)*eta/C.Camera.FOV_deg(2);
P.InFOV=d(3,:)>0 & abs(xi)<=C.Camera.FOV_deg(1)/2+1e-10 & abs(eta)<=C.Camera.FOV_deg(2)/2+1e-10;
if strcmp(C.Camera.Shape,'circular'), P.InFOV=P.InFOV & acosd(max(-1,min(1,d(3,:))))<=C.Camera.FOV_deg(1)/2+1e-10; end
end
