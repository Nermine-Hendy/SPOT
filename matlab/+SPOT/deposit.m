function a=deposit(C,x,y,value,active)
%DEPOSIT Flux-conserving bilinear sampling on the active detector, Eq.68.
sz=[C.Camera.Ny C.Camera.Nx]; a=zeros(sz);
x=x(:); y=y(:); value=value(:);
use=isfinite(x)&isfinite(y)&isfinite(value)&value~=0 & x>=1-1e-8 & x<=sz(2)+1e-8 & y>=1-1e-8 & y<=sz(1)+1e-8;
x=min(sz(2),max(1,x(use))); y=min(sz(1),max(1,y(use))); value=value(use);
if isempty(x), return; end
x0=floor(x); y0=floor(y); fx=x-x0; fy=y-y0;
xx=[x0 x0+1 x0 x0+1]; yy=[y0 y0 y0+1 y0+1];
w=[(1-fx).*(1-fy),fx.*(1-fy),(1-fx).*fy,fx.*fy];
valid=xx>=1 & xx<=sz(2) & yy>=1 & yy<=sz(1);
idx=sub2ind(sz,min(sz(1),max(1,yy)),min(sz(2),max(1,xx)));
w=w.*valid.*active(idx); den=sum(w,2); w=w./max(den,realmin);
vals=w.*value; a=reshape(accumarray(idx(:),vals(:),[prod(sz) 1],@sum,0),sz);
end
