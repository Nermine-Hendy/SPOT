function [out,S]=optics(C,rate,active)
%OPTICS One linear operator for scene AND every target truth map.
% Input-position-weighted PSF banks conserve interior flux; detector-edge
% losses are physical. Output-position kernel selection is deliberately absent.
[x,y]=meshgrid(1:C.Camera.Nx,1:C.Camera.Ny);
r=hypot(x-(C.Camera.Nx+1)/2,y-(C.Camera.Ny+1)/2); r=r/max(r(active));
r=min(1,r); tau=1; vig=ones(size(rate));
if C.Optics.Enabled
    tau=C.Optics.Transmission;
    if ~isempty(C.Optics.ElementTransmission), tau=prod(C.Optics.ElementTransmission); end
    if C.Optics.VignettingEnabled, vig=1-(1-C.Optics.EdgeTransmission)*r.^C.Optics.RadialExponent; end
end
S.Integrated=rate*C.Time.Exposure_s; S.Throughput=S.Integrated*tau;
S.Vignetted=S.Throughput.*vig; out=S.Vignetted;
if C.Optics.Enabled && C.Optics.PSFEnabled
    if strcmp(C.Optics.PSFModel,'gaussian') || C.Optics.SigmaEdgeAdd_pix==0
        out=blur(out,C.Optics.SigmaCenter_pix);
    else
        out=zeros(size(rate)); nodes=linspace(0,1,C.Optics.PSFBanks);
        for k=1:numel(nodes)
            w=max(0,1-abs(r-nodes(k))*(numel(nodes)-1));
            if ~any(S.Vignetted(w>0),'all'), continue; end
            out=out+blur(S.Vignetted.*w,C.Optics.SigmaCenter_pix+C.Optics.SigmaEdgeAdd_pix*nodes(k));
        end
    end
end
out=max(0,out).*active; S.PSF=out; S.Vignetting=vig; S.ThroughputFactor=tau;
end
function out=blur(a,sigma)
if sigma<=1e-12, out=a; return; end
q=-ceil(4*sigma):ceil(4*sigma); h=exp(-.5*(q/sigma).^2); h=h/sum(h);
out=conv2(h',h,a,'same');
end
