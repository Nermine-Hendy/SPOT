function [t,w]=timeGrid(startUTC,exposure,step,extra)
%TIMEGRID Endpoint-complete grid with an exact midpoint and normalized trapezoids.
if nargin<4, extra=[]; end
q=unique([0:step:exposure exposure exposure/2 reshape(extra,1,[])]);
q=q(q>=0 & q<=exposure);
q=q([true diff(q)>max(64*eps(exposure),1e-12*exposure)]);
q(1)=0; q(end)=exposure; [~,middle]=min(abs(q-exposure/2)); q(middle)=exposure/2;
t=startUTC+seconds(q(:));
d=diff(q); w=[d(1)/2 (d(1:end-1)+d(2:end))/2 d(end)/2]/exposure;
assert(abs(sum(w)-1)<1e-10,'SPOT:Integration','Quadrature weights do not sum to one.');
end
