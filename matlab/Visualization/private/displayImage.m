function [a,lim]=displayImage(a,C)
x=sort(a(isfinite(a))); p=C.Plot.DisplayPercentiles;
q=1+(numel(x)-1)*p/100; lim=interp1(1:numel(x),x,q,'linear');
if lim(2)<=lim(1), lim(2)=lim(1)+1; end
a=max(0,min(1,(double(a)-lim(1))/diff(lim)));
end
