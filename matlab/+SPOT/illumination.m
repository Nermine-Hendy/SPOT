function I=illumination(C,r,sun)
%ILLUMINATION Geometric cylindrical eclipse or finite angular-disk overlap.
RE=C.Constants.EarthRadius_m; rr=vecnorm(r); s=sun./vecnorm(sun);
behind=sum(r.*s,1)<0; distance=vecnorm(r-s.*sum(r.*s,1));
I=double(~(behind & distance<RE)); I(rr<=RE)=0;
if strcmp(C.Visibility.EclipseModel,'cylinder'), return; end
for k=1:size(r,2)
    if rr(k)<=RE, continue; end
    toSun=sun(:,k)-r(:,k); ds=norm(toSun); a=asin(min(1,C.Constants.SunRadius_m/ds)); b=asin(RE/rr(k));
    d=acos(max(-1,min(1,dot(-r(:,k)/rr(k),toSun/ds))));
    if d>=a+b, overlap=0;
    elseif d<=abs(b-a), overlap=pi*min(a,b)^2;
    else
        x=max(-1,min(1,(d*d+a*a-b*b)/(2*d*a))); y=max(-1,min(1,(d*d+b*b-a*a)/(2*d*b)));
        overlap=a*a*acos(x)+b*b*acos(y)-.5*sqrt(max(0,(-d+a+b)*(d+a-b)*(d-a+b)*(d+a+b)));
    end
    I(k)=max(0,min(1,1-overlap/(pi*a*a)));
end
end
