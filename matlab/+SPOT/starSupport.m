function keep=starSupport(C,directions,cameraToGCRF)
%STARSUPPORT Conservative union of cones enclosing every sampled camera FoV.
% This only limits working memory. Exact FoV/visibility tests still follow.
% The rectangle's corner angle bounds all rays in the rectangular FoV.
halfAngle=atand(norm(tand(C.Camera.FOV_deg/2)));
if strcmp(C.Camera.Shape,'circular'), halfAngle=C.Camera.FOV_deg(1)/2; end
limit=cosd(halfAngle)-1e-10; keep=false(1,size(directions,2));
for k=1:size(cameraToGCRF,3)
    keep=keep | cameraToGCRF(:,3,k)'*directions>=limit;
end
end
