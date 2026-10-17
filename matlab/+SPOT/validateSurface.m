function validateSurface(T)
%VALIDATESURFACE Reject invalid materials, geometry and target rotations.
S=T.Surface; nf=size(S.NormalBody,2);
assert(size(S.NormalBody,1)==3 && nf>=1 && all(isfinite(S.NormalBody(:))) && all(abs(vecnorm(S.NormalBody)-1)<1e-8),'SPOT:Surface','Normals must be unit body-frame 3-by-N vectors.');
for name={'Area_m2','DiffuseReflectance','SpecularReflectance','SpecularExponent'}
    a=S.(name{1}); assert(isrow(a) && numel(a)==nf && all(isfinite(a)) && all(a>=0),'SPOT:Surface','Each surface field must contain one finite nonnegative value per face.');
end
assert(all(S.Area_m2>0) && all(S.DiffuseReflectance+S.SpecularReflectance<=1),'SPOT:Surface','Positive areas and diffuse+specular <=1 are required.');
assert(all(S.SpecularExponent>=1 & S.SpecularExponent<=1e5),'SPOT:Surface','Specular exponent must be in [1,1e5]; ideal mirrors are unsupported.');
A=T.Attitude;
assert(any(strcmp(A.Mode,{'fixed','spin','samples'})),'SPOT:TargetAttitude','Invalid target attitude mode.');
assert(isdatetime(A.EpochUTC) && isscalar(A.EpochUTC) && ~isnat(A.EpochUTC) && ~isempty(A.EpochUTC.TimeZone),'SPOT:TargetAttitude','Epoch must be timezone-aware.');
rotation(A.BodyToGCRF0);
assert(numel(A.SpinAxisBody)==3 && all(isfinite(A.SpinAxisBody)) && norm(A.SpinAxisBody)>0 && isscalar(A.SpinRate_deg_s) && isfinite(A.SpinRate_deg_s),'SPOT:TargetAttitude','Invalid spin axis/rate.');
if strcmp(A.Mode,'samples')
    t=A.TimeUTC;
    assert(isdatetime(t) && ~isempty(t.TimeZone) && numel(t)>=2 && all(~isnat(t)) && all(seconds(diff(t(:)))>0) && size(A.BodyToGCRF,3)==numel(t),'SPOT:TargetAttitude','Attitude samples must have increasing UTC times.');
    for k=1:numel(t)
        rotation(A.BodyToGCRF(:,:,k));
        if k>1
            q=(trace(A.BodyToGCRF(:,:,k-1)'*A.BodyToGCRF(:,:,k))-1)/2;
            assert(q>cosd(179),'SPOT:TargetAttitude','Add attitude knots: successive rotations must be <179 degrees.');
        end
    end
end
end
function rotation(R)
assert(isequal(size(R),[3 3]) && all(isfinite(R(:))) && norm(R'*R-eye(3),'fro')<1e-8 && det(R)>0,'SPOT:TargetAttitude','BodyToGCRF must be a proper rotation.');
end
