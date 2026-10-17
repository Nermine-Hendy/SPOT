%RUN_SPOT_ALLTARGETGLINTS Moderate reflection peaks for all selected targets.
% Requires C and E from an example's editable configuration section.
% Run this script after defining those settings and before calling RUN_SPOT.
% This script applies the overrides, configures glints and runs the acquisition.
%
% These settings define a controlled glint demonstration, not calibrated
% material properties: each eligible target's +z face is aligned to reflect
% sunlight toward the observer at the chosen time in the FIRST exposure.
% Diffuse coefficients, face areas, target trajectories, exposure and camera
% settings are retained. Different targets need not have equal irradiance.

%% SETTINGS: edit these values
maxSpecularToDiffuse = 0.40; % maximum instantaneous specular/diffuse ratio
specularExponent = 40;      % specular lobe shape
spinRate_deg_s = 10;        % rotation speed; controls the glint duration
glintTimeFraction = 0.50;   % 0=start, 0.5=middle, 1=end of first exposure

%% APPLY AND RUN
assert(exist('C','var') == 1 && exist('E','var') == 1, ...
    'SPOT:GlintExample', 'Run the example settings first to define C and E.');
validateattributes(maxSpecularToDiffuse, {'numeric'}, {'scalar','finite','nonnegative','<=',1});
validateattributes(specularExponent, {'numeric'}, {'scalar','finite','>=',1,'<=',1e5});
validateattributes(spinRate_deg_s, {'numeric'}, {'scalar','finite'});
validateattributes(glintTimeFraction, {'numeric'}, {'scalar','finite','>=',0,'<=',1});

C = SPOT_ApplyExample(C,E);
C = SPOT.validate(C);
% Inputs are now resolved, including any selected CSV rows and seed overrides.
% Keep those exact states/rows when RUN_SPOT validates the edited configuration.
% SelectedSensor/SelectedTargets and per-target source metadata stay recorded.
C.Input.UseSensorCSV = false;
C.Input.UseTargetCSV = false;
C = rmfield(C, intersect(fieldnames(C), {'ExampleTarget','ExamplePointing'}));

tGlint = C.Time.StartUTC + seconds(glintTimeFraction*C.Time.Exposure_s);
G = SPOT.sample(C,tGlint,tGlint);

for j = 1:numel(C.Targets)
    assert(strcmp(C.Targets(j).Photometry.Model, 'surface_brdf'), ...
        'SPOT:GlintExample', 'Target %d must use surface_brdf photometry.', C.Targets(j).ID);

    % Select a modest specular coefficient relative to EACH face's existing
    % diffuse reflectance. This bounds the actual model, not just its plot:
    % fs/fd <= rho_s*(p+2)/(2*rho_d) <= maxSpecularToDiffuse.
    rd = C.Targets(j).Surface.DiffuseReflectance;
    C.Targets(j).Surface.SpecularExponent(:) = specularExponent;
    rs = 2*maxSpecularToDiffuse*rd/(specularExponent+2);
    C.Targets(j).Surface.SpecularReflectance = min(rs, max(0,1-rd));

    H = G.Target{j};
    if ~H.Visible || H.Illumination <= 0
        warning('SPOT:GlintNotVisible', ...
            'Target %d is not visible and sunlit at the glint time; its attitude is retained.', ...
            C.Targets(j).ID);
        continue
    end
    normals = C.Targets(j).Surface.NormalBody;
    assert(any(vecnorm(normals-[0;0;1]) < 1e-8), 'SPOT:GlintExample', ...
        'Target %d needs a +z face; this example uses the default six-face box.', C.Targets(j).ID);

    % Directions from this target to the Sun and observer.
    s = G.SunPosition_m-H.Position_m;
    s = s/norm(s);
    v = -H.LOS_GCRF;
    n = s+v;
    if norm(n) < 1e-10
        warning('SPOT:GlintGeometry', ...
            'Target %d has degenerate reflection geometry; its attitude is retained.', C.Targets(j).ID);
        continue
    end
    n = n/norm(n);

    ref = [0;0;1];
    if abs(dot(ref,n)) > 0.95, ref = [0;1;0]; end
    bx = cross(ref,n); bx = bx/norm(bx);
    by = cross(n,bx);

    C.Targets(j).Attitude.Mode = 'spin';
    C.Targets(j).Attitude.EpochUTC = tGlint;
    C.Targets(j).Attitude.BodyToGCRF0 = [bx by n];
    C.Targets(j).Attitude.SpinAxisBody = [0;1;0];
    C.Targets(j).Attitude.SpinRate_deg_s = spinRate_deg_s;
    C.Targets(j).Attitude.TimeUTC = datetime.empty(0,1);
    C.Targets(j).Attitude.BodyToGCRF = zeros(3,3,0);
end

C.Plot.Photometry = true;
Out = RUN_SPOT(C);

% Report the largest instantaneous ratio in the recorded first exposure.
% The same coefficient bound applies to subsequent exposures as well.
for j = 1:numel(Out.Exposures{1}.Targets)
    B = Out.Exposures{1}.Targets(j).Samples.Reflection;
    use = B.DiffuseIrradiance_W_m2 > 0;
    if any(use)
        ratio = max(B.SpecularIrradiance_W_m2(use)./B.DiffuseIrradiance_W_m2(use));
        fprintf('Target %d: maximum recorded specular/diffuse = %.3f (limit %.3f)\n', ...
            Out.Exposures{1}.Targets(j).ID, ratio, maxSpecularToDiffuse);
    else
        fprintf('Target %d: no positive diffuse samples in this exposure.\n', ...
            Out.Exposures{1}.Targets(j).ID);
    end
end