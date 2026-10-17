%DIAGNOSIS Print peak diffuse and specular irradiance from Out.
% Run a single-acquisition example before invoking this script.
for j = 1:numel(Out.Exposures{1}.Targets)
    T = Out.Exposures{1}.Targets(j);
    B = T.Samples.Reflection;
    fprintf('%s: peak diffuse %.3e, peak specular %.3e W/m^2\n', ...
        T.Name, max(B.DiffuseIrradiance_W_m2), ...
        max(B.SpecularIrradiance_W_m2));
end