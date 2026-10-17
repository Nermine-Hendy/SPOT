function C=SPOT_ApplyExample(C,E)
%SPOT_APPLYEXAMPLE Apply the editable section of an example to a configuration.
% Target overrides are applied AFTER CSV selection, inside SPOT.validate.
required={'Exposure_s','Cadence_s','NumExposures','SaveResults','SavePlots', ...
    'ShowPlots','PlotExposures','PlotProducts','PlotPhotometry','Target'};
assert(isstruct(E) && all(isfield(E,required)),'SPOT:Example','Use the complete settings block from an example.');
C.Time.Exposure_s=E.Exposure_s; C.Time.Cadence_s=E.Cadence_s;
C.Time.NumExposures=E.NumExposures;
for name={'SaveResults','SavePlots','ShowPlots','PlotPhotometry'}
    validateattributes(E.(name{1}),{'logical','numeric'},{'scalar','binary'});
end
C.Output.SaveResults=logical(E.SaveResults); C.Output.Plot=logical(E.ShowPlots);
C.Plot.Save=logical(E.SavePlots); C.Plot.Exposures=E.PlotExposures;
C.Plot.Photometry=logical(E.PlotPhotometry);
names={'camera','detector','radec','stages','dataset'};
fields={'Camera','Detector','RADEC','Stages','Dataset'};
products=E.PlotProducts;
if isstring(products) || ischar(products), products=cellstr(products); end
assert(iscell(products) && all(ismember(products,names)), ...
    'SPOT:Example','PlotProducts supports camera, detector, radec, stages and dataset.');
for j=1:numel(names), C.Plot.(fields{j})=any(strcmp(products,names{j})); end
if isfield(E,'Pointing')
    names=fieldnames(E.Pointing);
    allowed={'AzEl_deg','Roll_deg','OffNadir_deg','SpaceAzimuth_deg'};
    assert(all(ismember(names,allowed)),'SPOT:Example','Unsupported pointing override.');
    % Apply after CSV sensor selection so example-level edits always win.
    C.ExamplePointing=E.Pointing;
end
C.ExampleTarget=E.Target;
end
