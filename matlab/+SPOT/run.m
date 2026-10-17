function Out=run(C,exportProducts)
%RUN End-to-end SPOT acquisition, with one configuration and one pipeline.
% Out=SPOT.run(SPOT_Config('ground'));  % also exports files
% Out=SPOT.run(C,false);                     % compute in memory only
if nargin==0, C=SPOT_Config; end
if ischar(C) || isstring(C), C=SPOT_Config(C); end
if nargin<2, exportProducts=true; end
C=SPOT.validate(C); stars=SPOT.loadStars(C);
Out.SchemaVersion=C.SchemaVersion; Out.Configuration=SPOT.serializable(C,false);
Out.Provenance=struct('MATLAB',version,'Computer',computer,'CreatedUTC',datetime('now','TimeZone','UTC'), ...
    'Radiometry','target surface irradiance referenced to effective-band zero point; backgrounds remain configured phenomenological radiances', ...
    'DetectorCalibration',C.Radiometry.DetectorCalibration,'RadiometricGain_e_per_unit',C.Detector.RadiometricGain_e_per_unit, ...
    'Frame','GCRF','EarthOrientation',C.Frames.Model,'RandomSeed',C.RandomSeed);
Out.Exposures=cell(1,C.Time.NumExposures);
for ie=1:C.Time.NumExposures
    start=C.Time.StartUTC+seconds((ie-1)*C.Time.Cadence_s); mid=start+seconds(C.Time.Exposure_s/2);
    [G,w,requested]=SPOT.integrateGeometry(C,start); B=SPOT.background(C,start,stars);
    shape=[C.Camera.Ny C.Camera.Nx]; targetRate=zeros(shape); maps=cell(1,numel(C.Targets));
    template=struct('ID',0,'Name','','Class','','InputTrajectory','','Appearance','', 'Morphology',struct(), ...
        'DisplacementBound_pix',0,'RelativeFluxChange',0, ...
        'MeanReceivedIrradiance_W_m2',0,'MeanApparentMagnitude',Inf, ...
        'AnyVisible',false,'VisibleFraction',0,'IlluminatedFraction',0,'MeanVisibleRate',0, ...
        'InputDepositedRate',0,'IntegratedOpticalSignal',0,'Samples',struct(), ...
        'HasAnnotation',false,'PeakOpticalSignal',0,'MaskPixelCount',0,'MaskThreshold',0, ...
        'BoundingBox_xywh',[NaN NaN 0 0],'BoxCenter_xy',[NaN NaN]);
    records=repmat(template,1,numel(C.Targets));
    for j=1:numel(C.Targets)
        H=SPOT.photometry(C,G,j); H.RA_deg=mod(atan2d(H.LOS_GCRF(2,:),H.LOS_GCRF(1,:)),360);
        H.DEC_deg=asind(max(-1,min(1,H.LOS_GCRF(3,:))));
        displacement=hypot(max(H.x)-min(H.x),max(H.y)-min(H.y));
        variation=(max(H.RadiometricRate)-min(H.RadiometricRate))/max(max(H.RadiometricRate),realmin);
        M=SPOT.morphology(C,H,w); kind=M.Kind;
        % One renderer, including compact targets, partial visibility and light curves.
        rate=SPOT.deposit(C,H.x,H.y,w.*H.RadiometricRate.*H.Visible,B.Truth.Active);
        maps{j}=SPOT.optics(C,rate,B.Truth.Active); targetRate=targetRate+rate;
        records(j).ID=C.Targets(j).ID; records(j).Name=C.Targets(j).Name; records(j).Class=C.Targets(j).Class;
        records(j).InputTrajectory=C.Targets(j).Trajectory.Type;
        records(j).Appearance=kind; records(j).Morphology=M;
        records(j).DisplacementBound_pix=displacement;
        records(j).RelativeFluxChange=variation; records(j).AnyVisible=any(H.Visible);
        records(j).VisibleFraction=sum(w.*H.Visible); records(j).IlluminatedFraction=sum(w.*H.Illumination);
        records(j).MeanVisibleRate=sum(w.*H.Visible.*H.RadiometricRate);
        records(j).MeanReceivedIrradiance_W_m2=records(j).MeanVisibleRate*C.Radiometry.ZeroMagnitudeIrradiance_W_m2;
        if records(j).MeanVisibleRate>0, records(j).MeanApparentMagnitude=-2.5*log10(records(j).MeanVisibleRate); end
        records(j).InputDepositedRate=sum(rate,'all'); records(j).IntegratedOpticalSignal=sum(maps{j},'all');
        H.TimeUTC=G.TimeUTC; H.Weights=w;
        % Full samples preserve gaps; this filtered centerline must not be joined across gaps.
        H.VisibleSampleIndex=find(H.Visible); records(j).Samples=H;
    end
    scene=B.StarRate+B.DiffuseRate+targetRate;
    [signal,optical]=SPOT.optics(C,scene,B.Truth.Active);
    detector=SPOT.detector(C,signal,ie,B.Truth.Active);
    [labels,records]=SPOT.labels(C,maps,records);
    E=struct('Index',ie,'StartUTC',start,'MidUTC',mid,'EndUTC',start+seconds(C.Time.Exposure_s), ...
        'RequestedTimeUTC',requested,'TargetTimeUTC',G.TimeUTC,'TargetWeights',w, ...
        'BackgroundTimeUTC',B.TimeUTC,'BackgroundWeights',B.Weights, ...
        'Observer',G.Observer,'CameraToGCRF',G.CameraToGCRF,'Truth',B.Truth, ...
        'StarTruth',B.StarTruth,'SkyRate',B.SkyRate,'EarthRate',B.EarthRate, ...
        'StarRate',B.StarRate,'DiffuseRate',B.DiffuseRate,'TargetRate',targetRate, ...
        'SceneRate',scene,'OpticalStages',optical,'DetectorStages',detector, ...
        'Image_DN',detector.Image_DN,'Labels',labels,'Targets',records, ...
        'StarMetadata',B.StarCatalogueMetadata,'StarsUsed',B.StarsUsed, ...
        'RequestedRegime',C.Dataset.RequestedRegime);
    E.TotalTargetOpticalSignal=zeros(shape);
    for j=1:numel(maps), E.TotalTargetOpticalSignal=E.TotalTargetOpticalSignal+maps{j}; end
    if C.Output.KeepTargetMaps, E.TargetOpticalSignal=maps; else, E.TargetOpticalSignal={}; end
    E.Sampling=struct('Adaptive',C.Integration.Adaptive,'MaxStep_pix',C.Integration.MaxStep_pix, ...
        'TargetSampleCount',numel(G.TimeUTC),'BackgroundSampleCount',numel(B.TimeUTC), ...
        'PhotometryAdaptive',C.Integration.PhotometryAdaptive,'PhotometryRelTol',C.Integration.PhotometryRelTol, ...
        'TimeConvergenceCertified',false,'Note','Spatial and midpoint intrinsic-flux refinement with spin/lobe step cap; run SPOT.convergence for final scientific scenarios.');
    Out.Exposures{ie}=E;
    if C.Output.Verbose, fprintf('SPOT exposure %d/%d: %s, %d target samples, %d annotated targets.\n',ie,C.Time.NumExposures,E.Truth.GeometryRegime,numel(G.TimeUTC),nnz([records.HasAnnotation])); end
end
if exportProducts, Out.ExportDirectory=SPOT.exportResults(Out,C); end
end
