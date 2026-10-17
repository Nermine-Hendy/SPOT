function R=convergence(C)
%SPOT_CONVERGENCE Compare two noiseless runs with time and spatial steps halved.
% Use on final scenarios: visibility discontinuities can require much finer grids.
C.Detector.PoissonEnabled=false; C.Detector.ReadEnabled=false; C.Output.Verbose=false; C.Output.Plot=false;
A=SPOT.run(C,false); C.Integration.Step_s=C.Integration.Step_s/2;
C.Integration.MaxStep_pix=C.Integration.MaxStep_pix/2;
C.Integration.PhotometryRelTol=C.Integration.PhotometryRelTol/2;
C.Integration.PhotometryAbsTol=C.Integration.PhotometryAbsTol/2;
C.Integration.MaxBodyRotation_deg=C.Integration.MaxBodyRotation_deg/2;
C.Integration.SpecularSamplesPerHWHM=2*C.Integration.SpecularSamplesPerHWHM;
C.Integration.MaxSamples=2*C.Integration.MaxSamples;
C.Integration.MaxBackgroundSamples=2*C.Integration.MaxBackgroundSamples;
B=SPOT.run(C,false); R=struct('RelativeL1',zeros(1,numel(A.Exposures)),'TargetFluxRelativeChange',{{}},'MaskIoU',zeros(1,numel(A.Exposures)));
for k=1:numel(A.Exposures)
    a=A.Exposures{k}; b=B.Exposures{k};
    R.RelativeL1(k)=sum(abs(a.OpticalStages.PSF-b.OpticalStages.PSF),'all')/max(sum(abs(b.OpticalStages.PSF),'all'),realmin);
    R.TargetFluxRelativeChange{k}=abs([a.Targets.IntegratedOpticalSignal]-[b.Targets.IntegratedOpticalSignal])./max([b.Targets.IntegratedOpticalSignal],realmin);
    R.MaskIoU(k)=nnz(a.Labels.Binary & b.Labels.Binary)/max(nnz(a.Labels.Binary | b.Labels.Binary),1);
    if ~any(a.Labels.Binary | b.Labels.Binary,'all'), R.MaskIoU(k)=1; end
end
R.Note='No universal tolerance is imposed. Assess image and target flux differences and mask IoU against your scientific use case.';
end
