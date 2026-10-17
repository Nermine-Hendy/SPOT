function D=detector(C,signal,exposureIndex,active)
%DETECTOR Electron-domain stages followed by read noise, bias and optional ADC.
% Separate RNG streams keep source catalogue changes from changing noise seeds.
seed=mod(C.RandomSeed+700001+7919*exposureIndex,2^32);
rsSignal=RandStream('mt19937ar','Seed',seed);
rsDark=RandStream('mt19937ar','Seed',mod(seed+1000003,2^32));
rsRead=RandStream('mt19937ar','Seed',mod(seed+2000003,2^32));
old=RandStream.getGlobalStream; restore=onCleanup(@() RandStream.setGlobalStream(old));
RandStream.setGlobalStream(rsSignal);
D.ExpectedSignal_e=signal*C.Detector.RadiometricGain_e_per_unit;
if ~isempty(C.Detector.PRNUMap), D.ExpectedSignal_e=D.ExpectedSignal_e.*C.Detector.PRNUMap; end
D.ExpectedDark_e=zeros(size(signal));
if C.Detector.DarkEnabled, D.ExpectedDark_e(:)=C.Detector.DarkCurrent_e_s*C.Time.Exposure_s; end
if C.Detector.PoissonEnabled
    assert(exist('poissrnd','file')==2,'SPOT:Dependency','Shot noise requires Statistics and Machine Learning Toolbox (poissrnd).');
    D.SignalElectrons=poissrnd(D.ExpectedSignal_e);
    RandStream.setGlobalStream(rsDark); D.DarkElectrons=poissrnd(D.ExpectedDark_e);
else
    D.SignalElectrons=D.ExpectedSignal_e; D.DarkElectrons=D.ExpectedDark_e;
end
D.BeforeFullWell_e=D.SignalElectrons+D.DarkElectrons; D.FullWellSaturated=false(size(signal));
D.AfterFullWell_e=D.BeforeFullWell_e;
if C.Detector.FullWellEnabled
    D.FullWellSaturated=D.BeforeFullWell_e>C.Detector.FullWell_e;
    D.AfterFullWell_e=max(0,min(C.Detector.FullWell_e,D.BeforeFullWell_e));
end
D.ReadNoise_e=zeros(size(signal));
if C.Detector.ReadEnabled, D.ReadNoise_e=C.Detector.ReadNoise_e*randn(rsRead,size(signal)); end
D.AfterRead_e=D.AfterFullWell_e+D.ReadNoise_e;
D.Analog_DN=D.AfterRead_e/C.Detector.ConversionGain_e_per_DN;
if C.Detector.BiasEnabled, D.Analog_DN=D.Analog_DN+C.Detector.Bias_DN; end
D.ADCSaturated=false(size(signal)); D.Image_DN=D.Analog_DN;
if C.Detector.ADCEnabled
    hi=2^C.Detector.ADCBits-1; D.ADCSaturated=D.Analog_DN>hi;
    D.Image_DN=round(max(0,min(hi,D.Analog_DN)));
end
% With ADC disabled, retain D0 exactly, including negative read-noise excursions
% (manuscript following Eq.131). Clipping belongs to the enabled ADC stage.
D.ValidPixel=active;
if ~isempty(C.Detector.BadPixelMask), D.ValidPixel=D.ValidPixel & ~logical(C.Detector.BadPixelMask); end
% Bad/unused pixels remain numerically available and explicitly flagged.
clear restore
end
