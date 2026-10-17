function H=photometry(C,G,j)
%PHOTOMETRY Surface prediction or explicit legacy photometry, then atmosphere.
H=SPOT.intrinsicPhotometry(C,G,j);
[H.AtmosphericTransmission,H.ScintillationGain]=SPOT.atmosphere(C,G.TimeUTC,H.Elevation_deg,G.IsGround,C.Targets(j).ID);
H.RadiometricRate=H.IntrinsicRate.*H.AtmosphericTransmission.*H.ScintillationGain;
H.ReceivedIrradiance_W_m2=H.RadiometricRate*C.Radiometry.ZeroMagnitudeIrradiance_W_m2;
H.ApparentMagnitude=inf(size(H.RadiometricRate)); positive=H.RadiometricRate>0;
H.ApparentMagnitude(positive)=-2.5*log10(H.RadiometricRate(positive));
assert(all(isfinite(H.RadiometricRate)) && all(H.RadiometricRate>=0),'SPOT:Photometry','Invalid target radiometric rate.');
end
