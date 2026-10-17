function d=SPOT_ConfigSeed(mode,row)
d=struct('AngleMode',mode,'Angles_deg',row(1:2),'Range_m',row(3), ...
    'TransverseSpeed_mps',row(4),'Direction_deg',row(5),'RadialSpeed_mps',row(6));
end

