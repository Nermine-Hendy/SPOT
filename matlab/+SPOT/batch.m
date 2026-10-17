function index=batch(configurations)
%SPOT_BATCH Run a cell array of complete unified configurations; fail clearly.
assert(iscell(configurations),'SPOT:Batch','Pass a cell array of SPOT_Config structures.');
rows=cell(numel(configurations),4);
for k=1:numel(configurations)
    C=configurations{k}; Out=RUN_SPOT(C);
    rows(k,:)={string(C.Dataset.SampleID),string(C.Dataset.Split),exportDirectory(Out),numel(Out.Exposures)};
end
index=cell2table(rows,'VariableNames',{'sample_id','split','directory','num_exposures'});
end

function folder=exportDirectory(Out)
folder=""; if isfield(Out,'ExportDirectory'), folder=string(Out.ExportDirectory); end
end
