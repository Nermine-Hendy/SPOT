function [rows,index]=selectPoolRows(T,S,count,idColumn)
%SELECTPOOLROWS Stable ordered, seeded random (without replacement), or ID selection.
validateattributes(count,{'numeric'},{'scalar','integer','nonnegative'});
index=(1:height(T))';
if ~isempty(S.FilterColumn)
    assert(ismember(S.FilterColumn,T.Properties.VariableNames),'SPOT:Pool','Unknown filter column.');
    index=index(string(T.(S.FilterColumn))==string(S.FilterValue));
end
switch S.Mode
    case 'ordered'
        validateattributes(S.StartRow,{'numeric'},{'scalar','integer','positive'});
        assert(S.StartRow-1+count<=numel(index),'SPOT:Pool','Not enough rows in filtered pool.');
        index=index(S.StartRow:S.StartRow+count-1);
    case 'random'
        validateattributes(S.Seed,{'numeric'},{'scalar','integer','nonnegative','<',2^32});
        assert(count<=numel(index),'SPOT:Pool','Not enough rows for sampling without replacement.');
        rs=RandStream('mt19937ar','Seed',S.Seed); index=index(randperm(rs,numel(index),count));
    case 'ids'
        ids=string(S.IDs(:)); assert(numel(ids)==count && numel(unique(ids))==count,'SPOT:Pool','IDs must be unique and match requested count.');
        [found,where]=ismember(ids,string(T.(idColumn)(index)));
        assert(all(found),'SPOT:Pool','Requested IDs are missing from filtered pool.'); index=index(where);
    otherwise, error('SPOT:Pool','Selection must be ordered, random, or ids.');
end
rows=T(index,:);
end
