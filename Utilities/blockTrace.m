function tr = BlockTrace(T1,T2)
% BlockTrace Computes the trace(T1*T2') where T1 and T2 are sparse
% matrices of dimensions (dimT1 x m) and (dimT2 x m)

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 09-Nov-2013 15:14:26

dimT1 = size(T1,1);
tr = 0;
if nargin < 2
    for i = 1:dimT1
        row = T1(i,:);
        tr = tr + sum(nonzeros(T1*row').^2);
    end
else
    dimT2 = size(T2,1);
    if dimT1 <= dimT2
        for i = 1:dimT1
            row = T1(i,:);
            tr = tr + sum(nonzeros(T2*row').^2);
        end
    else
        for i = 1:dimT2
            row = T2(i,:);
            tr = tr + sum(nonzeros(T1*row').^2);
        end
    end
end
end