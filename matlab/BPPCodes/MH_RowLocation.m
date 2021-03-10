function [RowLocation, ColumnLocation]...
    = MH_RowLocation(X, trueBlockLocation, RowLocation, ColumnLocation, alpha)

for ii=1:size(RowLocation,1)
    EachClusterCount = makeEachClusterCount(X, trueBlockLocation, RowLocation, ColumnLocation);
    currentProb = calc_DirichletLikelihood(EachClusterCount, alpha);
    
    currentRowLocation = RowLocation;
    RowLocation(ii) = rand(1,1);
    EachClusterCount = makeEachClusterCount(X, trueBlockLocation, RowLocation, ColumnLocation);
    newProb = calc_DirichletLikelihood(EachClusterCount, alpha);
    
    AcceptReject = rand(1,1);
    if log(AcceptReject+eps) > newProb - currentProb
        RowLocation = currentRowLocation;
    end
end

for ii=1:size(ColumnLocation,1)
    EachClusterCount = makeEachClusterCount(X, trueBlockLocation, RowLocation, ColumnLocation);
    currentProb = calc_DirichletLikelihood(EachClusterCount, alpha);
    
    currentColumnLocation = ColumnLocation;
    ColumnLocation(ii) = rand(1,1);
    EachClusterCount = makeEachClusterCount(X, trueBlockLocation, RowLocation, ColumnLocation);
    newProb = calc_DirichletLikelihood(EachClusterCount, alpha);
    
    AcceptReject = rand(1,1);
    if log(AcceptReject+eps) > newProb - currentProb
        ColumnLocation = currentColumnLocation;
    end
end