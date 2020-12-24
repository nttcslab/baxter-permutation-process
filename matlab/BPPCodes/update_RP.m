function [BP, BetaRandom, trueBlockLocation, UniformRandom]...
    = update_RP(BetaRandom, UniformRandom, X, RowLocation, ColumnLocation, alpha)

[BP, trueBlockLocation] = reconstructRP(BetaRandom, UniformRandom);

for ii=1:size(UniformRandom,1)
    currentBP = BP;
    currentBetaRandom = BetaRandom;
    currenttrueBlockLocation = trueBlockLocation;
    currentUniformRandom = UniformRandom;
    EachClusterCount = makeEachClusterCount(X, trueBlockLocation, RowLocation, ColumnLocation);
    currentProb = calc_DirichletLikelihood(EachClusterCount, alpha);
    
    UniformRandom(ii) = rand;
    [BP, trueBlockLocation] = reconstructRP(BetaRandom, UniformRandom);
    EachClusterCount = makeEachClusterCount(X, trueBlockLocation, RowLocation, ColumnLocation);
    newProb = calc_DirichletLikelihood(EachClusterCount, alpha);
    
    AcceptReject = rand(1,1);
    if log(AcceptReject+eps) > newProb - currentProb
        BP = currentBP;
        BetaRandom = currentBetaRandom;
        trueBlockLocation = currenttrueBlockLocation;
        UniformRandom = currentUniformRandom;
    end
end

for ii=1:size(BetaRandom,1)
    currentBP = BP;
    currentBetaRandom = BetaRandom;
    currenttrueBlockLocation = trueBlockLocation;
    currentUniformRandom = UniformRandom;
    EachClusterCount = makeEachClusterCount(X, trueBlockLocation, RowLocation, ColumnLocation);
    currentProb = calc_DirichletLikelihood(EachClusterCount, alpha);
    
    BetaRandom(ii) = betarnd(1,1);
    [BP, trueBlockLocation] = reconstructRP(BetaRandom, UniformRandom);
    EachClusterCount = makeEachClusterCount(X, trueBlockLocation, RowLocation, ColumnLocation);
    newProb = calc_DirichletLikelihood(EachClusterCount, alpha);
    
    AcceptReject = rand(1,1);
    if log(AcceptReject+eps) > newProb - currentProb
        BP = currentBP;
        BetaRandom = currentBetaRandom;
        trueBlockLocation = currenttrueBlockLocation;
        UniformRandom = currentUniformRandom;
    end
end