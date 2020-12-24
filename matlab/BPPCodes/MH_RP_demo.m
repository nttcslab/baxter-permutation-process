function [BP, BetaRandom, trueBlockLocation, UniformRandom]...
    = MH_RP_demo(X, BP, BetaRandom, trueBlockLocation, UniformRandom, RowLocation, ColumnLocation, alpha, ENC)

currentBP = BP;
currentBetaRandom = BetaRandom;
currenttrueBlockLocation = trueBlockLocation;
currentUniformRandom = UniformRandom;
EachClusterCount = makeEachClusterCount(X, trueBlockLocation, RowLocation, ColumnLocation);
currentProb = calc_DirichletLikelihood(EachClusterCount, alpha)...
    - (size(UniformRandom,1)-ENC*log(size(UniformRandom,1)));

if rand(1)>0.5
    [BP, UniformRandom] = func_EvolvBP(BP, UniformRandom);
    LabelMatrix = func_BP2FP(BP);
    BetaRandom = [BetaRandom; betarnd(1,1)];
    trueBlockLocation = func_EvolvRP(BP, LabelMatrix, BetaRandom(end), trueBlockLocation);
else
    if size(BP,1)>1
        LabelMatrix = func_BP2FP(BP);
        trueBlockLocation = func_deleteRP(BP, LabelMatrix, trueBlockLocation);
        [BP, UniformRandom] = func_DeleteMax(BP, UniformRandom);
        BetaRandom(end) = [];
    end
end
EachClusterCount = makeEachClusterCount(X, trueBlockLocation, RowLocation, ColumnLocation);
newProb = calc_DirichletLikelihood(EachClusterCount, alpha)...
    - (size(UniformRandom,1)-ENC*log(size(UniformRandom,1)));

AcceptReject = rand(1,1);
if log(AcceptReject+eps) > newProb - currentProb
    BP = currentBP;
    BetaRandom = currentBetaRandom;
    trueBlockLocation = currenttrueBlockLocation;
    UniformRandom = currentUniformRandom;
end

    
