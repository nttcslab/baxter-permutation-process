function [BP, LabelMatrix, UniformRandom] = FPupdate(X, CRPv, CRPh, alpha, BP, LabelMatrix, UniformRandom)
currentBP = BP;
currentLabelMatrix = LabelMatrix;
currentUniformRandom = UniformRandom;
currentProb = calc_DirichletLikelihood(X, LabelMatrix(CRPv, CRPh), alpha);

n = size(BP,1);

if rand(1)>0.5
    [BP, UniformRandom] = func_EvolvBP(BP, UniformRandom);
else
    if size(BP,1)>1
        [BP, UniformRandom] = func_DeleteMax(BP, UniformRandom);
        CRPv(CRPv==n) = randi(n-1, size(CRPv(CRPv==n)));
        CRPh(CRPh==n) = randi(n-1, size(CRPh(CRPh==n)));
    end
end
LabelMatrix = func_BP2FP(BP);
newProb = calc_DirichletLikelihood(X, LabelMatrix(CRPv, CRPh), alpha);

AcceptReject = rand(1,1);
if log(AcceptReject+eps) > newProb - currentProb
    BP = currentBP;
    LabelMatrix = currentLabelMatrix;
    UniformRandom = currentUniformRandom;
end