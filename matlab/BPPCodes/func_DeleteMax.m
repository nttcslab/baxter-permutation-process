

function [newBP, updateUniformRandom] = func_DeleteMax(BP, UniformRandom)

[A, I] = max(BP);
newBP = BP;
newBP(I) = [];
updateUniformRandom = UniformRandom(1:end-1);