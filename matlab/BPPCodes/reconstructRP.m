function [BP, trueBlockLocation] = reconstructRP(BetaRandom, UniformRandom)

BP = [1];
trueBlockLocation = zeros(1, 4);
trueBlockLocation(1,:) = [0 1 0 1];
for ii=1:size(UniformRandom,1)-1
    BP = temp_EvolvBP(BP, UniformRandom(1:ii), UniformRandom(ii+1));
    LabelMatrix = func_BP2FP(BP);
    trueBlockLocation = func_EvolvRP(BP, LabelMatrix, BetaRandom(ii), trueBlockLocation);
end

function newBP = temp_EvolvBP(BP, UniformRandom, NewUniformRandom)

n = size(BP,1);

[Left2RightMaxima, Right2LeftMaxima, Left2RightLocation, Right2LeftLocation]...
    = ExtractMaxima(BP);

if NewUniformRandom < UniformRandom(Left2RightMaxima(1))
    BP = [BP(1:Left2RightLocation(1)-1); n+1; BP(Left2RightLocation(1):end)];
end
if size(Left2RightMaxima,1)>1
    for ii=2:size(Left2RightMaxima,1)
        if UniformRandom(Left2RightMaxima(ii-1)) <= NewUniformRandom...
                && NewUniformRandom < UniformRandom(Left2RightMaxima(ii))
            BP = [BP(1:Left2RightLocation(ii)-1); n+1; BP(Left2RightLocation(ii):end)];
        end
    end
end
if UniformRandom(Right2LeftMaxima(1)) < NewUniformRandom
    BP = [BP(1:Right2LeftLocation(1)); n+1; BP(Right2LeftLocation(1)+1:end)];
end
if size(Right2LeftMaxima,1)>1
    for ii=2:size(Right2LeftMaxima,1)
        if  NewUniformRandom <= UniformRandom(Right2LeftMaxima(ii-1))...
                && UniformRandom(Right2LeftMaxima(ii)) < NewUniformRandom
            BP = [BP(1:Right2LeftLocation(ii)); n+1; BP(Right2LeftLocation(ii)+1:end)];
        end
    end
end

newBP = BP;
updateUniformRandom = [UniformRandom; NewUniformRandom];