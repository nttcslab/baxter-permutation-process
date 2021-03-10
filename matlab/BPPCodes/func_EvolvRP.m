function newtrueBlockLocation = func_EvolvRP(BP, LabelMatrix, BetaRandom, trueBlockLocation)


Block = makeBP(LabelMatrix);
ii = size(BP,1)-1;
trueBlockLocation = [trueBlockLocation; zeros(1,4)];
if Block(ii, 2) < Block(ii+1, 2)
    aboveBlockIndex = LabelMatrix(Block(ii+1, 1)-1, Block(ii+1,3):Block(ii+1,4) );
    temp = trueBlockLocation(aboveBlockIndex, :);
    minLength = min(temp(:, 2)-temp(:, 1));
    cutLength = (1-BetaRandom) * minLength;
    trueBlockLocation(aboveBlockIndex, 2) = trueBlockLocation(aboveBlockIndex, 2) - cutLength;
    trueBlockLocation(ii+1, 1) = max( trueBlockLocation(aboveBlockIndex, 2) );
    trueBlockLocation(ii+1, 2) = max( trueBlockLocation(aboveBlockIndex, 2) ) + cutLength;
    trueBlockLocation(ii+1, 3) = min( trueBlockLocation(aboveBlockIndex, 3) );
    trueBlockLocation(ii+1, 4) = max( trueBlockLocation(aboveBlockIndex, 4) );
else
    leftBlockIndex = LabelMatrix(Block(ii+1, 1):Block(ii+1, 2), Block(ii+1,3)-1);
    temp = trueBlockLocation(leftBlockIndex, :);
    minLength = min(temp(:, 4)-temp(:, 3));
    cutLength = (1-BetaRandom) * minLength;
    trueBlockLocation(leftBlockIndex, 4) = trueBlockLocation(leftBlockIndex, 4) - cutLength;
    trueBlockLocation(ii+1, 1) = min( trueBlockLocation(leftBlockIndex, 1) );
    trueBlockLocation(ii+1, 2) = max( trueBlockLocation(leftBlockIndex, 2) );
    trueBlockLocation(ii+1, 3) = min( trueBlockLocation(leftBlockIndex, 4) );
    trueBlockLocation(ii+1, 4) = min( trueBlockLocation(leftBlockIndex, 4) )+ cutLength;
end
newtrueBlockLocation = trueBlockLocation;