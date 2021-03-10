function newtrueBlockLocation = func_deleteRP(BP, LabelMatrix, trueBlockLocation)

Block = makeBP(LabelMatrix);

ii = size(BP,1)-1;
if Block(ii, 2) < Block(ii+1, 2)
    aboveBlockIndex = LabelMatrix(Block(ii+1, 1)-1, Block(ii+1,3):Block(ii+1,4) );
    cutLength = trueBlockLocation(ii+1, 2) - trueBlockLocation(ii+1, 1);
    trueBlockLocation(aboveBlockIndex, 2) = trueBlockLocation(aboveBlockIndex, 2) + cutLength;
else
    leftBlockIndex = LabelMatrix(Block(ii+1, 1):Block(ii+1, 2), Block(ii+1,3)-1);
    cutLength = trueBlockLocation(ii+1, 4) - trueBlockLocation(ii+1, 3);
    trueBlockLocation(leftBlockIndex, 4) = trueBlockLocation(leftBlockIndex, 4) + cutLength;
end
newtrueBlockLocation = trueBlockLocation;
newtrueBlockLocation(ii+1, :) = [];