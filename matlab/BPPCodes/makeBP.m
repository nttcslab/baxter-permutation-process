function Block = makeBP(LabelMatrix)

Block = zeros(size(max(LabelMatrix(:)),1), 4);

temp = LabelMatrix;
for ii=1:max(LabelMatrix(:))
    temptemp=temp;
    temptemp(temptemp~=ii)=0;
    minR = find( sum(temptemp,2)>0 ,1);
    maxR = find(sum(temptemp,2)>0, 1, 'last');
    minC = find(sum(temptemp,1)>0, 1);
    maxC = find(sum(temptemp,1)>0, 1, 'last');
    Block(ii,:) = [minR maxR minC maxC];
end