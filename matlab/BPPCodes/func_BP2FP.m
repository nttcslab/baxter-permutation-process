function LabelMatrix = func_BP2FP(BP)

n = size(BP,1);
LabelMatrix = zeros( size(BP,1), size(BP,1) );

Block = zeros(n, 4);
Block(BP(1),:) = [1 n 1 n];
LabelMatrix(:,:) = BP(1);

for ii=2:size(BP,1)
    if BP(ii-1)<BP(ii)
        temp = Block( BP(ii-1), : );
        Block( BP(ii), : ) = [temp(1) temp(2) ii temp(4)];
        LabelMatrix( temp(1):temp(2), ii:temp(4) ) = BP(ii);
        Block( BP(ii-1), : ) = [temp(1) temp(2) temp(3) ii-1];
        temp = Block( BP(ii), : );
        rowindex = Block( BP(ii), 2 );
        while rowindex < n
            if LabelMatrix( rowindex+1, temp(4) ) < LabelMatrix( rowindex, temp(4) )
                Block( BP(ii), 2 ) = rowindex + 1;
                LabelMatrix( rowindex+1, temp(3):temp(4) ) = LabelMatrix( rowindex, temp(4) );
                rowindex = rowindex + 1;
            else
                rowindex = n;
            end
        end
    else
        temp = Block( BP(ii-1), : );
        Block( BP(ii), : ) = [temp(1) n-ii+1 temp(3) temp(4)];
        LabelMatrix( temp(1):(n-ii+1), temp(3):temp(4) ) = BP(ii);
        Block( BP(ii-1), : ) = [n-ii+2 temp(2) temp(3) temp(4)];
        temp = Block( BP(ii), : );
        columnindex = Block( BP(ii), 3 );
        while columnindex > 1
            if LabelMatrix( temp(1), columnindex-1 ) > LabelMatrix( temp(1), columnindex ) 
                Block( BP(ii), 3 ) = columnindex - 1;
                LabelMatrix( temp(1):temp(2), columnindex-1 ) = LabelMatrix( temp(1), columnindex );
                columnindex = columnindex - 1;
            else
                columnindex = 1;
            end
        end
    end
end


for ii=1:size(BP,1)
    temp = LabelMatrix==ii;
    temptemp = sum(temp,2);
    Block(ii, 1) = find(temptemp>0, 1,'first');
    Block(ii, 2) = find(temptemp>0, 1,'last');
    temptemp = sum(temp,1);
    Block(ii, 3) = find(temptemp>0, 1,'first');
    Block(ii, 4) = find(temptemp>0, 1,'last');
end

%imagesc(LabelMatrix)