function [Left2RightMaxima, Right2LeftMaxima, Left2RightLocation, Right2LeftLocation]...
    = ExtractMaxima(BP)

n = size(BP,1);
Left2RightMaxima = BP(1);
Left2RightLocation = 1;
Right2LeftMaxima = BP(end);
Right2LeftLocation = size(BP,1);
for ii=2:n
    if BP(ii) > Left2RightMaxima(end)
        Left2RightMaxima = [Left2RightMaxima; BP(ii)];
        Left2RightLocation = [Left2RightLocation; ii];
    end
    if BP(n+1-ii) > Right2LeftMaxima(end)
        Right2LeftMaxima = [Right2LeftMaxima; BP(n+1-ii)];
        Right2LeftLocation = [Right2LeftLocation; n+1-ii];
    end
end
