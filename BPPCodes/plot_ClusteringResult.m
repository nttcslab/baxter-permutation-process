function plot_ClusteringResult(X, trueBlockLocation, RowLocation, ColumnLocation)

[sortedRowLocation,I_ver] = sort(RowLocation);
[sortedColumnLocation,I_hor] = sort(ColumnLocation);
imagesc(-X(I_ver,I_hor));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap(gray);
%colormap(winter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on;
%temp = A(sort(Loc_ver),sort(Loc_hor));
for ii=1:size(trueBlockLocation)
    ExtractedRow = find( (trueBlockLocation(ii,1)<sortedRowLocation)...
        .* (sortedRowLocation <= trueBlockLocation(ii,2)) );
    ExtractedColumn = find( (trueBlockLocation(ii,3)<sortedColumnLocation)...
        .* (sortedColumnLocation <= trueBlockLocation(ii,4)) );
    minR = min(ExtractedRow)-0.5;
    maxR = max(ExtractedRow)+0.5;
    minC = min(ExtractedColumn)-0.5;
    maxC = max(ExtractedColumn)+0.5;
    if numel(minR)~=0 && numel(maxR)~=0 && numel(minC)~=0 && numel(maxC)~=0
        rectangle('position', [minC,minR,maxC-minC,maxR-minR],'LineWidth',3,'EdgeColor','r');
    end
end
hold off;