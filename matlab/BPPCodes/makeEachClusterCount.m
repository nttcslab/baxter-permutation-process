function EachClusterCount = makeEachClusterCount(X, trueBlockLocation, RowLocation, ColumnLocation)

ObservationRange = max(X(:));
EachClusterCount = zeros(ObservationRange, size(trueBlockLocation,1));

for ii = 1:size(trueBlockLocation)
    ExtractedRow = find( (trueBlockLocation(ii,1)<RowLocation)...
        .* (RowLocation <= trueBlockLocation(ii,2)) );
    ExtractedColumn = find( (trueBlockLocation(ii,3)<ColumnLocation)...
        .* (ColumnLocation <= trueBlockLocation(ii,4)) );
    ExtractedData = X(ExtractedRow, ExtractedColumn);
    for jj=1:ObservationRange
        EachClusterCount(jj,ii) = sum(sum( (ExtractedData==jj) ));
    end
end

