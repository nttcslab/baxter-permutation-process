function plot_trueBlockLocation(trueBlockLocation)

for ii=1:size(trueBlockLocation,1)
    minR = trueBlockLocation(ii,1);
    maxR = trueBlockLocation(ii,2);
    minC = trueBlockLocation(ii,3);
    maxC = trueBlockLocation(ii,4);
    rectangle('position', [minC,minR,maxC-minC,maxR-minR],'LineWidth',3,'EdgeColor','r');
end
set(gca,'YDir','reverse');