function result = plot_rectangles(data, A, Loc_ver, Loc_hor)
[bbb,I_ver] = sort(Loc_ver);
[bbb,I_hor] = sort(Loc_hor);
imagesc(-data(I_ver,I_hor));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap(gray);
%colormap(winter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on;
%temp = A(sort(Loc_ver),sort(Loc_hor));
temp = A(I_ver,I_hor);
for ii=1:max(A(:))
    temptemp=temp;
    temptemp(temptemp~=ii)=0;
    %minR = min(find(sum(temptemp,2)>0))-0.5;
    %maxR = max(find(sum(temptemp,2)>0))+0.5;
    %minC = min(find(sum(temptemp,1)>0))-0.5;
    %maxC = max(find(sum(temptemp,1)>0))+0.5;
    minR = find( sum(temptemp,2)>0 ,1)-0.5;
    maxR = find(sum(temptemp,2)>0, 1, 'last')+0.5;
    minC = find(sum(temptemp,1)>0, 1)-0.5;
    maxC = find(sum(temptemp,1)>0, 1, 'last')+0.5;
    if numel(minR)~=0 && numel(maxR)~=0 && numel(minC)~=0 && numel(maxC)~=0
    rectangle('position', [minC,minR,maxC-minC,maxR-minR],'LineWidth',3,'EdgeColor','r');
    end
end
hold off;
%daspect([1 1 1]);
result = 1;