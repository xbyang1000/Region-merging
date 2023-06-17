function [BW, Time] = OurRegionMerge(Im,label,BW)
BW1 = BW;
[row,col,~] = size(Im);
SPixelNum = length(unique(label(:)));
RGBvec = PixelSample(Im);
[FiveTuple,NeighborMatrix]= ObtainSample(SPixelNum,label,RGBvec,BW,row,col);
BaseRegionCenter = FiveTuple(:,3:5);

%[B1,B] = ToVisualization(FiveTuple,Im,BW); %B1 SuperPixel Segmentation; B: SuperPixel Label

tolerance =150;
OriLabel = label;  MultiLabel = [];

step = 1;
tic;
for tol = 1:step:tolerance                      %##
    
    if size(NeighborMatrix,1)<2;break;end    %##
    
    MerRegion = RegionMergeNeiborDis(NeighborMatrix,BaseRegionCenter,tol);
    
    %if (abs(tol-20)<0.001),pause;end  %61 for fire9.jpg  103 for 2-class plane and sky
    if (isempty(MerRegion)),continue;end
    
    [NewLabel,UniqueNewLabel,MultiLabel] = DealMergeInit(MerRegion,OriLabel,MultiLabel);
    OriLabel = NewLabel;
    PixelNum = length(UniqueNewLabel);
    BW = boundarymask(NewLabel);
    for i = 1:PixelNum
        BW( BW == UniqueNewLabel(i) ) = i;
        NewLabel( NewLabel== UniqueNewLabel(i) ) = i;
    end
%     figure(4);hold on;    
%     imshow(imoverlay(Im,BW,'cyan')); %final pixel labels stored in BW
    %pause;

    FiveTuple = [];NeighborMatrix = [];
    [FiveTuple,NeighborMatrix]= ObtainSample(PixelNum,NewLabel,RGBvec,BW,row,col);    
    BaseRegionCenter = FiveTuple(:,3:5);
end
Time = toc;
end

function [FiveTuple,NeighborMatrix]= ObtainSample(SPixelNum,label,RGBvec,BW,row,col)
FiveTuple = zeros(SPixelNum,5);
UniqueLabel = unique(label(:));
NeighborMatrix = zeros(SPixelNum,SPixelNum);
for i = 1:SPixelNum
    Label1 = (label(:) == UniqueLabel(i));
    Label1 = reshape(Label1,[row col]);
    Boundary = Label1 & BW;
    [x,y] = find(Boundary==1);
    NeighborRegion = GetNeighbor(label,[x y],UniqueLabel(i));
    NeighborMatrix(i,NeighborRegion) = 1;
    NeighborMatrix(NeighborRegion,i) = 1;
    
    [x,y] = find(label==UniqueLabel(i));
    Location = mean([x y]);
    ind = find(label(:)==UniqueLabel(i));
    TmpRegion = RGBvec(ind,:);
    FiveTuple(i,:) = [Location mean(TmpRegion)];
end
end
function [label,UniqueLabel,MultiLabel] = DealMergeInit(MerRegion,OriLabel,MultiLabel)
    label = OriLabel;
    UniqueLabel = (unique(OriLabel(:)))';
    n = size(MerRegion,1);
    k = size(MultiLabel,1);
    for i = 1:1:n
        R = MerRegion{i,1};
        for ii = 1:size(R,1)        
            RR = R(ii,:);
            RR = UniqueLabel(RR);
            [MultiLabel,FlagLabel] = BelongToRegion(RR,MultiLabel);
            for j = 2:length(FlagLabel),  label( label== FlagLabel(j)) = FlagLabel(1);   end        
        end
    end  
    UniqueLabel = (unique(label(:)));
end

function [MultiLabel,FlagLabel] = BelongToRegion(R,MultiLabel)
Flag = false;
FlagLabel = R;
n = size(MultiLabel,1);
for i=1:n
    Region = MultiLabel{i,1};
    Set = intersect(R,Region);
    if (~isempty(Set))
        MultiLabel{i,1} = sort(union(MultiLabel{i,1},R));
        FlagLabel = MultiLabel{i,1}; 
        Flag = true;
        break;
    end
end
if (~Flag)
    n = n+1;
    MultiLabel{n,1} = R;
end
end

function [B1,B] =ToVisualization(FiveTuple,Im,BW)
figure(2); hold on;
B1 = imoverlay(Im,BW,'cyan');


Location = FiveTuple(:,1:2);
str = cell(size(FiveTuple,1),1);
for i = 1:size(FiveTuple,1)
    str{i} = num2str(i);
end

hold on;
L = floor([Location(:,2) Location(:,1)]);
B= insertText(B1,L,str,'BoxOpacity',0.0,'TextColor','r','FontSize',8,'Anchorpoint','center');
%B = imoverlay(RGB,BW,'cyan');
imshow(B);%,'InitialMagnification',67);
%B = insertText(Im,L,str,'BoxOpacity',0.0,'TextColor','r','FontSize',8,'Anchorpoint','center');
%imwrite(B,'Tmp.bmp'); %write to image file
end

function NeighborRegion = GetNeighbor(label,location,i)
[row,col] = size(label);
x = location(:,1); y = location(:,2);
Loc = [x-1, y;x-1 y-1; x-1 y+1 ];
Loc = [Loc; x, y; x,y-1; x,y+1 ];
Loc = [Loc; x+1,y-1; x+1,y; x+1, y+1];
Loc(Loc(:,1)<1,1) = 1;
Loc(Loc(:,2)<1,2) = 1;
Loc(Loc(:,1)>row,1) = row;
Loc(Loc(:,2)>col,2) = col;

ind = sub2ind(size(label),Loc(:,1),Loc(:,2));
region = label(ind);
% region = [];
% for i=1:size(Loc,1)
%     region = union(region,label(Loc(i,1),Loc(i,2)));
% end

NeighborRegion = unique(region);
NeighborRegion = setdiff(NeighborRegion,i);
end

function RGBvec = PixelSample(Im)
R = Im(:,:,1);  G =  Im(:,:,2);   B = Im(:,:,3);
R = R(:);  G = G(:);   B = B(:);
RGBvec = double([R G B]);
end