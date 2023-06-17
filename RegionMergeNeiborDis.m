function MerRegion = RegionMergeNeiborDis(NeighborMatrix,BaseRegionCenter,tolerance)


%filename = 'polygon_small.tif';
%filename = 'fire1.jpg';
%filename = 'rectangle.tif';
% [row,col,~] = size(Im);
% 
% %% Get boundary points by the number of superpixels
% NeighborMatrix = zeros(SPixelNum,SPixelNum);
% BW = boundarymask(label);
% RGBvec = PixelSample(Im);
% 
% % figure(1);hold on;
% % imshow(imoverlay(Im,BW,'cyan'),'InitialMagnification',67);
% FiveTuple = zeros(SPixelNum,5);
% for i = 1:SPixelNum
%     Label1 = (label(:) == i);
%     Label1 = reshape(Label1,[row col]);
%     Boundary = Label1 & BW;
%     [x,y] = find(Boundary==1);
%     NeighborRegion = GetNeighbor(label,[x y],i);
%     NeighborMatrix(i,NeighborRegion) = 1;
%     NeighborMatrix(NeighborRegion,i) = 1;
%     
%     [x,y] = find(label==i);
%     Location = mean([x y]);
%     ind = find(label(:)==i);
%     TmpRegion = RGBvec(ind,:);
%     FiveTuple(i,:) = [Location median(TmpRegion)];
% end
% 
% BaseRegionCenter = FiveTuple(:,3:5);
%[B1,B] = ToVisualization(FiveTuple,Im,BW); 

MerRegion = GetMergeSet(NeighborMatrix,BaseRegionCenter,tolerance);
% ToShowFinalResult(MerRegion,label,BW,B1);
% ToShowMergeRegion(MerRegion,label,B,SPixelNum);


end
function ToShowMergeRegion(MerRegion,label,Im,SPixelNum)
FullSet = 1:1:SPixelNum;
SelectSet = [];
ImI = Im;
n = size(MerRegion,1);
for i = n:-1:1
    R = MerRegion{i,1};
    for ii = 1:size(R,1)
        LabelMatrix = zeros(size(label,1),size(label,2));
        RR = R(ii,:);
        SelectSet = union(SelectSet,RR);
        for j = 1:length(RR)
            TmpLabel = RR(j);
            LabelMatrix( label(:)==TmpLabel ) =1;
        end
    LabelMatrix = boundarymask(LabelMatrix);    
    ImI = imoverlay(ImI,LabelMatrix,'c'); 
    end

end
figure(4);hold on;
% imshow(ImI);
RestSet = setdiff(FullSet,SelectSet);

for i = 1:length(RestSet) 
    LabelMatrix = zeros(size(label,1),size(label,2));
    LabelMatrix( label(:)==RestSet(i) ) =1;
    LabelMatrix = boundarymask(LabelMatrix);
    ImI = imoverlay(ImI,LabelMatrix,'c'); 
end

hold on;
imshow(ImI);
end

function ToShowFinalResult(MerRegion,label,BW,Im)
%Im = imread('Tmp.bmp');
ImI = Im;
n = size(MerRegion,1);
figure(3);hold on;
ColorMap = 'gymrbgymrbgymrb';

for i = n:-1:1
%     LabelMatrix = zeros(size(label,1),size(label,2));
    R = MerRegion{i,1};
    for ii = 1:size(R,1)
        LabelMatrix = zeros(size(label,1),size(label,2));
        RR = R(ii,:);
        for j = 1:length(RR)
            TmpLabel = RR(j);
            LabelMatrix( label(:)==TmpLabel ) =1;
        end
    LabelMatrix = boundarymask(LabelMatrix);
    %LabelMatrix = LabelMatrix & BW;
    ImI = imoverlay(ImI,LabelMatrix,ColorMap(i));
    hold on;
    imshow(ImI);%,'InitialMagnification',67); 
   %pause;
    end
%     LabelMatrix = boundarymask(LabelMatrix);
%     LabelMatrix = LabelMatrix & BW;
%     ImI = imoverlay(ImI,LabelMatrix,ColorMap(i));
%     hold on;
%     imshow(ImI,'InitialMagnification',67);     
%     pause;
end
%     LabelMatrix = boundarymask(LabelMatrix);
%     LabelMatrix = LabelMatrix & BW;
%     ImI = imoverlay(ImI,LabelMatrix,ColorMap(1));
%     hold on;
%     imshow(ImI,'InitialMagnification',67);
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
B = insertText(Im,L,str,'BoxOpacity',0.0,'TextColor','r','FontSize',8,'Anchorpoint','center');
%imwrite(B,'Tmp.bmp'); %write to JPG file
end


function ConvexRegion = GetMergeSet(NeighborMatrix,BaseRegionCenter,tolerance)

n = size(BaseRegionCenter,1);
RecursionFlag = 1;
LastConvexRegion = [];
ConvexRegion = (1:1:n)';
TmpConvexRegion = ConvexRegion;
k = 0;
%% Obtain the maximum convex set, may be not unique
while(RecursionFlag)

        %TmpConvexRegion = ConvexSetMerge(NeighborMatrix,ConvexRegion,TmpConvexRegion,BaseRegionCenter,tolerance);
        
        TmpConvexRegion = RegionMerge(NeighborMatrix,ConvexRegion,TmpConvexRegion,BaseRegionCenter,tolerance);
        
%     for i = 1:size(TmpConvexRegion)
%       TmpConvexRegion(i,:)=sort(TmpConvexRegion(i,:));    
%     end
    [C,~,~] = unique(TmpConvexRegion,'rows'); %
    TmpConvexRegion = C;    
    if ( isempty(TmpConvexRegion) || TmpConvexRegion(1)==-1)
        RecursionFlag = 0;
    end
    if (RecursionFlag)
        k = k+1;        
        LastConvexRegion{k,1} = TmpConvexRegion;
    end
end

%% Backtracking k>0
%load Tmp1221;
ConvexRegion = [];
if (isempty(LastConvexRegion)), return; end
LastLayerRegion = LastConvexRegion{k,1};
MaxConvexRegion = ObtainMaxConvex(NeighborMatrix,BaseRegionCenter,LastLayerRegion);

ConvexRegion{1,1} = MaxConvexRegion;
LastConvexRegion(end) = [];


LastConvexRegion = RemoveRepeatSubset(ConvexRegion,LastConvexRegion);
k = size(LastConvexRegion,1);

layer = 1;
for i = k:-1:1
    if (isempty(LastConvexRegion)), continue; end
    CurrentLayerRegion = LastConvexRegion{i,1};
    if ( size(CurrentLayerRegion,1)<1 )        
        LastConvexRegion(end) = [];
        continue;
    end
    if (size(CurrentLayerRegion,1)==1) 
        layer = layer + 1;        
        ConvexRegion{layer,1} = CurrentLayerRegion;
    end
    if (size(CurrentLayerRegion,1)>1)
        MaxConvexRegion = ObtainMaxConvex(NeighborMatrix,BaseRegionCenter,CurrentLayerRegion);
        layer = layer + 1;        
        ConvexRegion{layer,1} = MaxConvexRegion;        
    end
    %LastConvexRegion = LastConvexRegion(1:k-1,1);
    LastConvexRegion = RemoveRepeatSubset(ConvexRegion,LastConvexRegion);
end
end

function LastConvexRegion = RemoveRepeatSubset(ConvexRegion,LastConvexRegion)
    for i = 1:size(ConvexRegion,1) %get cell, same dimension       
        R = ConvexRegion{i,1};
        for ii = 1:size(R,1)
            CurrentConvexRegion = R(ii,:);
            if (isempty(LastConvexRegion)), continue;end
             for j=1:size(LastConvexRegion,1) %get cell, same dimension
                 RR = LastConvexRegion{j,1};
                 
                 index = []; %FullSet = 1:1:size(RR,1);
                 for jj= 1:size(RR,1)
                     CandidateRegion = RR(jj,:);
                     SubsetResult = intersect(CurrentConvexRegion,CandidateRegion);
                     if (length(SubsetResult)<1) %disjoint
                         index = union(index,jj);
                     end
                 end
                 RR = RR(index',:);                  
                 LastConvexRegion{j,1} = RR; 
             end
        end
    end   

end


function TmpConvexRegion= ObtainMaxConvex(NeighborMatrix,BaseRegionCenter,LastLayerRegion)
TmpConvexSet = [];TmpConvexRegion = [];
n = size(LastLayerRegion,1);

DisMatrix = zeros(n,n)-1; 
DisjointSet = [];JointSet = [];
DisJnt = 0; Joint = 0;
DisjointIndex = [];

for i = 1:n-1
    for j = i+1:n
        set1 = LastLayerRegion(i,:);
        set2 = LastLayerRegion(j,:);
        TmpSet = intersect(set1,set2);
        RegionVertex1 = set1';  m1 = mean(BaseRegionCenter(RegionVertex1,:));
        RegionVertex2 = set2';  m2 = mean(BaseRegionCenter(RegionVertex2,:));
        DisMatrix(i,j) = sqrt((m1-m2)*(m1-m2)'); %DisMatrix(j,i) = DisMatrix(i,j);
        
        if (length(TmpSet)> 0),continue;end        
        
        if (length(TmpSet) <= 0) % neighor or disjoint, i.e., no same triangles
            DisJnt = DisJnt+1;
            DisjointSet(DisJnt,:) = [i j];
        else
            Joint = Joint+1;
            JointSet(Joint,:) = [i j];
        end
    end
end

if (isempty(DisjointSet)) % no disjoint maximum convex set, return the first one as the maximum convex set
     TmpConvexRegion = LastLayerRegion(1,:);
    return;
end

% if existing multiple maximum convex set,firstly put the longest two into temporary sets named TmpConvexSet and TmpConvexRegion
% tmp 
[x,y]  = find(DisMatrix == max(max(DisMatrix)));  %[3,5]
flag = ismember([x y],DisjointSet,'rows');
if (flag)    
    TmpConvexRegion = LastLayerRegion([x y]',:);
    %DisMatrix(x,y) = -2;
else
    % only one maximum convex set, return    
    TmpConvexRegion = LastLayerRegion(x,:);
    return;
end

ind = ismember(DisjointSet,[x y],'rows');
DisjointSet(ind,:) = [];   % remove [3,5] 

% restore distance matrix (Upper triangle) to full matrix
for i=1:n-1
    DisMatrix(i,i) = inf;
    for j=i+1:n
        DisMatrix(j,i) = DisMatrix(i,j);
    end
end
DisMatrix(n,n) = inf;

Current_Set = [x y];
Full_Set = DisjointSet(:);
Full_Set = (unique(Full_Set))';
Rest_Set = setdiff(Full_Set,Current_Set);



while(~isempty(Rest_Set))
   % obtain the most similar regions to selected ones from un-selected regions.
   row = Current_Set;
   col = Rest_Set;
   TmpDis = DisMatrix(row,col);
   minValue = min(min(TmpDis));
   if (minValue<0),continue;end
   [x,y] = find(TmpDis==minValue);
   Selected = row(x(1));
   UnSelected = col(y(1));
   Rest_Set = setdiff(Rest_Set,UnSelected);% remove unselected one since it has been investiaged
   SelectedPair = [Selected UnSelected];
   SelectedPair = sort(SelectedPair);% ascending order
   
   
   % Does it appear in disjoint set?
   ind = ismember(DisjointSet,SelectedPair,'rows');
   if (sum(ind)>0) %Yes
       Current_Set = union(Current_Set,UnSelected);      
      TmpConvexRegion = [TmpConvexRegion;LastLayerRegion(UnSelected,:)];
   end
end
end

% TmpConvexRegion = ConvexSetMerge(NeighborMatrix,ConvexRegion,TmpConvexRegion,BaseRegionCenter,tolerance);;
function TmpConvexRegion = ConvexSetMerge(NeighborMatrix,ConvexRegion,LastConvexRegion,BaseRegionCenter,tolerance)
TmpConvexRegion=[];
k=0;

for j = 1:size(LastConvexRegion,1)
    set2 = LastConvexRegion(j,:);
    tmpNeighborMatrix = NeighborMatrix(set2',:);
    [~,y] = find(tmpNeighborMatrix==1);
    neighbors = unique(y');
    neighbors = setdiff(neighbors,set2');
    Len = length(neighbors);
    for i = 1:Len
        set1 = ConvexRegion(neighbors(i),:);  
       if (~isempty(intersect(set1,set2))),continue;end
       if (isNeighbor(set1,set2,NeighborMatrix,BaseRegionCenter,tolerance))% && isSimilar (set1,set2,BaseRegionCenter,tolerance) )
           set1 = union (set1,set2);  
           k = k+1;           
           Convex2Region(k,:) = set1;
       end   
   end
end   
if k>0   
    TmpConvexRegion = Convex2Region;
end
    
end

function Flag = isSimilar (set1,set2,BaseRegionCenter,tolerance)
Flag = false;
%    for i=1:length(set1)
%        for j=1:length(set2)
% set1 = set1';
% set2 = set2';
% %TmpSet = union(set1,set2);
% X = BaseRegionCenter(set2',:);
% n = size(X,1);
% X_mean = mean(X);
% 
% VarValue = trace(X*X')-n*(X_mean*X_mean');
% VarValue = sqrt(VarValue);
Center1 = BaseRegionCenter(set1',:);
if (length(set2)>1)
     Center2 = mean(BaseRegionCenter(set2',:));
else
    Center2 = BaseRegionCenter(set2',:);
end
Center = Center1-Center2;
Center = sqrt(Center*Center');

           if (Center<tolerance)% && VarValue < 2*tolerance)
               Flag = true;
           end

end

function Flag =  isNeighbor(set1,set2,NeighborMatrix,BaseRegionCenter,tolerance)
Flag = false;
X1 = BaseRegionCenter(set2',:);
X  = BaseRegionCenter(set1,:);
Dis = pdist2(X1,X);
ind = find(Dis<tolerance);
set2ind = set2(ind');
if (sum(NeighborMatrix(set1,set2ind'))>0)    
     Flag = true;
end
% for i=1:length(set1)
%        for j=1:length(set2)
%            if (NeighborMatrix(set1(i),set2(j)) > 0 )
%                Flag = true;
%                return;
%            end
%        end
% end
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
