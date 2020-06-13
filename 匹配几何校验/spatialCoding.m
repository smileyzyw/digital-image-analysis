%% Set path and parameters
clear;
close all;
clc;

% src_1 = './test images/37967br1.jpg';  
% src_2 = './test images/791.jpg';

% src_1 = './test images/4.jpg';  
% src_2 = './test images/Apollo-266.jpg';

src_1 = './test images/771.jpg';  
src_2 = './test images/305.jpg';

% src_1 = './test images/Apollo-49.jpg';
% src_2 = './test images/Apollo-266.jpg';

% src_1 = './test images/Disney-00524.jpg';  
% src_2 = './test images/Disney-00550.jpg';

ext = '.sift'; % extension name of SIFT file
siftDim = 128;
maxAxis = 400;

%  Load image
im_1 = imread(src_1);
if max(size(im_1)) > maxAxis
    im_1 = imresize(im_1, maxAxis / max(size(im_1)));
end

im_2 = imread(src_2);
if max(size(im_2)) > maxAxis
    im_2 = imresize(im_2, maxAxis / max(size(im_2)));
end

%  Load SIFT feature from file
featPath_1 = [src_1, ext];
featPath_2 = [src_2, ext];

fid_1 = fopen(featPath_1, 'rb');
featNum_1 = fread(fid_1, 1, 'int32');
SiftFeat_1 = zeros(siftDim, featNum_1);
paraFeat_1 = zeros(4, featNum_1);
for i = 1 : featNum_1
    SiftFeat_1(:, i) = fread(fid_1, siftDim, 'uchar');
    paraFeat_1(:, i) = fread(fid_1, 4, 'float32');
end
fclose(fid_1);

fid_2 = fopen(featPath_2, 'rb');
featNum_2 = fread(fid_2, 1, 'int32');
SiftFeat_2 = zeros(siftDim, featNum_2);
paraFeat_2 = zeros(4, featNum_2);
for i = 1 : featNum_2
    SiftFeat_2(:, i) = fread(fid_2, siftDim, 'uchar');
    paraFeat_2(:, i) = fread(fid_2, 4, 'float32');
end
fclose(fid_1);

% Normalization
SiftFeat_1 = SiftFeat_1 ./ repmat(sqrt(sum(SiftFeat_1.^2)), size(SiftFeat_1, 1), 1);
SiftFeat_2 = SiftFeat_2 ./ repmat(sqrt(sum(SiftFeat_2.^2)), size(SiftFeat_2, 1), 1);

% Check match based on distances between SIFT descriptors across images
normVal = mean(sqrt(sum(SiftFeat_1.^2)));
matchInd = zeros(featNum_1, 1);
matchDis = zeros(featNum_1, 1);
validDis = [];
gridDisVec = [];
ic = 0;
for i = 1 : featNum_1
    tmpFeat = repmat(SiftFeat_1(:, i), 1, featNum_2);
    d = sqrt(sum((tmpFeat - SiftFeat_2).^2)) / normVal; % L2 distance
    matchDis(i) = min(d);
    [v, ind] = sort(d);
    if v(1) < 0.4               % 最小距离小于0.4，则认为构成一对匹配
        matchInd(i) = ind(1);   % match上的第二张图的index
        ic = ic + 1;
        validDis(ic, 1 : 3) = [v(1), v(2), v(1) / v(2)];
        tmp = (SiftFeat_1(:, i) - SiftFeat_2(:, ind(1))).^2;
        tmp2 = reshape(tmp(:), 8, 16);
        gridDisVec(ic, 1 : 16) = sqrt(sum(tmp2)); % 两个SIFT特征的16个格子对应的距离
    end
end
% figure; stem(matchDis); ylim([0, 1.2]); % 所有的点
% figure; stem(matchDis(matchInd > 0)); ylim([0, 1.2]); % 小于0.4的点

%% spatial coding 
matchLen = ic;
matchInd1 = find(matchInd > 0); % matched feature1的id
matchInd2 = matchInd(matchInd1);
Xmap1 = zeros(matchLen,matchLen);
Ymap1 = zeros(matchLen,matchLen);
Xmap2 = zeros(matchLen,matchLen);
Ymap2 = zeros(matchLen,matchLen);
for i = 1:ic
    for j = 1:ic
        if paraFeat_1(1,matchInd1(i)) >= paraFeat_1(1,matchInd1(j))
            Xmap1(i,j) = 1;
        else
            Xmap1(i,j) = 0;
        end
    end
end

for i = 1:ic
    for j = 1:ic
        if paraFeat_1(2,matchInd1(i)) >= paraFeat_1(2,matchInd1(j))
            Ymap1(i,j) = 1;
        else
            Ymap1(i,j) = 0;
        end
    end
end
for i = 1:ic
    for j = 1:ic
        if paraFeat_2(1,matchInd2(i)) >= paraFeat_2(1,matchInd2(j))
            Xmap2(i,j) = 1;
        else
            Xmap2(i,j) = 0;
        end
    end
end

for i = 1:ic
    for j = 1:ic
        if paraFeat_2(2,matchInd2(i)) >= paraFeat_2(2,matchInd2(j))
            Ymap2(i,j) = 1;
        else
            Ymap2(i,j) = 0;
        end
    end
end

%%
indRecord = ones(ic,2);
indRecord(:,1) = matchInd1;
indRecord(:,2) = matchInd2;
indVec = ones(ic,1);
Xmap1_temp = Xmap1;
Xmap2_temp = Xmap2;
Ymap1_temp = Ymap1;
Ymap2_temp = Ymap2;

Vx = xor(Xmap1_temp,Xmap2_temp);
Vy = xor(Ymap1_temp,Ymap2_temp);

for i = 1:floor(ic/2) % x,y轴都删，一次删2个点
    Sx = sum(Vx,2);
    if sum(Sx) ~= 0
        xid = find(Sx == max(Sx));
        indRecord(xid(1),:) = []; % 删除对应的匹配对
        Vx(xid(1),:) = [];
        Vx(:,xid(1)) = [];
        Vy(xid(1),:) = [];
        Vy(:,xid(1)) = [];
        Sy = sum(Vy);
        yid = find(Sy == max(Sy));
        indRecord(yid(1),:) = []; % 删除对应的匹配对
        Vx(yid(1),:) = [];
        Vx(:,yid(1)) = [];
        Vy(yid(1),:) = [];
        Vy(:,yid(1)) = [];
    end
end
[rr,cc] = size(indRecord);


% Sy = sum(Vy);

%% Show the local matching results on RGB image
[row, col, cn] = size(im_1);
[r2, c2, n2] = size(im_2);
imgBig = 255 * ones(max(row, r2), col + c2, 3); % 大ccanvas
imgBig(1 : row, 1 : col, :) = im_1;
imgBig(1 : r2, col + 1 : end, :) = im_2;
np = 40;
thr = linspace(0,2*pi,np) ;
Xp = cos(thr);
Yp = sin(thr);
paraFeat_2_new = paraFeat_2;
paraFeat_2_new(1, :) = paraFeat_2_new(1, :) + col; % 第二张图的列坐标
figure(3); imshow(uint8(imgBig)); axis on;
hold on;
matchCount = rr;
for i = 1 : ic
    xys = paraFeat_1(:, matchInd1(i));
    xys2 = paraFeat_2_new(:, matchInd2(i));
    figure(3);
    hold on;
%     I = find([matchInd1(i) matchInd2(i)] == indRecord);
    if find(indRecord(:,1) == matchInd1(i)) ~= 0    % 匹配特征（蓝色）
        plot(xys(1) + Xp * xys(3) * 6, xys(2) + Yp * xys(3) * 6, 'b'); % 画圆
        plot(xys2(1) + Xp * xys2(3) * 6, xys2(2) + Yp * xys2(3) * 6, 'b');
        hold on; plot([xys(1), xys2(1)], [xys(2), xys2(2)], '-b', 'LineWidth', 0.8); % 连线        
    else                                            % 非匹配特征（红色）
%         plot(xys(1) + Xp * xys(3) * 6, xys(2) + Yp * xys(3) * 6, 'r'); % 画圆
%         plot(xys2(1) + Xp * xys2(3) * 6, xys2(2) + Yp * xys2(3) * 6, 'r');
%         hold on; plot([xys(1), xys2(1)], [xys(2), xys2(2)], '-r', 'LineWidth', 0.8); % 连线    
    end
end
figure(3);
% title(sprintf('Total local matches : %d (%d-%d)', length(find(matchInd)), featNum_1 ,featNum_2));
hold off;

