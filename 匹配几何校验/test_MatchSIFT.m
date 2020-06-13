
%% Set path and parameters
clear;
close all;
clc;

src_1 = './test images/37967br1.jpg';  
src_2 = './test images/791.jpg';

% src_1 = './test images/4.jpg';  
% src_2 = './test images/Apollo-266.jpg';

% src_1 = './test images/771.jpg';  
% src_2 = './test images/305.jpg';

% src_1 = './test images/Apollo-49.jpg';
% src_2 = './test images/Apollo-266.jpg';


ext = '.sift'; % extension name of SIFT file
siftDim = 128;
maxAxis = 400;


%%  Load image
im_1 = imread(src_1);
if max(size(im_1)) > maxAxis
    im_1 = imresize(im_1, maxAxis / max(size(im_1)));
end

im_2 = imread(src_2);
if max(size(im_2)) > maxAxis
    im_2 = imresize(im_2, maxAxis / max(size(im_2)));
end


%%  Load SIFT feature from file
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


%%Normalization
SiftFeat_1 = SiftFeat_1 ./ repmat(sqrt(sum(SiftFeat_1.^2)), size(SiftFeat_1, 1), 1);
SiftFeat_2 = SiftFeat_2 ./ repmat(sqrt(sum(SiftFeat_2.^2)), size(SiftFeat_2, 1), 1);


%% Check match based on distances between SIFT descriptors across images
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
    if v(1) < 0.4            % 最小距离小于0.4，则认为构成一对匹配
        matchInd(i) = ind(1);
        ic = ic + 1;
        validDis(ic, 1 : 3) = [v(1), v(2), v(1) / v(2)];
        tmp = (SiftFeat_1(:, i) - SiftFeat_2(:, ind(1))).^2;
        tmp2 = reshape(tmp(:), 8, 16);
        gridDisVec(ic, 1 : 16) = sqrt(sum(tmp2));
    end
end
figure; stem(matchDis); ylim([0, 1.2]);
figure; stem(matchDis(matchInd > 0)); ylim([0, 1.2]);


%% Show the local matching results on RGB image
[row, col, cn] = size(im_1);
[r2, c2, n2] = size(im_2);
imgBig = 255 * ones(max(row, r2), col + c2, 3);
imgBig(1 : row, 1 : col, :) = im_1;
imgBig(1 : r2, col + 1 : end, :) = im_2;
np = 40;
thr = linspace(0,2*pi,np) ;
Xp = cos(thr);
Yp = sin(thr);
paraFeat_2(1, :) = paraFeat_2(1, :) + col;
figure(3); imshow(uint8(imgBig)); axis on;
hold on;
matchCount = 0;
for i = 1 : featNum_1
    if matchInd(i) > 0
        matchCount = matchCount + 1;
        xys = paraFeat_1(:, i);
        xys2 = paraFeat_2(:, matchInd(i));
        figure(3);
        hold on;
        plot(xys(1) + Xp * xys(3) * 6, xys(2) + Yp * xys(3) * 6, 'r');
        plot(xys2(1) + Xp * xys2(3) * 6, xys2(2) + Yp * xys2(3) * 6, 'r');
        hold on; plot([xys(1), xys2(1)], [xys(2), xys2(2)], '-b', 'LineWidth', 0.8);
        
%         figure(4);
%         clf;
%         subplot(311); stem(SiftFeat_1(:, i)); xlim([0, 128]); ylim([0, 0.5]);
%         title(sprintf('Feature pair %d', matchCount));
%         subplot(312); stem(SiftFeat_2(:, matchInd(i))); xlim([0, 128]); ylim([0, 0.5]);
%         tmp = SiftFeat_1(:, i) - SiftFeat_2(:, matchInd(i));
%         subplot(313); stem(tmp); xlim([0, 128]);
%         title('Difference between the above two features per dimension');
%         disVal = sum(sqrt(tmp.^2));
%         ylim([0, 0.2]);
    end
end
figure(3);
title(sprintf('Total local matches : %d (%d-%d)', length(find(matchInd)), featNum_1 ,featNum_2));
hold off;

