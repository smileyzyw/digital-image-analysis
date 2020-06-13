clear;
close all;
clc;
% 读图像和sift特征
srcFolderPath = 'E:\MATLAB2017b\projects\DIA\第四次作业-特征匹配\Image';
siftFolderPath = 'E:\MATLAB2017b\projects\DIA\第四次作业-特征匹配\sift';
srcAllFiles = dir(srcFolderPath);
siftAllFiles = dir(siftFolderPath);
imgNum = 1000;
img = cell(1,imgNum);
for i = 3 : length(srcAllFiles)
    fileName = srcAllFiles(i).name;
    imgPath = [srcFolderPath, '/', fileName];
    img{i-2} = imread(imgPath);
end

siftDim = 128;
siftCell = cell(1,imgNum); 
paraCell = cell(1,imgNum);
for i = 3 : length(siftAllFiles)
    fileName = siftAllFiles(i).name;
    siftPath = [siftFolderPath, '/', fileName];
    fid = fopen(siftPath, 'rb');
    featNum = fread(fid, 1, 'int32');                   % 本图SIFT特征的数目
    SiftFeat = zeros(siftDim, featNum); 
    paraFeat = zeros(4, featNum);
    for ii = 1 : featNum                                
        SiftFeat(:, ii) = fread(fid, siftDim, 'uchar'); % 128维描述子
        paraFeat(:, ii) = fread(fid, 4, 'float32');     % [x, y, scale, orientation]
    end
    fclose(fid);
    siftCell{i-2} = SiftFeat ./ repmat(sqrt(sum(SiftFeat.^2)), size(SiftFeat, 1), 1); 
    paraCell{i-2} = paraFeat;
end
% 特征拼接
imgNum = 40;
featNum = zeros(1,imgNum);
allFeat10 = cell(1,imgNum);
for i = 1 : imgNum
    allFeat10{1,i} = siftCell{i};
end
allFeat = cell2mat(allFeat10);
% 建立码本
cluNum = [1024 2048 4096 8192];
k1 = cluNum(1);
b = 2;
h = log2(k1);
opts = statset('Display','final','MaxIter',200);
subset = cell(1,2^(h+1)-1);
idxCell = cell(1,2^(h)-1);
cenCell = cell(1,2^(h)-1);
subset{1} = allFeat;
% hierarchical k-means
for i = 1:(2^h-1)
    [idx, centroid] = kmeans(subset{i}',b,'Options',opts);
    idxCell{i} = idx;
    cenCell{i} = centroid;
    subset{i*2} = subset{i}(:,find(idx == 1));
    subset{i*2+1} = subset{i}(:,find(idx == 2));
end
%% 保存/加载码本
% temp = cell2mat(cenCell');
% codebook = temp(end-k1+1:end,:);
% save bow_codebook_1024 codebook
load bow_codebook_1024 codebook
%% 量化
normMat = zeros(1,k1);
[row,col] = size(allFeat);
for i = 1:col
    temp = allFeat(:,i)';
    resdualMat = codebook - repmat(temp,k1,1);
    for j = 1:k1
        normMat(j) = norm(resdualMat(j,:),2);
    end
    qindx(i) = find(normMat == min(normMat)); % 对应量化码字index
%     res(i) = min(normMat);
end
%% 直方图
start = 1;
codeWord = cell(1,length(allFeat10));
for i = 1:length(allFeat10)
    [row,col] = size(allFeat10{i});    
    codeWord{i} = qindx(start:start+col-1);
    start = start + col;
end
xbins1 = 1:k1;
normCounts = cell(1,length(allFeat10));
for i = 1:length(allFeat10)
    [counts,centers] = hist(codeWord{i},xbins1);
    counts = counts./norm(counts,2); % L2 norm 
%     counts = counts./norm(counts,1); % L1 norm
%     bar(centers,counts)
%     xlabel('Feature Index')
%     ylabel('1-Normalized Value')
    normCounts{i} = counts;
end
distMat = zeros(length(allFeat10),length(allFeat10));
for i = 1:length(allFeat10)
    for j = 1:length(allFeat10)
        distMat(i,j) = norm(normCounts{i}-normCounts{j},2);        
    end
end
% 计算匹配精度
indx = zeros(length(allFeat10),4);
for i = 1:length(allFeat10)
    S = sort(distMat(i,:));
    for j = 1:4
        indx(i,j) = find(distMat(i,:)==S(j));
    end
end
result = ceil(indx ./4);
error = 0;
for i = 1:length(allFeat10)
    for j = 1:4
        if result(i,j) ~= result(i,1) % 每一列第一个都是自己
            error = error + 1;
        end 
    end
end
accuracy = 1 - error/(4*length(allFeat10));