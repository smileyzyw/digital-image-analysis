clear;
close all;
clc;
% ��ͼ���sift����
srcFolderPath = 'E:\MATLAB2017b\projects\DIA\���Ĵ���ҵ-����ƥ��\Image';
siftFolderPath = 'E:\MATLAB2017b\projects\DIA\���Ĵ���ҵ-����ƥ��\sift';
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
    featNum = fread(fid, 1, 'int32');                   % �ļ���SIFT��������Ŀ
    SiftFeat = zeros(siftDim, featNum); 
    paraFeat = zeros(4, featNum);
    for ii = 1 : featNum                                 
        SiftFeat(:, ii) = fread(fid, siftDim, 'uchar'); 
        paraFeat(:, ii) = fread(fid, 4, 'float32');     
    end
    fclose(fid);
    siftCell{i-2} = SiftFeat ./ repmat(sqrt(sum(SiftFeat.^2)), size(SiftFeat, 1), 1); 
    paraCell{i-2} = paraFeat;
end
%% ����ƴ��
imgNum = 40;
featNum = zeros(1,imgNum);
allFeat10 = cell(1,imgNum);
for i = 1 : imgNum
    allFeat10{1,i} = siftCell{i};
end
allFeat = cell2mat(allFeat10);
%% �����뱾
cluNum = [8 16 32 64];
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

%% �����뱾
temp = cell2mat(cenCell');
codebook = temp(end-k1+1:end,:);
save vlad_codebook_8 codebook
%% �����뱾
% load bow_codebook_1024 codebook
%% ����
normMat = zeros(1,k1);
[row,col] = size(allFeat);
for i = 1:col
    temp = allFeat(:,i)';
    resdualMat = codebook - repmat(temp,k1,1);
    for j = 1:k1
        normMat(j) = norm(resdualMat(j,:),2);
    end
    qindx(i) = find(normMat == min(normMat)); % ��Ӧ��������index
%     res(i) = min(normMat);
end
%% VLAD
start = 1;
codeWord = cell(1,length(allFeat10));
for i = 1:length(allFeat10)
    [row,col] = size(allFeat10{i});    
    codeWord{i} = qindx(start:start+col-1); % codeWord��ÿ��ͼ��������������
    start = start + col;
end
% ��ʼ��Ϊ0
resSum = cell(1,imgNum);
normResSum = cell(1,imgNum);
for i = 1:imgNum
    resSum{i} = zeros(128,k1);
    normResSum{i} = zeros(128,k1);
end

for i = 1:imgNum
    for j = 1:length(codeWord{i}) 
        idx_cobk = codeWord{i}(j); % ��ʾ�뱾����ĵڼ�������
        resSum{i}(:,idx_cobk) = allFeat10{i}(:,j)-codebook(idx_cobk,:)' + resSum{i}(:,idx_cobk);
    end
end

%L2 normalization
for i = 1:length(allFeat10)
    for j = 1:k1
       normResSum{i}(:,j) = resSum{i}(:,j)/(norm(resSum{i}(:,j),2)+1e-5);
    end   
end
distMat = zeros(length(allFeat10),length(allFeat10));
for i = 1:length(allFeat10)
    for j = 1:length(allFeat10)
        distMat(i,j)=norm((normResSum{i}-normResSum{j}),'fro');       
    end
end
% ����ƥ�侫��
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
        if result(i,j) ~= result(i,1) % ÿһ�е�һ�������Լ�
            error = error + 1;
        end 
    end
end
accuracy = 1 - error/(4*length(allFeat10));