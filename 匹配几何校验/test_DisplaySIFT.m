%% Display SIFT features of two images
%
clear; 
close all;
clc;

%% Load two images and their SIFT features
src_1 = './test images/Mona-Lisa-73.jpg';
src_2 = './test images/Mona-Lisa-489.jpg';
ext1 = '.sift';  % extension name of SIFT file
ext2 = '.sift';  % extension name of SIFT file
siftDim = 128;

%% load image 
im_1 = imread(src_1);
im_2 = imread(src_2);

%% load SIFT feature
% SIFT
% feature�ļ���ʽ��binary��ʽ����ͷ�ĸ��ֽڣ�int��Ϊ������Ŀ���������ΪSIFT�����ṹ�壬ÿ��SIFT��������128D�������ӣ�128���ֽڣ�
% �� [x, y, scale, orientation]��16�ֽڵ�λ�á��߶Ⱥ���������Ϣ��float��

featPath_1 = [src_1, ext1];
featPath_2 = [src_2, ext2];

fid_1 = fopen(featPath_1, 'rb');
featNum_1 = fread(fid_1, 1, 'int32'); % �ļ���SIFT��������Ŀ
SiftFeat_1 = zeros(siftDim, featNum_1); 
paraFeat_1 = zeros(4, featNum_1);
for i = 1 : featNum_1 % �����ȡSIFT����
    SiftFeat_1(:, i) = fread(fid_1, siftDim, 'uchar'); %�ȶ���128ά������
    paraFeat_1(:, i) = fread(fid_1, 4, 'float32');     %�ٶ���[x, y, scale, orientation]��Ϣ
end
fclose(fid_1);

fid_2 = fopen(featPath_2, 'rb');
featNum_2 = fread(fid_2, 1, 'int32'); % �ļ���SIFT��������Ŀ
SiftFeat_2 = zeros(siftDim, featNum_2);
paraFeat_2 = zeros(4, featNum_2);
for i = 1 : featNum_2 % �����ȡSIFT����
    SiftFeat_2(:, i) = fread(fid_2, siftDim, 'uchar'); %�ȶ���128ά������
    paraFeat_2(:, i) = fread(fid_2, 4, 'float32');     %�ٶ���[x, y, scale, orientation]��Ϣ
end
fclose(fid_1);

%% normalization
SiftFeat_1 = SiftFeat_1 ./ repmat(sqrt(sum(SiftFeat_1.^2)), size(SiftFeat_1, 1), 1);
SiftFeat_2 = SiftFeat_2 ./ repmat(sqrt(sum(SiftFeat_2.^2)), size(SiftFeat_2, 1), 1);

%% Display SIFT feature on RGB image
[row, col, cn] = size(im_1);
[r2, c2, n2] = size(im_2);
imgBig = 255 * ones(max(row, r2), col + c2, 3);
imgBig(1 : row, 1 : col, :) = im_1;
imgBig(1 : r2, col + 1 : end, :) = im_2; %% ������ͼ��ƴ����һ����ͼ����������

np = 40;
thr = linspace(0, 2*pi, np) ;
Xp = cos(thr); % ȷ��һ����λԲ
Yp = sin(thr);

paraFeat_2(1, :) = paraFeat_2(1, :) + col; % �ڶ���ͼ���е�SIFT feature��������Ҫ�޸�

figure(1); imshow(uint8(imgBig)); axis on;
hold on;
for i = 1 : featNum_1
    xys =  paraFeat_1(:, i);
    if xys(3) < 25 && xys(3) > 3   % ��������̫�ֻ࣬��ʾ�߶���(10, 20)�����ڵ�SIFT    
        figure(1);
        hold on;
        radius = xys(3) * 6;
        plot(xys(1) + Xp * radius, xys(2) + Yp * radius, 'g');
    end
end

for i = 1 : featNum_2
    xys2 = paraFeat_2(:, i);
    if xys2(3) < 25 && xys2(3) > 3  % ��������̫�ֻ࣬��ʾ�߶���(10, 20)�����ڵ�SIFT       
        figure(1);
        hold on;
        radius = xys2(3) * 6;
        plot(xys2(1) + Xp * radius, xys2(2) + Yp * radius, 'g');
    end
end
hold off;