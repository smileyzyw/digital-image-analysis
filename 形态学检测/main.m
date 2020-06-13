clear;
close all;
clc;

%% Binarize the input image
im = cell(1,5);
im_gray = cell(1,5);
level = zeros(1,5);
bw = cell(1,5);

im{1} = imread('test images/barcode_1.png');
im{2} = imread('test images/barcode_3.png');
im{3} = imread('test images/barcode_4.png');
im{4} = imread('test images/barcode_6.png');
im{5} = imread('test images/barcode_7.png');
for i = 1:5
    im_gray{i} = rgb2gray(im{i});
    level(i) = graythresh(im_gray{i});
    bw{i} = im2bw(im_gray{i}, level(i));
end
%% Date:2020/3/12; By zyw
% 腐蚀一次去除不相关的文字和图案,背景为目标，腐蚀膨胀反过来
se = strel('rectangle',[50,10]);
% 这里读的是第一张图像，修改数字1 -> 2,3,4,5运行其他几张图像
gray_img = im_gray{5};
bin_img = bw{5};
% imshow(bin_img)
BW2 = imdilate(bin_img,se);
[imgRow,imgCol] = size(BW2);
% 图像对半分别找出前景的区间
[up_row,up_col] = find(BW2(1:round(imgRow/2),:)==0);
up_max_row = max(up_row) + 26;
up_min_row = min(up_row) - 26;
[dw_row,dw_col] = find(BW2(round(imgRow/2):imgRow,:)==0);
dw_max_row = max(dw_row) + 26 + round(imgRow/2); %腐蚀
dw_min_row = min(dw_row) - 26 + round(imgRow/2);
% 取出条形码部分
temp = ones(imgRow,imgCol);
temp(up_min_row:up_max_row,:) = bin_img(up_min_row:up_max_row,:);
temp(dw_min_row:dw_max_row,:) = bin_img(dw_min_row:dw_max_row,:);
temp(:,round(imgCol-imgCol/5):imgCol) = ones(imgRow,round(imgCol/5)+1);
% subplot(1, 2, 1);imshow(temp)
% subplot(1, 2, 2);imshow(bin_img)
%% 开运算二值化图像去噪
% imshow(temp)
se = strel('rectangle',[2,2]);
denoiseImg = imerode(imdilate(temp,se),se);
% imshow(denoiseImg)
temp1 = imdilate(imerode(denoiseImg,strel('rectangle',[40,1])),strel('rectangle',[40,1]));
% imshow(temp1)
temp2 = denoiseImg - temp1;
% imshow(temp2)
temp3 = imdilate(imerode(temp2,strel('rectangle',[2,2])),strel('rectangle',[2,2]));
% temp3 = imerode(temp2,strel('rectangle',[1,2]));
% imshow(temp3) 
temp4 = imerode(imdilate(temp3,strel('rectangle',[10,1])),strel('rectangle',[10,1]));
% imshow(temp4) 
%% 连通域标记
L = bwlabel(temp4,8);
num = max(max(L));

imshow(gray_img);
hold on
% [r, c] = find(L==1);
% rc = [r c];
% y = min(rc(:,1)) - 5;
% x = min(rc(:,2)) - 5;
% h = max(rc(:,1)) - min(rc(:,1)) + 10;
% w = max(rc(:,2)) - min(rc(:,2)) + 10;
% rectangle('Position',[x,y,w,h],'EdgeColor','r','LineWidth',1)
for i = 1 : num
    [r, c] = find(L==i);
    rc = [r c];
    y = min(rc(:,1)) - 2;
    x = min(rc(:,2)) - 2;
    h = max(rc(:,1)) - min(rc(:,1)) + 4;
    w = max(rc(:,2)) - min(rc(:,2)) + 4;
    rectangle('Position',[x,y,w,h],'EdgeColor','r','LineWidth',1)
    hold on
end