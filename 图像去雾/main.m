clc;
clear
addpath('./test images');
img = cell(1,5);
for i = 1:length(img)
    img{i} = imread(strcat('haze',num2str(i),'.jpg'));
end
% imshow(img{1})
%% 暗通道
imgA = img{5};
% imgA = imread('land.jpg');
[row,col,channel] = size(imgA);
% darkImg = zeros(row,col);
minChannel = zeros(row,col);
% 取最小通道值
for i = 1 : row
    for j = 1 : col
        minChannel(i,j) = min(imgA(i,j,:));
    end
end
% 最小值滤波,滤波前先padding复制边缘
temSize = 15; % 模板size值为单数
region = ones(temSize,temSize); % 滤波模板区域
pad_minChannel = padarray(minChannel,[floor(temSize/2) floor(temSize/2)],'replicate'); % padding
pad_darkImg = ordfilt2(pad_minChannel,1,region); % 最小值滤波
darkImg = pad_darkImg(1+floor(temSize/2):floor(temSize/2)+row,1+floor(temSize/2):floor(temSize/2)+col); % reverse padding
% imshow(uint8(darkImg)) 
%% 计算大气光和透射率
pob = 0.001;
maxNum = round(pob * row * col); % 暗通道最大的0.1%个元素
sortEle = sort(reshape(darkImg,1,[]),'descend'); 
[max_row,max_col] = find(darkImg >= sortEle(maxNum)); % 大于标定制的元素位置
% 取原图中这些位置最亮的像素最为大气光
grayA = rgb2gray(imgA); 
[atm_row,atm_col] = find(grayA == max(diag(grayA(max_row,max_col)))); %计算亮度最大值的位置
atmsphy = double(imgA(atm_row(1),atm_col(1),:)); % 大气光：RGB格式
%% 计算透射率图tx
w0 = 0.9;                    
% t0 = 0.1;
t_x = 1 - w0 * darkImg/min(atmsphy);
% imshow(uint8(t_x))
%% 对tx导向滤波
p = t_x;
I = t_x;
r = 3;
eps = 0.1;

[hei, wid] = size(p);
N = boxfilter(ones(hei, wid), r); 
mean_I = boxfilter(I, r) ./ N;
mean_p = boxfilter(p, r) ./ N;
mean_Ip = boxfilter(I.*p, r) ./ N;
% this is the covariance of (I, p) in each local patch.
cov_Ip = mean_Ip - mean_I .* mean_p; 
mean_II = boxfilter(I.*I, r) ./ N;
var_I = mean_II - mean_I .* mean_I;
a = cov_Ip ./ (var_I + eps); 
b = mean_p - a .* mean_I; 
mean_a = boxfilter(a, r) ./ N;
mean_b = boxfilter(b, r) ./ N;
q = mean_a .* I + mean_b;
fil_tx = q;
% imshow(uint8(fil_tx))
%% 计算去雾图像
t0 = 0.01;
fil_tx(find(fil_tx < t0)) = t0;
hazeRemov = (double(imgA) - atmsphy.*ones(row,col,channel))./fil_tx + atmsphy.*ones(row,col,channel);
subplot(1,2,1); 
imshow(uint8(imgA))
subplot(1,2,2); 
imshow(uint8(hazeRemov))