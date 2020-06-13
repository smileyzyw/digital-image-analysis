function imDst = boxfilter(imSrc,r)
% w,h�ֱ�Ϊ�����˲����ĺ��Ӱ뾶
% w must <= (wid-1)/2
% h must <= (hei-1)/2
w = r;
h = r;
[hei, wid] = size(imSrc);
imDst = zeros(size(imSrc));
%Y�᷽����ۼ�
imCum = cumsum(imSrc, 1);
%���ȿ����ײ���H������
imDst(1:h+1, :) = imCum(1+h:2*h+1, :);
%�м�����
imDst(h+2:hei-h, :) = imCum(2*h+2:hei, :) - imCum(1:hei-2*h-1, :);
%β��h������
imDst(hei-h+1:hei, :) = repmat(imCum(hei, :), [h, 1]) - imCum(hei-2*h:hei-h-1, :);
%X�᷽����ۼ�
imCum = cumsum(imDst, 2);
%���ȿ����ײ��ģظ�����
imDst(:, 1:w+1) = imCum(:, 1+w:2*w+1);
%�����м�����
imDst(:, w+2:wid-w) = imCum(:, 2*w+2:wid) - imCum(:, 1:wid-2*w-1);
%����β��w������
imDst(:, wid-w+1:wid) = repmat(imCum(:, wid), [1, w]) - imCum(:, wid-2*w:wid-w-1);
end