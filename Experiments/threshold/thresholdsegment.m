%-- Author£ºWang, Peng. Fudan University 
%-- Date£º2019/11/28
%-- Fuction£ºBreast Tumor Segmentation of Ultra Images based in threshold
%-- Reference:
%   https://github.com/FubuFabian/Breast-Tumor-Segmentation-Matlab
%   https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1194626&tag=1
function SEG_BW = thresholdsegment(I)
    imGray = mat2gray(I);
 %%%%%%%% Preprocessing Image %%%%%%%%%%%%%
    imEq = histeq(imGray);
    imDiff = SRAD3D(imEq,25,1);
    imDiff = mat2gray(imDiff);
 %%%%%%%% Segmenation based in threshold%%%%%%%%%%%%%
    I= imDiff;
    [a,b]=size(I);
    Io=I;
    n=0.35;
    Io = imcomplement(im2bw(Io,n));
%     for p=1:1:a-1
%         for q=1:1:b-1
%             if (I(p,q)>n)
%                 Io(p,q)=0;
%             else
%                 Io(p,q)=1;
%             end
%         end
%     end
%%%%%%%% Extract the contour of tumor %%%%%%%%%
    count = size(Io);
    se = strel('square',2);
    img_contour = im2bw(Io);
    bw = imerode(img_contour,se);
    bw = imerode(bw,se);
    bw = imerode(bw,se);
    contour1 = bwperim(bw,8);
    contour2 = imfill(contour1,'holes');
    contour3 = bwareaopen(contour2,200,8);
    imLabel = bwlabel(contour3);                
    stats = regionprops(imLabel,'Area');    
    area = cat(1,stats.Area);
    index = find(area == max(area));        
    img = ismember(imLabel,index);          
    SEG_BW = logical(imfill(img,'holes'));
end