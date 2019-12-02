% CLAHE Contrast-limited adaptive histogram equalization.
%   J = CLAHE(I,TILES,CLIP) performs the contrast-limited adaptive histogram 
%   equalization to enhance the contrast of an image I, where TILES is a
%   1x2 vector indicating the number of contextual regions in x and y axes,
%   and CLIP is a constant in the range [0,1] to crop the histogram.
%
%   Example:
%   -------
%   I = imread('BUSreal.tif');
%   J = clahe(I,[2 2],0.5);
%   figure; 
%   subplot 121; imshow(I); title('Original Image');
%   subplot 122; imshow(J); title('Contrasted Image');
%
%   See also FUZZYENH HISTEQU SACE SIGMOIDFILT   
%
%
%   Reference:
%   ----------
%   K. Zuiderveld, "Contrast limited adaptive histograph equalization," 
%   in: P.S. Heckbert (Ed.), Graphic Gems IV, Academic Press, 
%   San Diego, USA, 1994, pp. 474-485.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   CLAHE Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function J = clahe(I,tiles,clip)
% Genera indices de cada elemento de la reticula
L = 256;
nx = tiles(1);
ny = tiles(2);
I = double(I);
[Y,X] = size(I);
sizex = fix(X/nx);
sizey = fix(Y/ny);
sx = 1:sizex:X; sx(nx+1) = X;
sy = 1:sizey:Y; sy(ny+1) = Y;
% Mejoramiento adaptativo
Mapa = zeros(ny,nx,256);
for i = 1:nx
   for j = 1:ny
       K = I(sy(j):sy(j+1)-1,sx(i):sx(i+1)-1);  % Region contextual
       Hg = accumarray(K(:)+1,ones(numel(K),1),[L 1],@sum,0);  % Histograma de la region
       Hc = cliphistogram(Hg,clip,L);
       fa = histequal(Hc,L);
       Mapa(j,i,:) = fa;       
   end
end
% interpolacion bilineal
J = bilinear(I,Mapa,nx,ny);
%------------------------------------------------------------------------
function fa = histequal(Hc,L)
Pf = cumsum(Hc/sum(Hc));
fa = (L-1)*Pf;
%------------------------------------------------------------------------
function Hc = cliphistogram(H,clip,L)
c = clip*max(H);
Hc = H;
idx = H>=c;
Hc(idx) = c;
Hc = Hc + ((sum(H)-sum(Hc))/L);
%------------------------------------------------------------------------
function CEImg = bilinear(I,Mapa,nx,ny)
[Y,X] = size(I);
sizex = fix(X/nx);
sizey = fix(Y/ny);
CEImg = zeros(Y,X);
xI = 1;
for i = 1:nx+1
    if i == 1
        subX = fix(sizex/2);
        xU = 1;
        xB = 1;
    elseif i == nx+1
        subX = X-xI+1;
        xU = nx;
        xB = nx;
    else
        subX = sizex;
        xU = i - 1;
        xB = i;
    end
    yI = 1;
    for j = 1:ny+1
        if j == 1
            subY = fix(sizey/2);
            yL = 1;
            yR = 1;
        elseif j == ny+1
            subY = Y-yI+1;
            yL = ny;
            yR = ny;
        else
            subY = sizey;
            yL = j - 1;
            yR = j;
        end
        UL = Mapa(yL,xU,:);
        UR = Mapa(yR,xU,:);
        BL = Mapa(yL,xB,:);
        BR = Mapa(yR,xB,:);
        subImg = I(yI:yI+subY-1,xI:xI+subX-1);
        subImg = interpolate(subImg,UL,UR,BL,BR,subX,subY);
        CEImg(yI:yI+subY-1,xI:xI+subX-1) = subImg;
        yI = yI + subY;
    end
    xI = xI + subX;
end
CEImg = uint8(CEImg);
%-----------------------------------------------------------------------
function [subImage] = interpolate(subBin,LU,LB,RU,RB,XSize,YSize)
subImage = zeros(size(subBin));
num = XSize * YSize;
for i = 0:YSize-1
    inverseI = YSize - i;
    for j = 0:XSize-1
        inverseJ = XSize - j;
        val = subBin(i+1,j+1)+1;
        subImage(i+1,j+1) = (inverseI*(inverseJ*LU(val) + j*RU(val)) + i*(inverseJ*LB(val) + j*RB(val)))/num;
    end
end