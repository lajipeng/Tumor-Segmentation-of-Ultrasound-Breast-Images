% SACE Adaptive contrast enhancement based on sigmoidal mapping function.
%   J = SACE(I,TILES,OPTION) performs the adaptive contrast enhancement 
%   based on sigmoidal mapping function to enhance the contrast of a breast
%   ultrasound image I, where TILES is a 1x2 vector indicating the number of 
%   contextual regions in x and y axes, and OPTION is a flag parameter such
%   that if OPTION=1 the output image is not filtered and if OPTION=2 the
%   output image is despeckled by the circular hybrid median filter.
%
%   NOTE: The functions "estpab" and "estmutualinfo" were compiled from the 
%   mutual information computation, conditional probability and entropy 
%   estimation for discrete/categorical variables package by Hanchuan Peng
%   at http://penglab.janelia.org/proj/mRMR/
%
%   Example:
%   -------
%   I = imread('BUSreal.tif');
%   J1 = sace(I,[2 2],1);
%   J2 = sace(I,[2 2],2);
%   figure; 
%   subplot 131; imshow(I); title('Original Image');
%   subplot 132; imshow(mat2gray(J1)); title('Contrasted Image');
%   subplot 133; imshow(mat2gray(J2)); title('Filtered Image');
%
%   See also CLAHE FUZZYENH HISTEQU SIGMOIDFILT  
%
%
%   Reference:
%   ----------
%   W. Gomez, W. C. A. Pereira, "A contrast enhancement method for 
%   improving the segmentation of breast lesions on ultrasonography," 
%   Computers in Medicine and Biology, vol.80, pp. 14-23, 2017.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   SACE Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function J = sace(go,tiles,opt)
% Parametros de entrada
nx = tiles(1);
ny = tiles(2);
L = 256;
go = double(go);
% Intensity variation map
g = ivm(go);
[Y,X] = size(g);
sizex = fix(X/nx);
sizey = fix(Y/ny);
sx = 1:sizex:X; sx(nx+1) = X;
sy = 1:sizey:Y; sy(ny+1) = Y;
% Procesamiento por regiones contextuales
MapSig = zeros(ny,nx,L);
for i = 1:length(sx)-1
   for j = 1:length(sy)-1 
        K = g(sy(j):sy(j+1)-1,sx(i):sx(i+1)-1);  % Region contextual
        [y,hg] = maxvar(K,L);
        chg = cumsum(hg/sum(hg));
        x = find(chg>=0.05,1,'first');
        z = find(chg<=0.95,1,'last');
        S = zeros(L,1);
        gf = linspace(0,L-1,L);
        % Rangos de la funcion
        ind1 = gf<=x;
        ind2 = (x<gf)&(gf<=y);
        ind3 = (y<gf)&(gf<=z);
        ind4 = z<gf;
        % Funcion sigmoidal
        S(ind1) = 0;
        S(ind2) = ((gf(ind2)-x).^2)./((y-x)*(z-x));
        S(ind3) = 1-((gf(ind3)-z).^2)./((z-y)*(z-x));
        S(ind4) = 1;
        MapSig(j,i,:) = floor((L-1)*S);
   end
end
% Interpolacion bilineal
if opt==1
    J = bilinear(go,MapSig,nx,ny);
elseif opt==2
    J = bilinear(g,MapSig,nx,ny);
end
J = fix(J);
%**********************************************************************
% Calcula maxima varianza (Otsu)
function [T,yi,mx] = maxvar(K,L)
Nk = numel(K(:));
yi = accumarray(K(:)+1,ones(Nk,1),[L 1],@sum,0);  % Histograma de la region
Aj = cumsum(yi)+eps;
Bj = cumsum((0:L-1)'.*yi)+eps;
An = Aj(end); Bn = Bj(end);
mu_j = Bj./Aj;
nu_j = (Bn-Bj)./((An-Aj)+eps);
sigma2 = Aj.*(An-Aj).*(mu_j-nu_j).^2;
[mx,ind] = max(sigma2);
T = ind-1;
%**********************************************************************
% Interpolacion bilineal segun CLAHE
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
%-----------------------------------------------------------------------
function [subImage] = interpolate(subBin,LU,LB,RU,RB,XSize,YSize)
subImage = zeros(size(subBin));
num = XSize * YSize;
for i = 0:YSize-1
    inverseI = YSize - i;
    for j = 0:XSize-1
        inverseJ = XSize - j;
        val = subBin(i+1,j+1)+1;
        subImage(i+1,j+1) = fix((inverseI*(inverseJ*LU(val) + j*RU(val)) + i*(inverseJ*LB(val) + j*RB(val)))/num);
    end
end
%***********************************************************************
function [actI,w] = ivm(I)
e = 0.01;            % Umbral de minima diferencia
w = 5;
preMI = Inf;
flag = 1;
I = double(I);
[Mi,Ni] = size(I);
mn = min(Mi/2,Ni/2);
preI = double(I);
while flag
    actI = isf(I,w);
    MI = mutualinfo(preI,actI);
    preI = actI;
    if abs(MI-preMI) < e % Verifica condicion de paro
       flag = 0;
    else
        w = w+2;
        preMI = MI;
    end
    if w > mn
        flag = 0;
        w = w-2;
    end
end
%***********************************************************************
function h = mutualinfo(vec1,vec2)
vec1 = vec1(:);
vec2 = vec2(:);
[p12, p1, p2] = estpab(vec1,vec2);
h = estmutualinfo(p12,p1,p2);