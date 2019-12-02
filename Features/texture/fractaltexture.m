% FRACTALTEXTURE Fractal features of gray-level image.
%   [X,FEATS] = FRACTALTEXTURE(I,BW) computes eight fractal features
%   from the gray-level image I masked by the binary image BW: fractal dimension
%   (FD, first value) and fractal Brownian motion (FBM, 2nd to 8th values). 
%   X is a numeric vector with the feature values and FEATS is a cell vector 
%   with the name of the features in the same order as in X.
%
%   Example:
%   -------
%   load('BUS01.mat');   
%   [x,feats] = fractaltexture(I,Smanual);
%
%   See also CLXCURVE HISTFEATURES BGC
%
%
%   Reference:
%   ---------
%   D.-R. Chen, R.-F. Chang, C.-J. Chen, et al., "Classification of breast 
%   ultrasound images using fractal feature," Clin. Imaging, vol. 29, no. 4,
%   pp. 235-245, 2005.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   FRACTALTEXTURE Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = fractaltexture(US,BW)
se = strel('ball',3,3);
US = imclose(US,se);
US = histequ(US);
[y,x] = find(BW);
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
I = double(US(ymin:ymax,xmin:xmax));
[M,N] = size(I);
K = 8;
i = 1:K;
id = zeros(1,K);
for k = i
    % Primera diferencia: Vertical
    yp = [ones(1,k) 1:M-k];
    D1 = abs(I-I(yp,:))/(M*(M-k));
    % Segunda diferencia: Horizontal
    xp = [ones(1,k) 1:N-k];
    D2 = abs(I-I(:,xp))/(N*(N-k));
    % Tercera diferencia: Diagonal
    dp = round(k/sqrt(2));
    xp = [ones(1,dp) 1:N-dp];
    yp = [ones(1,dp) 1:M-dp];
    D3 = abs(I-I(yp,xp))/((M-dp)*(N-dp));
    % Cuarta diferencia: Diagonal asimetrica
    xp = [dp+1:N N*ones(1,dp)];
    yp = [ones(1,dp) 1:M-dp]; 
    D4 = abs(I-I(yp,xp))/((M-dp)*(N-dp));
    % Calcula id(k)
    id(k) = sum(sum(D1+D2+D3+D4))/4;
end
% Salidas
y = log(id)-log(id(1));
p = polyfit(log(i),y,1);
y(1) = [];
x = [p(1) y];
feats = {'FD','FBM1','FBM2','FBM3','FBM4','FBM5','FBM6','FBM7'};