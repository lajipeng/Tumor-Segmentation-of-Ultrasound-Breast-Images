% ADF Anisotropic diffusion filtering.
%   J = ADF(I,KAPPA,DT,MAXITER,OPTION) applies the anisotropic diffusion 
%   filtering on a gray scale image I, where KAPPA is a constant that controls 
%   the conduction, DT is the time step constant, MAXITER is the maximum number 
%   of iterations for diffusion, and OPTION selects between the two conduction 
%   functions proposed by Perona & Malik:
%
%      1 - c(x,y,t) = exp(-(nablaI/kappa).^2),
%          privileges high-contrast edges over low-contrast ones. 
%      2 - c(x,y,t) = 1./(1 + (nablaI/kappa).^2),
%          privileges wide regions over smaller ones. 
% 
%   Example:
%   -------
%   I = imread('BUSreal.tif');
%   J = adf(I,3,0.25,500,2);
%   figure; 
%   subplot 121; imshow(I); title('Original Image');
%   subplot 122; imshow(mat2gray(J)); title('Filtered Image');
%
%   See also ADLG CHMF ISF ISFAD
%
%
%   Reference:
%   ----------
%   P. Perona, J. Malik, "Scale-Space and Edge Detection Using Anisotropic 
%   Diffusion," IEEE Transactions on Pattern Analysis and Machine Intelligence, 
%   vol. 12, no. 7, pp. 629-639, 1990.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   ADF Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2014, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function J = adf(J,kappa,dt,iter,opt)
J = double(J);
[Mi,Ni] = size(J);
yp = [1, 1:Mi-1];
ym = [2:Mi, Mi];
xp = [1, 1:Ni-1];
xm = [2:Ni, Ni];
for t = 1:iter
    [Gxp,Gxm,Gyp,Gym] = grads(J,xp,xm,yp,ym);
    cGxp = diffcoef(Gxp,kappa,opt);
    cGxm = diffcoef(Gxm,kappa,opt);
    cGyp = diffcoef(Gyp,kappa,opt);
    cGym = diffcoef(Gym,kappa,opt);
    J = J + (dt/4)*(cGxp.*Gxp + cGyp.*Gyp + cGxm.*Gxm + cGym.*Gym);
end
%---------------------------------------------------------------------
function [Gxp,Gxm,Gyp,Gym] = grads(I,xp,xm,yp,ym)
Gxp = I(:,xp) - I;
Gxm = I(:,xm) - I;
Gyp = I(yp,:) - I;
Gym = I(ym,:) - I;
%---------------------------------------------------------------------
function c = diffcoef(I,kappa,opt)
if opt == 1
    c = exp(-(abs(I)/kappa).^2);
elseif opt == 2
    c = 1./(1+((abs(I)/kappa).^2));
end