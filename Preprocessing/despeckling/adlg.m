% ADLG Anisotropic diffusion guided by Log-Gabor filters.
%   J = ADLG(I,DT,MAXITER) performs the anisotropic difussion
%   filtering guided by Log-Gabor filters for despeckling a breast 
%   ultrasound image I, where DT is the time step constant, and MAXITER is 
%   the maximum number of iterations for diffusion.
%
%   Example:
%   -------
%   I = imread('BUSreal.tif');
%   J = adlg(I,0.25,1000);
%   figure; 
%   subplot 121; imshow(I); title('Original Image');
%   subplot 122; imshow(mat2gray(J)); title('Filtered Image');
%
%   See also ADF CHMF ISF ISFAD
%
%
%   Reference:
%   ----------
%   W. Gomez, W. C. A. Pereira, A. F. C. Infantosi, "Breast ultrasound 
%   despeckling using anisotropic diffusion guided by texture descriptors," 
%   Ultrasound in Medicine and Biology, vol. 40, no. 11, pp. 2609-2621, 2014.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   ADLF Version 1.0 (Matlab R2012b Unix)
%   April 2014
%   Copyright (c) 2014, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function I = adlg(I,dt,niter)
% Validations
if nargin < 2
    dt = 0.25;
    niter = 1000;
end
if nargin < 3
    niter = 1000;
end
% Smoothing to homogeneize similar textures
I = hybmedfilt2(I,11);
% Normalizacion to the range [0,1]
I = imnormalize(I); 
% Creating the bank of Log-Gabor filters (Equation 7)
nscale  = 3;  % Scales
norient = 24; % Orientations
LGab = LogGabor(I, nscale, norient, 1, 3.0, .55); % 2 octaves
% Filtering the input image to extract texture information (Equation 8)
T = gabortextures(I,LGab);
% Texture-gradient responses (Equation 2)
GT = gradients(T);
% Texture-gradient magnitudes, adaptive kappa in N8 (Equation 10)
RT = magnitudes(GT);
% Anisotropic diffusion process
Iact = round(255*I);
for t = 1:niter
    % Intensity-gradient responses (Equation 2)
    GI = gradients(I);
    % Diffusion coefficient (Equation 9)
    C = diffcoeffs(GI,RT);
    % Divergencies
    D = divergences(GI,C);
    % Solution update (Equation 12)
    I = updates(I,D,dt);
    % Stop criterion by measuring MSE
    if ~mod(t,10)
        Iant = Iact;
        Iact = round(255*I);
        err = mserr(Iact,Iant); % MSE (Equation 16)
        if err<0.5
            disp(['Iteration stop: ' num2str(t) ' with MSE: ' num2str(err)]);
            break;
        end
    end
end
%*************************************************************************
%*************************************************************************
function OUT = updates(I,IN1,dt)
% Divergences in N8
D_W_N = IN1.D_W_N;
D_N_E = IN1.D_N_E;
D_E_S = IN1.D_E_S;
D_S_W = IN1.D_S_W;
D_NW_NE = IN1.D_NW_NE;
D_NE_SE = IN1.D_NE_SE;
D_SE_SW = IN1.D_SE_SW;
D_SW_NW = IN1.D_SW_NW;
% Equation 12
J = I + dt*(D_W_N + D_N_E + D_E_S + D_S_W + D_NW_NE + D_NE_SE + D_SE_SW + D_SW_NW);
% Normalizacion to the range [0,1]
OUT = imnormalize(J);
%*************************************************************************
%*************************************************************************
function OUT = divergences(IN1,IN2)
% Intensity-gradient responses in N8
I_W_N(:,:,1) = IN1.W_N.D1; I_W_N(:,:,2) = IN1.W_N.D2;
I_N_E(:,:,1) = IN1.N_E.D1; I_N_E(:,:,2) = IN1.N_E.D2;
I_E_S(:,:,1) = IN1.E_S.D1; I_E_S(:,:,2) = IN1.E_S.D2;
I_S_W(:,:,1) = IN1.S_W.D1; I_S_W(:,:,2) = IN1.S_W.D2;
I_NW_NE(:,:,1) = IN1.NW_NE.D1; I_NW_NE(:,:,2) = IN1.NW_NE.D2;
I_NE_SE(:,:,1) = IN1.NE_SE.D1; I_NE_SE(:,:,2) = IN1.NE_SE.D2;
I_SE_SW(:,:,1) = IN1.SE_SW.D1; I_SE_SW(:,:,2) = IN1.SE_SW.D2;
I_SW_NW(:,:,1) = IN1.SW_NW.D1; I_SW_NW(:,:,2) = IN1.SW_NW.D2;
% Conduction coefficients in N8
C_W_N = IN2.C_W_N;
C_N_E = IN2.C_N_E;
C_E_S = IN2.C_E_S;
C_S_W = IN2.C_S_W;
C_NW_NE = IN2.C_NW_NE;
C_NE_SE = IN2.C_NE_SE;
C_SE_SW = IN2.C_SE_SW;
C_SW_NW = IN2.C_SW_NW;
% Computes divergences in N8
D_W_N = sum(C_W_N.*I_W_N,3);
D_N_E = sum(C_N_E.*I_N_E,3);
D_E_S = sum(C_E_S.*I_E_S,3);
D_S_W = sum(C_S_W.*I_S_W,3);
D_NW_NE = sum(C_NW_NE.*I_NW_NE,3);
D_NE_SE = sum(C_NE_SE.*I_NE_SE,3);
D_SE_SW = sum(C_SE_SW.*I_SE_SW,3);
D_SW_NW = sum(C_SW_NW.*I_SW_NW,3);
% Outputs
OUT.D_W_N = D_W_N;
OUT.D_N_E = D_N_E;
OUT.D_E_S = D_E_S;
OUT.D_S_W = D_S_W;
OUT.D_NW_NE = D_NW_NE;
OUT.D_NE_SE = D_NE_SE;
OUT.D_SE_SW = D_SE_SW;
OUT.D_SW_NW = D_SW_NW;
%*************************************************************************
%*************************************************************************
function OUT = diffcoeffs(IN1,IN2)
% Intensity-gradient responses in N8
I_W_N(:,:,1) = abs(IN1.W_N.D1).^2; I_W_N(:,:,2) = abs(IN1.W_N.D2).^2;
I_N_E(:,:,1) = abs(IN1.N_E.D1).^2; I_N_E(:,:,2) = abs(IN1.N_E.D2).^2;
I_E_S(:,:,1) = abs(IN1.E_S.D1).^2; I_E_S(:,:,2) = abs(IN1.E_S.D2).^2;
I_S_W(:,:,1) = abs(IN1.S_W.D1).^2; I_S_W(:,:,2) = abs(IN1.S_W.D2).^2;
I_NW_NE(:,:,1) = abs(IN1.NW_NE.D1).^2; I_NW_NE(:,:,2) = abs(IN1.NW_NE.D2).^2;
I_NE_SE(:,:,1) = abs(IN1.NE_SE.D1).^2; I_NE_SE(:,:,2) = abs(IN1.NE_SE.D2).^2;
I_SE_SW(:,:,1) = abs(IN1.SE_SW.D1).^2; I_SE_SW(:,:,2) = abs(IN1.SE_SW.D2).^2;
I_SW_NW(:,:,1) = abs(IN1.SW_NW.D1).^2; I_SW_NW(:,:,2) = abs(IN1.SW_NW.D2).^2;
% Adaptive kappa in N8
R_W_N = IN2.R_W_N;
R_N_E = IN2.R_N_E;
R_E_S = IN2.R_E_S;
R_S_W = IN2.R_S_W;
R_NW_NE = IN2.R_NW_NE;
R_NE_SE = IN2.R_NE_SE;
R_SE_SW = IN2.R_SE_SW;
R_SW_NW = IN2.R_SW_NW; 
% Computes diffusion coefficient in N8 (Equation 9)
s = 1/100;
C_W_N = s./(s+(I_W_N./R_W_N)); C_W_N(C_W_N<0) = 0; C_W_N(C_W_N>1)=1;
C_N_E = s./(s+(I_N_E./R_N_E)); C_N_E(C_N_E<0) = 0; C_N_E(C_N_E>1)=1;
C_E_S = s./(s+(I_E_S./R_E_S)); C_E_S(C_E_S<0) = 0; C_E_S(C_E_S>1)=1;
C_S_W = s./(s+(I_S_W./R_S_W)); C_S_W(C_S_W<0) = 0; C_S_W(C_S_W>1)=1; 
C_NW_NE = s./(s+(I_NW_NE./R_NW_NE)); C_NW_NE(C_NW_NE<0) = 0; C_NW_NE(C_NW_NE>1)=1;
C_NE_SE = s./(s+(I_NE_SE./R_NE_SE)); C_NE_SE(C_NE_SE<0) = 0; C_NE_SE(C_NE_SE>1)=1;
C_SE_SW = s./(s+(I_SE_SW./R_SE_SW)); C_SE_SW(C_SE_SW<0) = 0; C_SE_SW(C_SE_SW>1)=1;
C_SW_NW = s./(s+(I_SW_NW./R_SW_NW)); C_SW_NW(C_SW_NW<0) = 0; C_SW_NW(C_SW_NW>1)=1;
% Outputs
OUT.C_W_N = C_W_N;
OUT.C_N_E = C_N_E;
OUT.C_E_S = C_E_S;
OUT.C_S_W = C_S_W;
OUT.C_NW_NE = C_NW_NE;
OUT.C_NE_SE = C_NE_SE;
OUT.C_SE_SW = C_SE_SW;
OUT.C_SW_NW = C_SW_NW;
%*************************************************************************
%*************************************************************************
function OUT = gradients(I)
I = double(I);
% Indices for finite differences
[Mi,Ni,~] = size(I);
S = [1,1:Mi-1];
N = [2:Mi, Mi];
E = [1,1:Ni-1];
W = [2:Ni, Ni];
% Finite differences (gradients) in N8 (Equation 2)
W_N.D1 = I(:,W,:)-I; W_N.D2 = I(N,:,:)-I;
N_E.D1 = I(N,:,:)-I; N_E.D2 = I(:,E,:)-I;
E_S.D1 = I(:,E,:)-I; E_S.D2 = I(S,:,:)-I;
S_W.D1 = I(S,:,:)-I; S_W.D2 = I(:,W,:)-I;
NW_NE.D1 = I(N,W,:)-I; NW_NE.D2 = I(N,E,:)-I;
NE_SE.D1 = I(N,E,:)-I; NE_SE.D2 = I(S,E,:)-I;
SE_SW.D1 = I(S,E,:)-I; SE_SW.D2 = I(S,W,:)-I;
SW_NW.D1 = I(S,W,:)-I; SW_NW.D2 = I(N,W,:)-I;
% Outputs
OUT.W_N = W_N;
OUT.N_E = N_E;
OUT.E_S = E_S;
OUT.S_W = S_W;
OUT.NW_NE = NW_NE;
OUT.NE_SE = NE_SE;
OUT.SE_SW = SE_SW;
OUT.SW_NW = SW_NW;
%*************************************************************************
%*************************************************************************
function F = gabortextures(I,Gab)
[Mi, Ni, z] = size(Gab);
TF = fft2(I); 
TF(1,1) = 0;
% Filter the input images with Log-Gabor
F = zeros(Mi,Ni,z);
for i = 1:z
    F(:,:,i) = real(ifft2(TF.*Gab(:,:,i))); % Equation 8
end
%*************************************************************************
%*************************************************************************
function OUT = magnitudes(IN)
% Intensity-gradient responses in N8
W_N = IN.W_N;
N_E = IN.N_E;
E_S = IN.E_S;
S_W = IN.S_W;
NW_NE = IN.NW_NE;
NE_SE = IN.NE_SE;
SE_SW = IN.SE_SW;
SW_NW = IN.SW_NW;
% Magnitude responses in N8 (Equation 10)
R_W_N(:,:,1) = sqrt(sum(abs(W_N.D1).^2,3)); R_W_N(:,:,2) = sqrt(sum(abs(W_N.D2).^2,3));
R_N_E(:,:,1) = sqrt(sum(abs(N_E.D1).^2,3)); R_N_E(:,:,2) = sqrt(sum(abs(N_E.D2).^2,3));
R_E_S(:,:,1) = sqrt(sum(abs(E_S.D1).^2,3)); R_E_S(:,:,2) = sqrt(sum(abs(E_S.D2).^2,3));
R_S_W(:,:,1) = sqrt(sum(abs(S_W.D1).^2,3)); R_S_W(:,:,2) = sqrt(sum(abs(S_W.D2).^2,3));
R_NW_NE(:,:,1) = sqrt(sum(abs(NW_NE.D1.^2),3)); R_NW_NE(:,:,2) = sqrt(sum(abs(NW_NE.D2).^2,3));
R_NE_SE(:,:,1) = sqrt(sum(abs(NE_SE.D1.^2),3)); R_NE_SE(:,:,2) = sqrt(sum(abs(NE_SE.D2).^2,3));
R_SE_SW(:,:,1) = sqrt(sum(abs(SE_SW.D1.^2),3)); R_SE_SW(:,:,2) = sqrt(sum(abs(SE_SW.D2).^2,3));
R_SW_NW(:,:,1) = sqrt(sum(abs(SW_NW.D1.^2),3)); R_SW_NW(:,:,2) = sqrt(sum(abs(SW_NW.D2).^2,3));
% Outputs
OUT.R_W_N = R_W_N+eps;
OUT.R_N_E = R_N_E+eps;
OUT.R_E_S = R_E_S+eps;
OUT.R_S_W = R_S_W+eps;
OUT.R_NW_NE = R_NW_NE+eps;
OUT.R_NE_SE = R_NE_SE+eps;
OUT.R_SE_SW = R_SE_SW+eps;
OUT.R_SW_NW = R_SW_NW+eps;
%*************************************************************************
%*************************************************************************
function J = imnormalize(I)
I = double(I);
mn = min(I(:));
mx = max(I(:));
J = (I-mn)./(mx-mn); % Normalization [0,1]
%*************************************************************************
%*************************************************************************
% Kovesi's implementation of Log-Gabor filters
% Copyright (c) 2001-2010 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/~pk/research/matlabfns/PhaseCongruency/Docs/convexpl.html
function LG = LogGabor(IO, nscale, norient, minWaveLength, mult, sigmaOnf)
[rows, cols] = size(IO);
if mod(cols,2)
    xrange = (-(cols-1)/2:(cols-1)/2)/(cols-1);
else
    xrange = (-cols/2:(cols/2-1))/cols; 
end
if mod(rows,2)
    yrange = (-(rows-1)/2:(rows-1)/2)/(rows-1);
else
    yrange = (-rows/2:(rows/2-1))/rows; 
end
[x,y] = meshgrid(xrange, yrange);
radius = sqrt(x.^2 + y.^2);
theta = atan2(y,x);
radius = ifftshift(radius);
theta  = ifftshift(theta);
radius(1,1) = 1;
sintheta = sin(theta);
costheta = cos(theta);
clear x; clear y; clear theta;
lp = lowpassfilter(xrange,yrange,.4,10);
logGabor = cell(1,nscale);
for s = 1:nscale
    wavelength = minWaveLength*mult^(s-1);
    fo = 1.0/wavelength;
    logGabor{s} = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2));  
    logGabor{s} = logGabor{s}.*lp;
    logGabor{s}(1,1) = 0;
    L = sqrt(sum(logGabor{s}(:).^2));
    logGabor{s} = logGabor{s}./L;        
end
% The main loop
LG = zeros(rows, cols, norient*nscale);
k = 1;
for o = 1:norient
    angl = (o-1)*pi/norient;
    wavelength = minWaveLength;
    ds = sintheta * cos(angl) - costheta * sin(angl);
    dc = costheta * cos(angl) + sintheta * sin(angl);
    dtheta = abs(atan2(ds,dc));
    dtheta = min(dtheta*norient/2,pi);
    spread = (cos(dtheta)+1)/2;        
    for s = 1:nscale
        filter = logGabor{s} .* spread;
        L = sqrt(sum(real(filter(:)).^2 + imag(filter(:)).^2 ))/sqrt(2);
        LG(:,:,k) = filter./L;
        wavelength = wavelength * mult;
        k = k+1;
    end
end
%------------------------------------------------------------------------
function f = lowpassfilter(xrange, yrange, cutoff, n)
[x,y] = meshgrid(xrange, yrange);
radius = sqrt(x.^2 + y.^2);
f = ifftshift( 1.0 ./ (1.0 + (radius ./ cutoff).^(2*n)) );
%*************************************************************************
%*************************************************************************
function m = mserr(X,Y)
m = mean2((X-Y).^2); % Equation 16
%*************************************************************************
%*************************************************************************
% Hybrid median filter
function f = hybmedfilt2(g, m)
m2 = round(m/2);
D = eye(m)+fliplr(eye(m)); D(m2,m2) = 1;
R = zeros(m); R(m2,:)=1; R(:,m2) = 1;
f1 = ordfilt2(g, fix(median(1:sum(D(:)))), D, 'symmetric');
f2 = ordfilt2(g, fix(median(1:sum(R(:)))), R, 'symmetric');
M = cat(3,g,f1,f2);
f = uint8(median(double(M),3));