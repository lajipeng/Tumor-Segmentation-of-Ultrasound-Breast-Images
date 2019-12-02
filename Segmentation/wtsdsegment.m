% WTSDSEGMENT Semiautomatic segmentation of breast lesions using watershed transform.
%   S = WTSDSEGMENT(I) computes the lesion segmentation using the watershed transform, 
%   where I is the breast ultrasound image. The constraint Gaussian variances are 
%   introduced manually by marking four points to indicate the lesion limits.
%
%   S = WTSDSEGMENT(I,B) B is a reference binary image used to calculate
%   the constraint Gaussian variances. It could be the radiologist manual
%   outlining. I and B must be of the same size.
%
%   Example:
%   -------
%   load('BUS01.mat');
%   S = wtsdsegment(I,Smanual);
%   delintumor(I,Smanual); title('Radiologist outlining')
%   delintumor(I,S); title('Computerized segmentation')
%
%   See also AUTOSEGMENT DELINTUMOR HORSCH
%
%
%   Reference:
%   ---------
%   W. Gomez, L. Leija, A. V. Alvarenga, A. F. C. Infantosi, W. C. A. Pereira, 
%   "Computerized lesion segmentation of breast ultrasound based on 
%   marker-controlled watershed transformation," Medical Physics,
%   vol. 37, no. 1, pp: 82-95, 2010.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (BraZil)
%   WTSDSEGMENT Version 3.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function BW2 = wtsdsegment(I,B)
if nargin < 2
    B = [];
end
% Image preprocessing
[I0,Gstr,rect] = getROI(I,B);
[Mi,Ni] = size(I0);
% Despeckle and contrast enhacement
%   You can put here your own desplecking and contrast
%   enhacement methods
I1a = fuzzyenh(I0);
I1  = isf(I1a,15);
% Radial gradients
G2 = radgrad(I1,Gstr);	
% Constraint Gaussian function
I3 = gaussfun(I1,Gstr);
% ARD loop
ARD = zeros(1,63);
BWs = false(Mi,Ni,63);
for t = 1:63
      BW = imfill(I3>=t,'holes');
      BW = mcwt(G2,BW,Gstr);
      ind = bwperim(BW);
      ARD(t) = mean(G2(ind));
      BWs(:,:,t) = BW;
end
% Best margin
[~,ind] = nanmax(ARD);
BW = BWs(:,:,ind);
BW = imfill(BW,'holes');
% Return to original image size
BW2 = false(size(I));
xw2min = rect(1);
xw2max = rect(2);
yh2min = rect(3);
yh2max = rect(4);
BW2(yh2min:yh2max,xw2min:xw2max) = BW;
BW2 = smoothboundary(BW2);
%*************************************************************
function BW2 = mcwt(grI,BW,Str)
se1 = Str.strel1;
se2 = Str.strel2;
xc = Str.xc;
yc = Str.yc;
mint = imerode(BW,se1);
mext = ~imdilate(BW,se2); 
mrk = mint|mext;
grI2 = imimposemin(grI,mrk);
C = watershed(grI2);
% Segmentation
if numel(unique(C(:))) > 3
    BW2 = ones(size(C));
else
    BW2 = C==C(yc,xc);
end
%*************************************************************
function F2 = radgrad(I1,Gstr)
I1 = double(I1);
[Y,X] = size(I1);
[gx,gy] = gradient(I1);
[x,y] = meshgrid(1:X,1:Y);
x = x - Gstr.xc;
y = y - Gstr.yc;
n = sqrt(x.^2+y.^2);
xu = x./n;
yu = y./n;
F1 = gx.*xu + gy.*yu;
F2 = sqrt(abs(F1));
id1 = isnan(F1);
F2(id1) = 0;
F2 = F2/max(F2(:));
ind = bwperim(ones(Y,X));
F2(ind) = 0;
%***********************************************************************
function [I2,G,rect] = getROI(I,B)
[yo,xo] = size(I);
if isempty(B)
    sel = 'N';
    while sel == 'N'
        h1 = figure('Name','Mark lesion limits',...
                    'NumberTitle','off',...
                    'WindowStyle','modal',...
                    'Resize','on',...
                    'Toolbar','none',...
                    'Menubar','none');
        imshow(I,[]);
        [xx, yy] = ginput(4);
        hold on; 
        line([min(xx) min(xx)],[1 yo],'Color','r','LineWidth',2,'LineStyle',':');
        line([max(xx) max(xx)],[1 yo],'Color','r','LineWidth',2,'LineStyle',':');
        line([1 xo],[min(yy) min(yy)],'Color','r','LineWidth',2,'LineStyle',':');
        line([1 xo],[max(yy) max(yy)],'Color','r','LineWidth',2,'LineStyle',':');
        pause(0.5);
        selection = questdlg('The lesion limits are correct?',...
                             'Lesion limits',...
                             'Yes','No','Yes');
        if strcmp(selection,'Yes')
            sel = 'Y';
        else
            sel = 'N';
        end
        close(h1); pause(0.1);
    end
else
    [yy,xx] = find(B);
end
    
xx = round(xx);
yy = round(yy);
xmin = min(xx);
ymin = min(yy);
ymax = max(yy);
xmax = max(xx);
w = abs(xmax - xmin); % Width
h = abs(ymax - ymin); % Heigth
w2 = fix(w/4);
h2 = fix(h/4);
xw2min = xmin-w2;
xw2max = xmax+w2;
yh2min = ymin-h2;
yh2max = ymax+h2;
if xw2min < 1
    xw2min = 1;
    lx = xmin;
else
    lx = w2;
end
if xw2max > xo
    xw2max = xo;
    ux = xo-xmax;
else
    ux = w2;
end
if yh2min < 1
    yh2min = 1;
    ly = ymin;
else
    ly = h2;
end
if yh2max > yo
    yh2max = yo;
    uy = yo-ymax;
else
    uy = h2;
end

I2 = I(yh2min:yh2max,xw2min:xw2max);
rect = [xw2min xw2max yh2min yh2max];

[nyo,nxo] = size(I2);
nxmin = 1+lx;
nxmax = nxo-ux;
nymin = 1+ly;
nymax = nyo-uy;
xc = fix((nxmax+nxmin)/2);
yc = fix((nymax+nymin)/2);
G.width = w;
G.height = h;
G.xc = xc;
G.yc = yc;
G.strel1 = strel('disk',3,0);
G.strel2 = strel('disk',9,0);
%***********************************************************************
function I2 = gaussfun(I1,str)
% Crea funcion Gaussiana delimitadora
[Mi,Ni] = size(I1);
w = str.width; % ancho
h = str.height; % alto
xc = str.xc;
yc = str.yc;
K = [(w/2) 0;0 (h/2)].^2;
[x,y] = meshgrid(1:Ni,1:Mi);
G = exp(-0.5*(((x-xc).^2)/K(1,1) + ((y-yc).^2)/K(2,2)))/(2*pi*sqrt(det(K)));
G = (G-min(G(:)))/(max(G(:))-min(G(:))); % Funcion Gaussiana
In = 1-(I1/max(I1(:)));
In = In.*G;
I2 = fix(64*(In-min(In(:)))/(max(In(:))-min(In(:))));
%***********************************************************************
function out = smoothboundary(S0)
BW = imfill(S0,'holes');
B = bwboundaries(BW);
x = B{1}(:,1);
y = B{1}(:,2);
S = x+1j*y;
F = fft(S);
N = numel(F);
ZP = floor(N/2)+1;
FShift = fftshift(F); 
numF = round(0.05*N);
FZeros = zeros(1,N);
FZeros((ZP-numF):(ZP+numF)) = FShift((ZP-numF):(ZP+numF));
G = ifftshift(FZeros);
K = ifft(G);
x = real(K);
y = imag(K);
out = roipoly(S0,y,x);