function ModelFilterCompare(subset1,subset2)
%-- Author：Wang, Peng. Fudan University 
%-- Date：2019/11/30
%-- Fuction：Compare the preprocessing among diffirent model
%% loading and initializating image
load('C:/Github/医学成像技术/project2/Tumor-Segmentation-of-Ultrasound-Breast-Images/Input.mat')
load('C:/Github/医学成像技术/project2/Tumor-Segmentation-of-Ultrasound-Breast-Images/GT.mat')
original_image = Input{subset1};original_image1 = Input{subset2};
%% Extract the contour of ground truth
GroundTrue = logical(imfill(GT{subset1},'holes'));
[B1,L1] = bwboundaries(GroundTrue,'noholes');
boundary1 = B1{1}; 
GroundTrue1 = logical(imfill(GT{subset2},'holes'));
[B2,L2] = bwboundaries(GroundTrue1,'noholes');
boundary2 = B2{1}; 
%% The prepocessing of autosegment
I0 = original_image;
I1 = clahe(I0,[2 2],0.5);
I2 = isf(I1,11);
autosegment = 1-I2/max(I2(:));
I0 = original_image1;
I1 = clahe(I0,[2 2],0.5);
I2 = isf(I1,11);
autosegment1 = 1-I2/max(I2(:));
%% The prepocessing of horsch
% Despeckle
[I0,Gstr,rect] = getROI(original_image,GroundTrue);
I1 = medfilt2(original_image,[10 10]);
I1 = double(I1);
% Radial gradients
G2 = radgrad(I1,Gstr);	
% Constraint Gaussian function
horschsegment = gaussfun(I1);

% Despeckle
[I0,Gstr,rect] = getROI(original_image1,GroundTrue1);
I1 = medfilt2(original_image,[10 10]);
I1 = double(I1);
% Radial gradients
G2 = radgrad(I1,Gstr);	
% Constraint Gaussian function
horschsegment1 = gaussfun(I1,Gstr);
%% The prepocessing of wtsdsegment
[I0,Gstr,rect] = getROI(original_image,GroundTrue);
I1a = fuzzyenh(original_image);
I1  = isf(I1a,15);
% Radial gradients
G2 = radgrad(I1,Gstr);	
% Constraint Gaussian function
wtsdsegment = gaussfun(I1,Gstr);

[I0,Gstr,rect] = getROI(original_image1,GroundTrue1);
I1a = fuzzyenh(original_image1);
I1  = isf(I1a,15);
% Radial gradients
G2 = radgrad(I1,Gstr);	
% Constraint Gaussian function
wtsdsegment1 = gaussfun(I1,Gstr);
%% The prepocessing of threshold
imGray = mat2gray(original_image);
imEq = histeq(imGray);
imDiff = SRAD3D(imEq,25,1);
threshold = mat2gray(imDiff);

imGray = mat2gray(original_image1);
imEq = histeq(imGray);
imDiff = SRAD3D(imEq,25,1);
threshold1 = mat2gray(imDiff);
%% imshow
figure;
subplot 251; imshow(original_image,'InitialMagnification', 'fit'); xlabel('(a)');
subplot 252; imshow(threshold,'InitialMagnification', 'fit'); xlabel('(b)');
hold on; 
plot(boundary1(:, 2), boundary1(:, 1), 'r', 'LineWidth', 0.5); 
subplot 253; imshow(horschsegment,[],'InitialMagnification', 'fit'); xlabel('(c)');
hold on; 
plot(boundary1(:, 2), boundary1(:, 1), 'r', 'LineWidth', 0.5); 
subplot 254; imshow(autosegment,'InitialMagnification', 'fit'); xlabel('(d)');
hold on; 
plot(boundary1(:, 2), boundary1(:, 1), 'r', 'LineWidth', 0.5); 
subplot 255; imshow(wtsdsegment,[],'InitialMagnification', 'fit'); xlabel('(e)');
hold on; 
plot(boundary1(:, 2), boundary1(:, 1), 'r', 'LineWidth', 0.5); 

subplot 256; imshow(original_image1,'InitialMagnification', 'fit'); xlabel('(f)');
subplot 257; imshow(threshold1,'InitialMagnification', 'fit'); xlabel('(g)');
hold on; 
plot(boundary2(:, 2), boundary2(:, 1), 'r', 'LineWidth', 0.5); 
subplot 258; imshow(horschsegment1,[],'InitialMagnification', 'fit'); xlabel('(h)');
hold on; 
plot(boundary2(:, 2), boundary2(:, 1), 'r', 'LineWidth', 0.5); 
subplot 259; imshow(autosegment1,'InitialMagnification', 'fit'); xlabel('(i)');
hold on; 
plot(boundary2(:, 2), boundary2(:, 1), 'r', 'LineWidth', 0.5); 
subplot(2,5,10); imshow(wtsdsegment1,[],'InitialMagnification', 'fit'); xlabel('(j)');
hold on; 
plot(boundary2(:, 2), boundary2(:, 1), 'r', 'LineWidth', 0.5); 
end
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
end
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
end