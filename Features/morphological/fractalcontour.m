% FRACTALCONTOUR Fractal dimension of the contour of a lesion.
%   [X,FEATS] = FRACTALCONTOUR(BW) computes two methods of the fractal 
%   dimension of the contour of a breast lesion BW: ruler and box methods. 
%   X is the numeric value of the features and FEATS is the name of the
%   features in the same order as in X.
%
%   [X,FEATS] = FRACTALCONTOUR(BW,FEATURES) computes computes a specific set 
%   of morphological features according with FEATURES. The set of valid 
%   strings includes (case insensitive):
%
%       'all'   - Entire feature set.
%       'ruler' - Ruler method.
%       'box1'   - Box method (d={1/4,...,1/128}).
%       'box2'   - Box method (d={1/5,...,1/160}).
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS02.mat');   
%   [x,feats] = fractalcontour(BW);
%   % Same as [x,feats] = fractalcontour(BW,'all')
%
%   Example 2: Compute the box method
%   ---------------------------------
%   load('BUS02.mat');   
%   [x,feats] = fractalcontour(BW,'box');
%
%   See also FOURIERFACTOR FOURIERSHAPE POLYMODEL
%
%
%   Reference:
%   ---------
%   R. M. Rangayyan, T. M. Nguyen, "Fractal analysis of contours
%   of breast masses in mammograms," Journal of Digital Imaging, 
%   vol. 20, no. 3, pp. 223-237, 2007.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   FRACTALCOUNTOUR Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Arturo Rodriguez Cristerna, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = fractalcontour(BW,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'ruler';'box1';'box2'};      
idxStats = RequestedInputs(officialStats,opts{:});
%*********************************************************************
Prps = regionprops(BW,'Centroid');
xc = Prps.Centroid(1);
yc = Prps.Centroid(2);
Baux = bwboundaries(BW);
C = Baux{1};
minCyz=min(C(:,1)); maxCyz=max(C(:,1));
minCxz=min(C(:,2)); maxCxz=max(C(:,2));
width=minCxz-maxCxz;
height=minCyz-maxCyz;
normBW=max([width, height]);
CDF= sqrt(((C(:,1)-yc)./normBW).^2+((C(:,2)-xc)./normBW).^2); %normalization in base to the maximum side of the lesion
CDF=CDF/max(CDF); %normalization of the radius
idx=1:size(C,1);
idx=idx./max(idx); %normalization of the lenght of the contour
%---- fractal compass (ruler method)
%size of compass
szBR=0.050; m=1;
while szBR<=0.2
    rsz(m)=szBR;
    m=m+1; szBR=szBR+0.025;
end
for msz=1:size(rsz,2) %counting how many compasses of size rsz(msz) are in the contour
    countrsz(msz)=0;
    n=1;
    while n<size(CDF,1)
        nnew=n;
        for m=n+1:size(CDF,1)
            if rsz(msz)>=sqrt( (idx(m)-idx(n))^2 + (CDF(m)-CDF(n))^2 )
                nnew=m-1;
            end
        end
        n=nnew+1;
        countrsz(msz)=countrsz(msz)+1;
    end
end
%compass (ruler) dimension
[pC,~] = polyfit(log(1./rsz),log(countrsz),1);
%---- fractal box (box counting method) paper
aux1 = fractal_box(CDF,idx,4,128);
aux2 = fractal_box(CDF,idx,5,160);
% Features
aux_x=[pC(1), aux1, aux2];
aux_f = {'fractalD_ruler' 'fractalD_box1' 'fractalD_box2'};
x = aux_x(idxStats);
feats = aux_f(idxStats);

%*********************************************************************
function aux = fractal_box(CDF,idx,szBBo,szBBo2)
m=1;
while szBBo<=szBBo2
    bsz(m)=1/szBBo;
    m=m+1; szBBo=szBBo*2;
end
for msz=1:size(bsz,2)
    qL=1/bsz(msz); %quantization level to be used in both axes    
    yQ=round(CDF*qL);%quantizing the centroid distance function
    xQ=round(idx*qL);%quantizing the contour points
    for m=1:size(CDF,1)
        crM(m)=yQ(m)+xQ(m)*qL;%column row of the matriz of points of the centroid distance function
    end    
    UB=unique(crM);%unique boxes that have at least one value
    countbsz(msz)=size(UB,2);
end
%box dimension
[pB,~] = polyfit(log(1./bsz),log(countbsz),1);
aux = pB(1);
