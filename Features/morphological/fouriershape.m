% FOURIERSHAPE Fourier-based shape features of a lesion.
%   [X_CSCF,FEATS_CSCF,X_CDF,FEATS_CDF] = FOURIERSHAPE(BW) computes two kind 
%   of functions to the shape signature of the lesion BW: centroid distance 
%   function (CDF) and centroid shape context function (CSCF). X is the numeric 
%   value of the features and FEATS is the name of the features in the same 
%   order as in X.
%
%   [X,FEATS] = FOURIERSHAPE(BW,FEATURES) computes computes a specific set 
%   of morphological features according with FEATURES. The set of valid 
%   strings includes (case insensitive):
%
%       'all'   - Entire feature set.
%       'cscf'  - Centroid shape context function.
%       'cdf'   - Centroid distance function.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS02.mat');   
%   [x_cscf,feats_cscf,x_cdf,feats_cdf] = fouriershape(BW);
%   % Same as [x_cscf,feats_cscf,x_cdf,feats_cdf] = fouriershape(BW,'all')
%
%   Example 2: Compute the centroid distance function
%   -------------------------------------------------
%   load('BUS02.mat');   
%   [x_cdf,feats_cdf] = fouriershape(BW,'cdf');
%
%   See also FOURIERFACTOR FRACTALCONTOUR POLYMODEL
%
%
%   Reference:
%   ---------
%   J. H. Lee, Y. K. Seong, C. H. Chang, et al., "Fourier-based shape feature 
%   extraction technique for computer-aided B-mode ultrasound diagnosis of 
%   breast tumor," in 34th Annual International Conference of the IEEE
%   EMBS, 2012, pp. 6551-6554.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   FOURIERSHAPE Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Arturo Rodriguez Cristerna, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function varargout = fouriershape(BW,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'CSCF';'CDF'};      
idxStats = RequestedInputs(officialStats,opts{:});
%*********************************************************************
Prps = regionprops(BW,'Centroid');
xc = Prps.Centroid(1);
yc = Prps.Centroid(2);
Baux = bwboundaries(BW);
C = Baux{1};
CDF= sqrt((C(:,1)-yc).^2+(C(:,2)-xc).^2);
fCDF=fft(CDF);
fCDF=padarray(fCDF,[25,0],0,'post');
fCDF(1)=0;
fNCDF=abs(fCDF./fCDF(2));
C(:,1) = C(:,1)-yc;
C(:,2) = C(:,2)-xc;
[t,r] = cart2pol(C(:,1),C(:,2));
td=discretize (t,12);%discretize theta
lrd=discretize (log(r)./log(2),5);%discretize log radio
bins=12;
lrd0=hist(td(lrd==0),bins);
lrd1=hist(td(lrd==1),bins);
lrd2=hist(td(lrd==2),bins);
lrd3=hist(td(lrd==3),bins);
lrd4=hist(td(lrd==4),bins);
CSCF=[lrd0 lrd1 lrd2 lrd3 lrd4];
fCSCF=fft(CSCF);
fCSCF(1)=0;
fNCSCF=abs(fCSCF./fCSCF(2));
%index 1 is gone because is zero
%index 2 is gone because is one
xfNCDF = fNCDF(3:25)';
xfNCSCF= fNCSCF(3:60);
featfNCDF={};
featfNCSCF={};
for n=3:25
    featfNCDF={featfNCDF{:} strcat('fNCDF_',num2str(n))};
end
for n=3:60
    featfNCSCF={featfNCSCF{:} strcat('fNCSCF_',num2str(n))};
end
if numel(idxStats) == 2
    varargout{1} = xfNCSCF;
    varargout{2} = featfNCSCF;
    varargout{3} = xfNCDF;
    varargout{4} = featfNCDF;
elseif idxStats == 1
    varargout{1} = xfNCSCF;
    varargout{2} = featfNCSCF;
elseif idxStats == 2
    varargout{1} = xfNCDF;
    varargout{2} = featfNCDF;
end
%*********************************************************************
function d = discretize (data,bins)
maxd=max(data);
mind=min(data);
d=round((bins-1)*((data-mind)/(maxd-mind)));