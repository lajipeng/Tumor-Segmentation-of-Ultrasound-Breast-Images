% GLCM Gray-level co-occurrence matrix features.
%   [X,FEATS,CM] = GLCM(I,BW,Q,DIST,OPTION1,OPTION2) computes 20-by-distances  
%   glcm features from the gray-level image I masked by the binary image BW,  
%   where Q is a positive integer to indicate the quantization level, DIST   
%   is a vector specifying the distances between the pixel-of-interest and 
%   its neighbor. If OPTION=0 computes the glcm from the whole image. If 
%   OPTION1=1 computes the glcm only from the internal region of the lesion.
%   OPTION2 is a string that indicates the statitistic used to merge 
%   the features at the same distance overall the orientations: 'mean', 'mad',
%   and 'range'. X is a numeric vector with the feature values, FEATS is a 
%   cell vector with the name of the features in the same order as in X, and 
%   CM are the glcm matrices.
%
%   [X,FEATS] = glcm(I,BW,Q,DIST,OPTION1,OPTION2,FEATURES) computes a specific 
%   set of texture features according with FEATURES. The set of valid strings  
%   includes (case insensitive):
%
%       'autoc' - Autocorrelation
%       'contr' - Contrast
%       'corrm' - Correlation
%       'corrh' - Haralick's correlation
%       'cprom' - Cluster Prominence
%       'cshad' - Cluster Shade
%       'dissi' - Dissimilarity
%       'energ' - Energy
%       'entro' - Entropy
%       'homom' - Homogeneity
%       'maxpr' - Maximum probability
%       'sosvh' - Sum of squares (Variance)
%       'savgh' - Sum average
%       'svarh' - Sum variance
%       'senth' - Sum entropy
%       'dvarh' - Difference variance
%       'denth' - Difference entropy
%       'inf1h' - Information measure of correlation 1
%       'inf2h' - Informaiton measure of correlation 2
%       'indnc' - Inverse difference normalized
%       'idmnc' - Inverse difference moment normalized
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS02.mat');
%   D = [1 2 4 8]; % Distances in pixels
%   [x,feats] = glcm(I,BW,64,D,1,'mean');
%
%   Example 2: Compute contrast, correlation and entropy
%   ----------------------------------------------------
%   load('BUS02.mat');
%   D = [1 2 4 8]; % Distances in pixels
%   [x,feats] = glcm(I,BW,64,D,1,'mean','contr','corrm','entro');
%
%   See also AUTOCORR AUTOCOV BDIP_BVLC LAWSENERGY
%
%
%   References:
%   ----------
%   Haralick, R.M., K. Shanmugan, I. Dinstein, "Textural features for 
%   image classification", IEEE Transactions on Systems, Man, and 
%   Cybernetics, vol. SMC-3, pp. 610-621, 1973.
% 
%   Haralick, R.M., L.G. Shapiro. Computer and Robot Vision: Vol. 1, 
%   Addison-Wesley, pp. 459, 1992.
%
%   W. Gomez, W. Pereira, A.F.C. Infantosi, "Analysis of co-occurrence texture 
%   statistics as a function of gray-level quantization for classifying breast 
%   ultrasound," IEEE Trans. Med. Imaging, vol. 31, no. 10, pp. 1889-1899, 2012.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   GLCM Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x2,feats2,CM] = glcm(US,BW,Q,D,opt1,opt2,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'autoc';'contr';'corrm';'corrh';'cprom';'cshad';...
                 'dissi';'energ';'entro';'homom';'maxpr';'sosvh';...
                 'savgh';'svarh';'senth';'dvarh';'denth';'inf1h';...
                 'inf2h';'indnc';'idmnc'}; 
idxStats = RequestedInputs(officialStats,opts{:});
idxStats = sort(idxStats,'ascend');
%*********************************************************************
warning off;
opt1 = logical(opt1);
opt2 = lower(opt2);
if strcmpi(opt2,'mean') || strcmpi(opt2,'mad') || strcmpi(opt2,'range')
    fun = str2func(opt2);
else
    error('Only the statistics mean, mean absolute deviation (mad), and range are accepted.');
end
[y,x] = find(BW);
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
J = quantization(US(ymin:ymax,xmin:xmax),Q); % Intensity quantization
if opt1
    B = BW(ymin:ymax,xmin:xmax);
    J(~B) = NaN;
end
theta = [0 1; -1 1; -1 0; -1 -1]; % Four orientations 0, 45, 90 and 135
D = unique(D);
Offset = [];
for i = 1:numel(D)
    Offset = cat(1,Offset,D(i)*theta);
end
CMs = graycomatrix(J,'Of',Offset,'NumLevels',Q,'GrayLimits',[1 Q]);
% Promedia las caracteristicas a la misma distancia
[size_CM_1,size_CM_2,~] = size(CMs);
d = max(abs(Offset),[],2);
nd = unique(d);
size_CM_3 = numel(nd);
CM = cell(size_CM_3,4);
pixd = zeros(1,size_CM_3);
for i = 1:size_CM_3
    idx = d == nd(i);
    pixd(i) = nd(i);
    CMaux = CMs(:,:,idx);
    for j = 1:4
        CM(i,j) = {CMaux(:,:,j)};
    end
end
% checked
out.autoc = zeros(1,size_CM_3); % Autocorrelation: [2] 
out.contr = zeros(1,size_CM_3); % Contrast: matlab/[1,2]
out.corrm = zeros(1,size_CM_3); % Correlation: matlab
out.corrh = zeros(1,size_CM_3); % Correlation Haralick
out.cprom = zeros(1,size_CM_3); % Cluster Prominence: [2]
out.cshad = zeros(1,size_CM_3); % Cluster Shade: [2]
out.dissi = zeros(1,size_CM_3); % Dissimilarity: [2]
out.energ = zeros(1,size_CM_3); % Energy: matlab / [1,2]
out.entro = zeros(1,size_CM_3); % Entropy: [2]
out.homom = zeros(1,size_CM_3); % Homogeneity: matlab
out.maxpr = zeros(1,size_CM_3); % Maximum probability: [2]
out.sosvh = zeros(1,size_CM_3); % Sum of squares: Variance [1]
out.savgh = zeros(1,size_CM_3); % Sum average [1]
out.svarh = zeros(1,size_CM_3); % Sum variance [1]
out.senth = zeros(1,size_CM_3); % Sum entropy [1]
out.dvarh = zeros(1,size_CM_3); % Difference variance [4]
out.denth = zeros(1,size_CM_3); % Difference entropy [1]
out.inf1h = zeros(1,size_CM_3); % Information measure of correlation1 [1]
out.inf2h = zeros(1,size_CM_3); % Informaiton measure of correlation2 [1]
out.indnc = zeros(1,size_CM_3); % Inverse difference normalized (INN) [3]
out.idmnc = zeros(1,size_CM_3); % Inverse difference moment normalized [3]

% Indices
[i,j] = meshgrid(1:size_CM_1,1:size_CM_2);
idx1 = (i+j)-1;
idx2 = abs(i-j)+1;
ii = (1:(2*size_CM_1-1))';
jj = (0:size_CM_1-1)';
for k = 1:size_CM_3 % number of CMs
    autoc = zeros(1,4); 
    contr = zeros(1,4);
    corrm = zeros(1,4);
    corrh = zeros(1,4);
    cprom = zeros(1,4);
    cshad = zeros(1,4);
    dissi = zeros(1,4);
    energ = zeros(1,4);
    entro = zeros(1,4);
    homom = zeros(1,4);
    maxpr = zeros(1,4);
    sosvh = zeros(1,4);
    savgh = zeros(1,4);
    svarh = zeros(1,4);
    senth = zeros(1,4);
    dvarh = zeros(1,4);
    denth = zeros(1,4);
    inf1h = zeros(1,4);
    inf2h = zeros(1,4);
    indnc = zeros(1,4);
    idmnc = zeros(1,4);
    for l = 1:4
        CMaux = CM{k,l};
        % Normalize GLCM
        CM_sum = sum(CMaux(:));
        Pij = CMaux./CM_sum;  % Normalize each CM
        CM_mean = mean(Pij(:));    % compute mean after norm
        %
        p_x = squeeze(sum(Pij,2));
        p_y = squeeze(sum(Pij,1))';
        % 
        u_x = sum(sum(i.*Pij));
        u_y = sum(sum(j.*Pij));
        %
        p_xplusy  = zeros((2*size_CM_1 - 1),1); %[1]
        p_xminusy = zeros((size_CM_1),1);       %[1]
        for aux = 1:max(idx1(:))
           p_xplusy(aux) =  sum(Pij(idx1==aux));
        end
        for aux = 1:max(idx2(:))
           p_xminusy(aux) = sum(Pij(idx2==aux));
        end
        % Contrast
        contr(l) = sum(sum((abs(i-j).^2).*Pij));
        % Dissimilarity
        dissi(l) = sum(sum(abs(i-j).*Pij));
        % Energy
    	energ(l) = sum(sum(Pij.^2));
        % Entropy
        entro(l) = -sum(sum(Pij.*log(Pij+eps)));
        % Homogeneity Matlab
        homom(l) = sum(sum(Pij./(1+abs(i-j)))); 
        % Sum of squares: Variance
        sosvh(l) = sum(sum(Pij.*((j-CM_mean).^2)));
        % Inverse difference normalized
    	indnc(l) = sum(sum(Pij./(1+((abs(i-j).^2)./(size_CM_1.^2)))));
        % Inverse difference moment normalized
    	idmnc(l) = sum(sum(Pij./(1+(((i-j).^2)./(size_CM_1.^2)))));
        % Maximum probability
    	maxpr(l) = max(Pij(:));
        % Sum average
        savgh(l) = sum((ii+1).*p_xplusy);
        % Sum entropy 
        senth(l) = -sum(p_xplusy.*log(p_xplusy+eps));
        % Sum variance
    	svarh(l) = sum((((ii+1) - senth(l)).^2).*p_xplusy);
        % Difference entropy
    	denth(l) = -sum(p_xminusy.*log(p_xminusy+eps));
        % Difference variance
    	dvarh(l) = sum((jj.^2).*p_xminusy);
        % Computes correlation
        hxy1 = -sum(sum(Pij.*log(p_x*p_y' + eps)));
        hxy2 = -sum(sum((p_x*p_y').*log(p_x*p_y' + eps)));
        hx   = -sum(p_x.*log(p_x+eps));
        hy   = -sum(p_y.*log(p_y+eps));
        hxy  = entro(l);
        % Information measure of correlation 1 
    	inf1h(l) = (hxy-hxy1)/(max([hx,hy]));
        % Information measure of correlation 2
    	inf2h(l) = (1-exp(-2*(hxy2-hxy)))^0.5;
        % Cluster Prominence
    	cprom(l) = sum(sum(Pij.*((i+j-u_x-u_y).^4)));
        % Cluster Shade
    	cshad(l) = sum(sum(Pij.*((i+j-u_x-u_y).^3)));
        %
        s_x = sum(sum(Pij.*((i-u_x).^2)))^0.5;
        s_y = sum(sum(Pij.*((j-u_y).^2)))^0.5;
        autoc(l) = sum(sum(Pij.*(i.*j)));
        % Correlation Matlab;
        corm = sum(sum(Pij.*(i-u_x).*(j-u_y)));
        corrm(l) = corm/((s_x*s_y)+0.0001);
        % Correlation Haralick
        corrh(l) = ((sum(sum(Pij.*(i.*j))))-u_x*u_y)/(s_x*s_y);
    end
    out.autoc(k) = feval(fun,autoc);
    out.contr(k) = feval(fun,contr);
    out.corrm(k) = feval(fun,corrm);
    out.corrh(k) = feval(fun,corrh);
    out.cprom(k) = feval(fun,cprom);
    out.cshad(k) = feval(fun,cshad);
    out.dissi(k) = feval(fun,dissi);
    out.energ(k) = feval(fun,energ);
    out.entro(k) = feval(fun,entro);
    out.homom(k) = feval(fun,homom);
    out.maxpr(k) = feval(fun,maxpr);
    out.sosvh(k) = feval(fun,sosvh);
    out.savgh(k) = feval(fun,savgh);
    out.svarh(k) = feval(fun,svarh);
    out.senth(k) = feval(fun,senth);
    out.dvarh(k) = feval(fun,dvarh);
    out.denth(k) = feval(fun,denth);
    out.inf1h(k) = feval(fun,inf1h);
    out.inf2h(k) = feval(fun,inf2h);
    out.indnc(k) = feval(fun,indnc);
    out.idmnc(k) = feval(fun,idmnc);
end
% Salidas
nfeats = fieldnames(out);
numfea = size(nfeats,1);
x = zeros(1,numfea*size_CM_3);
feats = cell(1,numfea*size_CM_3);
a = 1;
for i = 1:numfea
    fea = cell(1,size_CM_3);
    for j = 1:size_CM_3
        fea(j) = {['glcm_' nfeats{i} '_D' num2str(pixd(j)) '_' opt2]};
    end
    x(a:a+(size_CM_3-1)) = getfield(out,nfeats{i});
    feats(a:a+(size_CM_3-1)) = fea;
    a = a+size_CM_3;
end

feats2 = [];
x2 = [];
for i = 1:numel(idxStats)
    idx = not(cellfun(@isempty,strfind(feats,officialStats{idxStats(i)})'));
    feats2 = cat(2,feats2,feats(idx));
    x2 = cat(2,x2,x(idx));
end