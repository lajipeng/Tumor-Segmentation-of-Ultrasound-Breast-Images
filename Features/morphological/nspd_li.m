% NSPD_LI Number of protuberances and depressions and lobulation index.
%   [X,FEATS] = NSPD_LI(BW) computes two morphological features from the
%   binary blob of a breast lesion BW: number of substantial protuberances 
%   and  depressions and lobulation index. X is a numeric vector with the 
%   feature values and FEATS is a cell vector with the name of the features 
%   in the same order as in X.
%
%   [X,FEATS] = NSPD_LI(BW,FEATURES) computes computes a specific set 
%   of morphological features according with FEATURES. The set of valid 
%   strings includes (case insensitive):
%
%       'all'   - Entire feature set.
%       'nspd'  - Number of substantial protuberances and depressions.
%       'li'    - Lobulation index.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS02.mat');   
%   [x,feats] = nspd_li(BW);
%   % Same as [x,feats] = nspd_li(BW,'all');
%
%   Example 2: Compute lobulation index
%   ----------------------------------
%   load('BUS02.mat');   
%   [x,feats] = nspd_li(BW,'li');
%
%   See also MARGCLASS SPICULATION
%
%
%   Reference:
%   ---------
%   C.-M. Chen, Y.-H. Chou, K.-C. Han, et al., "Breast lesions on sonograms: 
%   computer-aided diagnosis with nearly setting-independent features and 
%   artificial neural networks," Radiology, vol. 226, no. 2, pp. 504-514,2003.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   NSPD_LI Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = nspd_li(BW,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'nspd';'li'};
idxStats = RequestedInputs(officialStats,opts{:});
%*********************************************************************
% Propiedades del tumor
BW   = double(BW);
props = regionprops(BW,'ConvexHull','Area','Centroid','Perimeter');
xc = props.Centroid(1);
yc = props.Centroid(2);
% Extrae y parametriza contorno del tumor 
junk = bwboundaries(BW);
cBW  = junk{1};
yBW  = cBW(:,1); xBW = cBW(:,2);
% Empareja coordenadas del poligono convexo y del tumor
xCH = props.ConvexHull(:,1);
yCH = props.ConvexHull(:,2);
D1  = dist([xBW yBW],[xCH yCH]');
[~,ind] = min(D1);
xCH = xBW(ind);
yCH = yBW(ind);
% Crea poligino convexo
se = [0 1 0;1 1 1;0 1 0];
CH = roipoly(BW,xCH,yCH);
% Depresiones
dif = xor(BW,CH);
dif2 = imreconstruct(imerode(dif,se),dif);
if sum(dif2(:)) > 0
    % Remueve depresiones de area pequena
    propCH = regionprops(dif2,'Area');
    aCH = cat(1, propCH.Area);
    por = (aCH*100)/max(aCH);
    idx = por < round(0.5*median(por));
    aCH(idx) = [];
    dif2 = bwareaopen(dif2,min(aCH));
    % Etiquetas para recorrido
    L = bwlabel(dif2);
    dBW = imdilate(BW,se);
    junk = bwboundaries(dBW);
    cdBW  = junk{1};
    ydBW  = cdBW(:,1); xdBW = cdBW(:,2);
    ind = sub2ind(size(dBW),ydBW,xdBW);
    ets = L(ind); ets(ets==0)=[];
    ets = myunique(ets);
    % Identifica lineas rectas
    cDF  = bwperim(dif2);
    cCH  = bwperim(CH);
    rect = cCH&cDF;
    curv = ~rect&cDF;
    Lrect = zeros(size(L));
    Lcurv = zeros(size(L));
    Lrect(rect) = L(rect);
    Lcurv(curv) = L(curv);
    xP = zeros(numel(ets),1);
    yP = zeros(numel(ets),1);
    for i = 1:numel(ets)
        [yR,xR] = find(Lrect==ets(i));
        [yC,xC] = find(Lcurv==ets(i));
        nR = sqrt(sum([xR-xc yR-yc].^2,2));
        uR = [xR-xc yR-yc]./(repmat(nR,1,2)+eps);
        nC = sqrt(sum([xC-xc yC-yc].^2,2));
        uC = [xC-xc yC-yc]./(repmat(nC,1,2)+eps);
        D1 = dist(uC,uR');
        [~,ind1] = min(D1,[],2);
        xR2 = xR(ind1); yR2 = yR(ind1);
        D2 = sqrt((xR2-xC).^2+(yR2-yC).^2);
        [~,ind2] = max(D2);
        xP(i) = xC(ind2); yP(i) = yC(ind2);  
    end
    ets2 = [1:numel(ets) 1]';
    etseq = zeros(2*numel(ets2),1);
    etseq(1:2:end) = ets2;
    etseq(2:2:end) = ets2;
    etseq([1 end]) = [];
    etseq = reshape(etseq,2,numel(etseq)/2)';
    xPc = xP-round(xc); yPc = yP-round(yc);
    BW2 = BW;
    for i = 1:size(etseq,1)
        A = [xPc(etseq(i,1)) yPc(etseq(i,1))];
        B = [xPc(etseq(i,2)) yPc(etseq(i,2))];
        ang = real(acos(dot(A,B)/(norm(A)*norm(B)+eps)))*(180/pi);
        if ang < 120
            [x,y] = intline(xP(etseq(i,1)),xP(etseq(i,2)),yP(etseq(i,1)),yP(etseq(i,2)));
            indxy = sub2ind(size(BW),y,x);
            BW3 = ones(size(BW2));
            BW3(indxy) = 0; BW3 = imerode(BW3,ones(2));
            BW2 = BW2.*BW3;
        end
    end
    LBW2 = bwlabel(BW2);
    val = LBW2(round(yc),round(xc));
    LBW2(LBW2==val) = 0;
    propProt = regionprops(logical(LBW2),'Area');
    aCH = cat(1, propProt.Area);
    por = (aCH*100)/max(aCH);
    aCH(por <= 5) = [];
    nprot = numel(aCH); % Numero de protuberancias
    LBW3 = bwlabel(dif2);
    dEt = unique(LBW3);
    dEt(dEt==0) = [];
    ndepr = numel(dEt); % Numero de depresiones
    npd = nprot+ndepr;
    if nprot > 1
        LI = (max(aCH)-min(aCH))/(mean(aCH)); %indice de lobulacion
    else
        LI = 0;
    end
%      figure; 
%      imshow(dif2+LBW2,[]); hold on;
%      plot(xBW,yBW,'y');
%      plot(xc,yc,'y*')
else
   npd = 0;
   LI  = 0;
end
aux_x = [npd,LI];
aux_f = {'NSPD','LobIndex'};
x = aux_x(idxStats);
feats = aux_f(idxStats);
%*********************************************************************
function [a,ndx] = myunique(a)
[r,c]=size(a);
[s,si]=sort(a);
ds = false(r,c);
esit=diff(s)==0;
ds(1:end-1) = esit;
si(ds)=[];
ndx=sort(si);
a=a(ndx);
%*********************************************************************
function [x,y] = intline(x1, x2, y1, y2)
dx = abs(x2 - x1);
dy = abs(y2 - y1);
% Check for degenerate case.
if ((dx == 0) && (dy == 0))
  x = x1;
  y = y1;
  return;
end
flip = 0;
if (dx >= dy)
  if (x1 > x2)
    % Always "draw" from left to right.
    t = x1; x1 = x2; x2 = t;
    t = y1; y1 = y2; y2 = t;
    flip = 1;
  end
  m = (y2 - y1)/(x2 - x1);
  x = (x1:x2).';
  y = round(y1 + m*(x - x1));
else
  if (y1 > y2)
    % Always "draw" from bottom to top.
    t = x1; x1 = x2; x2 = t;
    t = y1; y1 = y2; y2 = t;
    flip = 1;
  end
  m = (x2 - x1)/(y2 - y1);
  y = (y1:y2).';
  x = round(x1 + m*(y - y1));
end
if (flip)
  x = flipud(x);
  y = flipud(y);
end