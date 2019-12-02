% AUTOSEGMENT Automatic segmentation of breast lesions using texture analysis.
%   [I,S] = AUTOSEGMENT(I) computes the automatic lesion segmentation (i.e., it 
%   does not require user-dependent parameters) using log-Gabor filtering
%   and texture analysis, where I is the breast ultrasound image. This method
%   was trained with 544 breast ultrasound images.
%
%   Example:
%   -------
%   load('BUS01.mat');
%   [I,S] = autosegment(I);
%   delintumor(I,Smanual); title('Radiologist outlining')
%   delintumor(I,S); title('Computerized segmentation')
%
%   See also MULTILOGGABOR WTSDSEGMENT HORSCH
%
%
%   Reference:
%   ---------
%   W. Gomez, B. A. Ruiz-Ortega, "New fully-automated method for segmentation
%   of breast lesions on ultrasound based on texture analysis," Ultrasound
%   in Medicine and Biology, vol. 42, no. 7, pp. 1637-1650, 2016.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (BraZil)
%   AUTOSEGMENT Version 1.0 (Matlab R2014a Unix)
%   December 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------
function [I0,BW] = autosegment(I)
if ~exist('train_data.mat','file')
    error('Training data is missing!');
else
    tr = load('train_data.mat');
end
% Compute Log-Gabor filtering and lattices
fprintf('Log-Gabor filtering and lattice division...\n');
[Data,Lattice] = buslattice(I);
% Compute texture features from lattices
fprintf('Computing texture features from lattice...\n');
XX = extractfeats(Data.Gabor,Lattice.WC,Lattice.WN,Lattice.CP,Lattice.MAP,Lattice.Window,Data.Levels);
% Lattice classification
fprintf('Lattice classification...\n');
Xtt = outliers(XX,[tr.lof;tr.uof]);
prb = gaussPDF(Lattice.Adjacency',tr.mu_,tr.isigma_);
Xtt = cat(2,Xtt,Lattice.Adjacency,prb);
Xtt = softmaxnorm(Xtt,[tr.m;tr.s])';
D = (Xtt-tr.b(:,ones(1,size(Xtt,2))))'*tr.W; % LDA scores
% Post-processing and segmentation
fprintf('Post-processing and contour refinement...\n');
I0 = Data.Original;
I1 = clahe(I0,[2 2],0.5);
I2 = isf(I1,11);
I3 = 1-I2/max(I2(:));
[Gx,Gy] = gradient(I2);
BW = postprocessing(D,I0,I3,Gx,Gy);
%*****************************************************************
function [Data,Lattice] = buslattice(I)
Q = 64; % Niveles de cuantizacion
W = 16;
% Procesamientos
Ig = multiloggabor(I,6,24,2,Q); % 6 escalas, 24 orientaciones, 2 octavas
[I,Ig,WC,WN,CP,MAP,ADJ] = imlattice(I,Ig,W); % Divide las imagenes en lattices
J = quantization(double(I),Q);        % Cuantiza la imagen original
% Estructuras para salvar
% Datos de imagenes
Data.Original = I;
Data.Gabor = Ig;
Data.Quanti = J;
Data.Levels = Q;
% Datos de lattice
Lattice.WC  = WC; % Center Window
Lattice.WN  = WN; % Neighborhod Windows
Lattice.CP  = CP; % Central point
Lattice.MAP = MAP; % Mapa de relacion de ventanas
Lattice.Adjacency = ADJ;
Lattice.Window = W;
%*****************************************************************
function [I,G,WC,WN,CP,MAP,ADJ] = imlattice(I,G,W)
% AJUSTA EL TAMANO DE LA IMAGEN
[Mi,Ni] = size(I);
Mb = 1:W:Mi;
Nb = 1:W:Ni;
I = I(1:Mb(end),1:Nb(end));
G = G(1:Mb(end),1:Nb(end),:);
% GENERA VENTANAS VECINAS
WN = ones(Mb(end),Nb(end));
WN(Mb,:) = 0;
WN(:,Nb) = 0;
WN = bwlabel(WN);
% GENERA VENTANAS CENTRALES
Mb2 = Mb+(W/2);
Mb2(end) = [];
Nb2 = Nb+(W/2);
Nb2(end) = [];
AUX = ones(Mb(end),Nb(end));
AUX(Mb2,:) = 0;
AUX(:,Nb2) = 0;
WC = zeros(Mb(end),Nb(end));
WC(Mb2(1):Mb2(end),Nb2(1):Nb2(end)) = AUX(Mb2(1):Mb2(end),Nb2(1):Nb2(end));
% GENERA PUNTOS DE CLASIFICACION
CP = bwmorph(WC,'shrink',Inf);
CP = bwlabel(CP);
WC = bwlabel(WC);
% GENERA EL MAPA DE VENTANAS
% [WC WN1 WN2 WN3 WN4]
Lc = unique(CP);
Lc(1) = []; % Elimina la etiqueta cero
Ls  = numel(Lc);
MAP = zeros(Ls,5);
for i = 1:Ls
    bw = WC == Lc(i);
    Lw = unique(WN(bw));
    Lw(1) = [];
    MAP(i,:) = [Lc(i) Lw'];
end
% INFORMACION DE ADYACENCIA
[y,x] = find(CP);
y = y/Mi; x = x/Ni;
ADJ = [x y];
%*****************************************************************
function XX = extractfeats(R,WC,WN,CP,MAP,W,Q)
% Offsets para la GLCM
o = [0 1;-1 1;-1 0;-1 -1];
Offset = [1*o;2*o;4*o;8*o];
% Tamanos de la imagen
[Mi,Ni,zi] = size(R);
Mb = numel(1:W:Mi);
Nb = numel(1:W:Ni);
% Arreglos vacios para almacenar descriptores
textWN = cell(Mb-1,Nb-1);
textWC = cell(Mb-2,Nb-2);
% Dimensionalidad
gc = 24*zi;
ac = 24*zi;
cr = zi;
ff = 8*zi;
cte = gc+ac+cr+ff;
% Extrae texturas
textWC = texturefeats(R,WC,textWC,Q,Offset); % Extrae texturas vecinos
textWN = texturefeats(R,WN,textWN,Q,Offset); % Extrae texturas centrales
MaxL = max(CP(:));
XX = zeros(MaxL,cte);
for i = 1:MaxL
    idx = MAP(i,:);
    XX(i,:) = nanmean([textWC{idx(1)};...
                       textWN{idx(2)};...
                       textWN{idx(3)};...
                       textWN{idx(4)};...
                       textWN{idx(5)}],1); % Promedia respuestas  
end
%*****************************************************************
function tWin = texturefeats(R,Win,tWin,Q,Offset)
zi = size(R,3);
MaxL = max(Win(:));
% Check if Parallel Computing Toolbox
res = isToolboxAvailable('Parallel Computing Toolbox','warning');
% If it is installed uses parfor for speed-up training
if res
    n = feature('numCores');        % Number of physical cores
    s = isempty(gcp('nocreate'));   % Check if pool is already open
    if s
        parpool('local',n); % Open parallel pool
    end
    parfor i = 1:MaxL
        tWin{i} = ith_texture(i,R,Win,zi,Q,Offset);
    end
    if s
        delete(gcp); % Close parallel pool
    end
else
    for i = 1:MaxL
        tWin{i} = ith_texture(i,R,Win,zi,Q,Offset);
    end
end
%*****************************************************************
function aux = ith_texture(i,R,Win,zi,Q,Offset)
[xc,yc] = find(Win==i);
xmin = min(xc); xmax = max(xc);
ymin = min(yc); ymax = max(yc);
RW  = R(xmin:xmax,ymin:ymax,:);
glcmW = zeros(24,zi);
acovW = zeros(24,zi);
acorrW= zeros(1,zi);
fracW = zeros(8,zi);
for j = 1:zi
    aux = RW(:,:,j); 
    glcmW(:,j)  = glcm(aux,Q,Offset);
    acovW(:,j)  = autocov(aux);
    acorrW(:,j) = autocorr(aux);
    fracW(:,j)  = fractfeat(aux);
end
aux = [acorrW;acovW;glcmW;fracW];
aux = aux(:)';
%*****************************************************************
function x = autocorr(US)
R = (double(US)+1).^2;
[NR,MR] = size(R);
% Autocorrelacion en profundidad
Cy = zeros(NR,MR);
for m = 1:MR
   for n = 1:NR
      Py = zeros(1,NR-n);
      for p = 1:NR-n
          Py(p) = R(n+p,m)*R(p,m);
      end
      Cy(n,m) = sum(Py);
   end
end
% Autocorrelacion lateral
Cy_bar = sum(Cy,2)+0.01;
% Autocorrelacion total
x = sum(Cy_bar./Cy_bar(1));

%*****************************************************************
function x = autocov(US)
I = double(US);
[M,N] = size(I);
f_bar = mean2(I); % Media imagen
% Computa matriz de autocovarianza de 5x5
sz = 5;
A = zeros(sz);
for dM = 1:sz
   for dN = 1:sz
        i = 1:(M - dM);
        j = 1:(N - dN);
        F1 = I(i,j)-f_bar;
        F2 = I(i+dM,j+dN) - f_bar;
        A(dM,dN) = (1/((M-dM)*(N-dN)))*sum(sum(F1.*F2));
   end
end
A = A+0.01;
% Salidas
N = (A./A(1,1))';
N = N(:)'; N(1) = [];
x = N';
%*****************************************************************
function x = fractfeat(US)
I = double(US);
I(1,1) = I(1,1)+1;
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
%*****************************************************************
function x = glcm(SI,Q,Offset)
CMs = graycomatrix(uint8(SI),'Of',Offset,'NumLevels',Q,'GrayLimits',[1 Q]);
% Promedia las caracteristicas a la misma distancia
[size_CM_1,size_CM_2,~] = size(CMs);
d = max(abs(Offset),[],2);
nd = unique(d);
size_CM_3 = numel(nd);
CM = cell(size_CM_3,4);%zeros(size_CM_1,size_CM_2,size_CM_3);
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
out.contr_avg = zeros(1,size_CM_3); % Contrast: matlab/[1,2]
out.corrm_avg = zeros(1,size_CM_3); % Correlation: matlab
out.entro_avg = zeros(1,size_CM_3); % Entropy: [2]
out.savgh_avg = zeros(1,size_CM_3); % Sum average [1]
out.senth_avg = zeros(1,size_CM_3); % Sum entropy [1]
out.homom_avg = zeros(1,size_CM_3); % Homogeneity: matlab
% Indices
[i,j] = meshgrid(1:size_CM_1,1:size_CM_2);
idx1 = (i+j)-1;
ii = (1:(2*size_CM_1-1))';
for k = 1:size_CM_3 % number CMs
    contr = zeros(1,4);
    entro = zeros(1,4);
    savgh = zeros(1,4);
    senth = zeros(1,4);
    corrm = zeros(1,4);
    homom = zeros(1,4);
    for l = 1:4
        CMaux = CM{k,l};
        % Normalize GLCM
        CM_sum = sum(CMaux(:));
        Pij = CMaux./CM_sum;  % Normalize each CM
        % 
        u_x = sum(sum(i.*Pij));
        u_y = sum(sum(j.*Pij));
        %
        s_x = sum(sum(Pij.*((i-u_x).^2)))^0.5;
        s_y = sum(sum(Pij.*((j-u_y).^2)))^0.5;
        %
        p_xplusy  = zeros((2*size_CM_1 - 1),1); %[1]
        for aux = 1:max(idx1(:))
           p_xplusy(aux) =  sum(Pij(idx1==aux));
        end
        % Contrast
        contr(l) = sum(sum((abs(i-j).^2).*Pij));
        % Entropy
        entro(l) = -sum(Pij(:).*log(Pij(:)+eps)); 
        % Sum average
        savgh(l) = sum((ii+1).*p_xplusy);
        % Sum entropy 
        senth(l) = -sum(p_xplusy.*log(p_xplusy+eps));
        % Correlation Matlab
        corm = sum(sum(Pij.*(i-u_x).*(j-u_y)));
        corrm(l) = corm/((s_x*s_y)+eps);
        % Homogeneity Matlab
        homom(l) = sum(sum(Pij./(1+abs(i-j))));  
    end
    % Contrast
    out.contr_avg(k) = mean(contr);
    % Entropy
    out.entro_avg(k) = mean(entro);  
    % Sum average
    out.savgh_avg(k) = mean(savgh);
    % Sum entropy 
    out.senth_avg(k) = mean(senth);
    % Correlation Matlab
    out.corrm_avg(k) = mean(corrm);
    % Homogeneity Matlab
    out.homom_avg(k) = mean(homom); 
end
% Salidas
nfeats = fieldnames(out);
numfea = size(nfeats,1);
x = zeros(1,numfea*size_CM_3);
a = 1;
for i = 1:numfea
    x(a:a+(size_CM_3-1)) = getfield(out,nfeats{i});
    a = a+size_CM_3;
end
x = x';
%*****************************************************************
function [Y,lof,uof] = outliers(X,stats)
if nargin == 1
    Q = prctile(X,[25 75],1);
    Q1 = Q(1,:);
    Q3 = Q(2,:);
    IQ = Q3-Q1;
    lof = Q1 - 12*IQ;
    uof = Q3 + 12*IQ;
elseif nargin == 2
    lof = stats(1,:);
    uof = stats(2,:);
end
N = size(X,1);
L = lof(ones(N,1),:);
U = uof(ones(N,1),:);
idL = X < L;
idU = X > U;
Y = X;
Y(idL) = L(idL);
Y(idU) = U(idU);
%*****************************************************************
function G = gaussPDF(x,m,iS)
N = size(x,2);
xm = bsxfun(@minus,x,m);
G = zeros(N,1);
for i = 1:N 
    G(i) = exp(-0.5*xm(:,i)'*iS*xm(:,i));
end
%*****************************************************************
function M = vect2mat(V,H,W,win)
Hr = (win+(H-2*win-1))/win;
Wr = (win+(W-2*win-1))/win;
M = reshape(V,Hr,Wr);
%*********************************************************************
function BW = postprocessing(D,I0,I3,Gx,Gy)
win = 16;
D = D/max(D);
[H,W] = size(I0);
RD = vect2mat(D,H,W,win);
ths = 0.2:0.05:0.8;
BWs  = zeros(H,W,numel(ths));
for i = 1:numel(ths)
    BW = RD >= ths(i);
    if sum(BW(:)) > 0
        BW = get_largest(BW);
        BW = imresize(BW,[H,W],'nearest');
        [G,ctrd] = gaussfun(BW);
        GI = fix(255*G.*I3);
        [Gxl,Gyl] = gradient(1-GI/max(GI(:)));
        B = ardfun(GI,Gxl,Gyl,H,W,ctrd);
        BWs(:,:,i)  = B+sum(BWs,3);
    else
        break;
    end
end
BWs = sum(BWs,3);
h = fspecial('gaussian',[3 3],1);
K = imfilter(BWs,h,'replicate');
K = K/max(K(:));
GI = fix(255*K.*I3);
out = zeros(1,255);
for i = 1:255
    BW = imfill(GI>=i,'holes');
    if sum(BW(:))>900
        [y,x] = find(BW);
        ctrd = [mean(x) mean(y)];
        out(i) = ard(bwperim(BW),Gx,Gy,ctrd);
    end
end
[~,th] = max(out);
BW = imfill(GI>=th,'holes');
%*********************************************************************
function B = get_largest(B)
LB = bwlabel(B,4);
ls = max(LB(:));
areas = zeros(1,ls);
for j = 1:ls
    areas(j) = sum(LB(:)==j);
end
if ~isempty(areas)
    [~,j] = max(areas);
    B = LB==j;
end
%*********************************************************************
function [G,ctrd] = gaussfun(BW)
[H,W] = size(BW);
[yb,xb] = find(BW);
[xi,yi] = meshgrid(1:W,1:H);
w = max(xb)-min(xb);
h = max(yb)-min(yb);
xc = fix(min(xb)+(w/2));
yc = fix(min(yb)+(h/2));
S = [(w/2)^2 0;0 (h/2)^2];
G = exp(-0.5*((xi-xc).^2/S(1,1) + (yi-yc).^2/S(2,2)))/(2*pi*sqrt(det(S)));
G = reshape(G,H,W);
G = G/max(G(:));
ctrd = [xc yc];
%*********************************************************************
function [B,mx] = ardfun(I,Gx,Gy,M,N,ctrd)
ards = zeros(1,255);
Bs = false(M,N,255);
for i = 1:255
    B  = imfill(I >= i,'holes');
    B = get_largest(B);
    if sum(B(:)) > 900;
        C = bwperim(B);
        ards(i) = ard(C,Gx,Gy,ctrd);
        Bs(:,:,i) = B;
    else
        continue;
    end
end
[mx,th] = max(ards);
BW = Bs(:,:,th);
B = double(BW);
%*********************************************************************
function out = ard(cont,Gx,Gy,ctrd)
xc = ctrd(1);
yc = ctrd(2);
[cy,cx] = find(cont); 
inds = sub2ind(size(cont),cy,cx);
N = numel(inds);
r = [cx-xc cy-yc];
norma = sqrt(r(:,1).^2+r(:,2).^2);
rn = r./norma(:,ones(1,2));
out = sum(diag([Gx(inds) Gy(inds)]'*rn))/N;
%*********************************************************************