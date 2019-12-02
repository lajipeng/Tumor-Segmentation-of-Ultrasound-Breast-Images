function [BW2,fb,params,fits] = pcnnsegment(I,S)
if nargin < 2
    S = [];
end
% Crop ROI
[I0,Gstr,rect] = getROI(I,S);
% Preprocessing
I1 = clahe(I0,[2 2],0.5);
I2 = isf(I1,15);
% Data structure
Data = datastr(I2,Gstr);
% Initialize population
[xi,fit,BWs] = initpop(Data);
% Differential evolution loop
tmax = Data.DE.tmax;
NP   = Data.DE.NP;
sze  = Data.Image.Size;
fits = zeros(NP,tmax);
for t = 1:tmax
    [F,CR] = deparams(t,tmax);
    % Mutation
    vi = mutation(xi,fit,F,Data);
    % Binomial crossover
    ui = crossover(xi,vi,CR,Data);
    % Evalute offspring
    BWaux = zeros(sze(1),sze(2),NP);
    fitaux = zeros(NP,1);
    parfor i = 1:NP
        [BWaux(:,:,i),fitaux(i)] = evalind(ui(i,:),Data);
    end
    % Selection
    [xi,fit,BWs] = selection(xi,ui,fit,fitaux,BWs,BWaux);
    % 
    disp(['Generation: ' num2str(t)]);
    fits(:,t) = fit;
end
% Get best solution
[fb,ib] = max(fit);
BW = logical(BWs(:,:,ib)-1);
params = xi(ib,:);
BW = postproc(BW,Data);
% Return to original image size
BW2 = false(size(I));
xw2min = rect(1);
xw2max = rect(2);
yh2min = rect(3);
yh2max = rect(4);
BW2(yh2min:yh2max,xw2min:xw2max) = BW;
%**********************************************************************
%************** BUS PROCESSING FUNCTIONS ******************************
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
G.strel0 = strel('disk',3,0);
G.strel1 = strel('disk',5,0);
G.strel2 = strel('disk',15,0);
%**********************************************************************
function I2 = gaussfun(I1,str)
[Mi,Ni] = size(I1);
w = str.width;  % ancho
h = str.height; % alto
xc = str.xc;
yc = str.yc;
K = [(w/2) 0;0 (h/2)].^2;
[x,y] = meshgrid(1:Ni,1:Mi);
G = exp(-(((x-xc).^2)/K(1,1) + ((y-yc).^2)/K(2,2)));
G = (G-min(G(:)))/(max(G(:))-min(G(:))); % Funcion Gaussiana
In = 1-(I1/max(I1(:)));
In = In.*G;
I2 = fix(64*(In-min(In(:)))/(max(In(:))-min(In(:))));
%**********************************************************************
function BW = postproc(BW,Data)
BW = imfill(BW,'holes');
BW = imopen(BW,Data.Image.Str.strel0);
A  = regionprops(BW,'Area'); 
As = cat(1,A.Area);
if ~isempty(As)
    BW = bwareaopen(BW,max(As)-1);
end
%**********************************************************************
%************** DIFFERENTIAL EVOLUTION FUNCTIONS **********************
function Data = datastr(I0,Gstr)
% Create synaptic matrices
Wr = 5;
W  = cell(1,Wr);
for r = 1:Wr
    W{r} = synapsis(r);
end
% Define limits for parameters
comun = [2 20;1 Wr]; % It, r  
r = [0.1 2.5;0.1 2.5;0.1 2;comun]; 
d = size(r,1);
% Objetive function
fun = @silindex;
% Constraint Gaussian function
I = gaussfun(I0,Gstr);
% Data
S  = double(I)+1;
Sn = S/max(S(:));
[M,N] = size(I);
% Population indices
NP = 20;
base = meshgrid(1:NP,1:NP);
base(logical(eye(NP))) = 0;
base = sort(base,2);
base(:,1) = [];
basej = meshgrid(1:d,1:NP);
% Distances
g = (1:65)/65;
[I1,I2] = meshgrid(g,g);
D = abs(I1-I2);
% Save data
Data.Image.Gray = S;
Data.Image.Norm = Sn;
Data.Image.Zeros = zeros(M,N);
Data.Image.Ones  = ones(M,N);
Data.Image.Size  = [M N];
Data.Image.Area  = M*N;
Data.Image.Mean  = mean2(Sn);
Data.Image.Str   = Gstr;
Data.Image.Dist  = D;
Data.Net.Synapsis = W;
Data.Net.Limits = r;
Data.Net.Params = d;
Data.Net.Name = @spcnn;
Data.DE.NP = NP;        % Population
Data.DE.tmax = 50;      % Max number of generations
Data.DE.Cost = fun;     % Objective function
Data.DE.Index = base;
Data.DE.j = basej;
%**********************************************************************
function [pop,fitness,BWs] = initpop(Data)
NP = Data.DE.NP;
sze = Data.Image.Size;
d = Data.Net.Params;
fitness = -1*inf(NP,1);
BWs = zeros(sze(1),sze(2),NP);
pop = zeros(NP,d);
pos = isinf(fitness);
spos = sum(pos);
limits = Data.Net.Limits;
l1 = limits(:,1)'; 
l2 = (limits(:,2)-limits(:,1))';
a = l1(ones(NP,1),:);
b = l2(ones(NP,1),:);
while spos > 0 % Repeat until all individuals have valid solutions
    % Create population randomly
    popaux = a + b.*rand(NP,d);
    popaux = trunk(popaux,limits([d-1 d],:));
    % Evaluate population
    j = find(pos);
    for i = 1:spos
        [BWs(:,:,j(i)),fitness(j(i))] = evalind(popaux(j(i),:),Data);
    end
    pos = isinf(fitness);
    spos = sum(pos);
    pop(~pos,:) = popaux(~pos,:);
end
%**********************************************************************
function [BW,fit] = evalind(Param,Data)
net  = Data.Net.Name;
cost = Data.DE.Cost;
str = param2struct(Param);
% Segment image with PCNN
BW = feval(net,Data,str);
% Evaluate fitness with CVI
fit = feval(cost,BW(:),Data);
%**********************************************************************
function [F,CR] = deparams(t,tmax)
F = 0.5*(1+rand);
CR = 0.5+(1-0.5)*(tmax-t)/tmax;
%**********************************************************************
function mpop = mutation(pop,fit,F,Data)
NP  = Data.DE.NP;
base = Data.DE.Index;
% Random indices
j = Shuffle(base,2);
% Apply mutation CURRENT-TO-BEST
[~,ib] = max(fit);
ib = ib*ones(NP,1);
mpop = pop + F*(pop(ib,:) - pop) + F*(pop(j(:,1),:) - pop(j(:,2),:));
%**********************************************************************
function vpop = crossover(pop,mpop,CR,Data)
NP = Data.DE.NP;
d = Data.Net.Params;
j = Data.DE.j;
rj = repmat(randi([1 d],NP,1),1,d);
vpop = pop;
idx = (rand(NP,d) < CR)|(j==rj);
vpop(idx) = mpop(idx);
vpop = bounce_back(pop,vpop,Data);
%**********************************************************************
function [Xi,fitness,BWs] = selection(Xi,Ui,fitness,fitaux,BWs,BWaux)
idx = fitaux >= fitness;
Xi(idx,:) = Ui(idx,:);
fitness(idx) = fitaux(idx);
BWs(:,:,idx) = BWaux(:,:,idx);
%**********************************************************************
function pop = bounce_back(popi,pop,Data)
limits = Data.Net.Limits; 
d = Data.Net.Params;
NP = Data.DE.NP;
LL = limits(:,1)';
UL = limits(:,2)';
LL = LL(ones(NP,1),:);
UL = UL(ones(NP,1),:);
idxU = pop>UL;
idxL = pop<LL;
rnd = rand(NP,d);
pop(idxU) = popi(idxU)+rnd(idxU).*(UL(idxU)-popi(idxU));
pop(idxL) = popi(idxL)+rnd(idxL).*(LL(idxL)-popi(idxL));
pop = trunk(pop,limits([d-1 d],:));
%**********************************************************************
function pop = trunk(pop,cmn)
N = size(pop,2);
pop(:,[N-1 N]) = fix(pop(:,[N-1 N]));
pop(:,N-1) = min(max(cmn(1,1),pop(:,N-1)),cmn(1,2)); % Trunk iterations iter = 2
pop(:,N)   = min(max(cmn(2,1),pop(:,N)),cmn(2,2));   % Trunk W radium   w = 2
%**********************************************************************
function f = silindex(BW,Data)
clrs = BW(:);
X  = Data.Image.Gray(:);
Nt = Data.Image.Area;
D  = Data.Image.Dist;
K = 2;	%number of clusters in this case is 2 (object / background)
% Validations about the number of clusters
if numel(unique(clrs)) ~= K
    f = -Inf;
    return;
else
    Nk = accumarray(clrs,ones(Nt,1),[K,1],@sum,0); % Total of elements in the cluster
    % validation about the size minimum of one cluster (100 pixels)
    if sum(Nk < 0.01*Nt) > 0
        f = -Inf;
        return;    
    end
end
S = zeros(1,K);
for i = 1:K
    j = setdiff(1:K,i);
    [ai,bi] = mDXX(X,D,clrs,i,j);
    S(i) = sum((bi-ai)./max([ai bi],[],2),1);
end
f = (1/Nt)*sum(S);
%********************************************************************
function [ai,bi] = mDXX(I,D,L,i,j)
% Intra cluster dispersion Ci
Ii = I(L==i);
Ui = unique(Ii);
XC2 = D(Ui,Ui);
Hi = accumarray(Ii,ones(numel(Ii),1),[65,1],@sum,0)';
Hi(Hi==0) = [];
Hi = Hi(ones(numel(Ui),1),:);
mi = sum(XC2.*Hi,2)/numel(Ii);
Si = nan(65,1);
Si(Ui) = mi;
ai = Si(Ii);
% Inter cluster dispersion Ci-Cj
Ij = I(L==j);
Uj = unique(Ij);
XC4 = D(Ui,Uj);
Hj = accumarray(Ij,ones(numel(Ij),1),[65,1],@sum,0)';
Hj(Hj==0) = [];
Hj = Hj(ones(numel(Ui),1),:);
mj = sum(XC4.*Hj,2)/numel(Ij);
Sj = nan(65,1);
Sj(Ui) = mj;
bi = Sj(Ii);
%**********************************************************************
function str = param2struct(Params)
str.B  = Params(1);
str.aT = Params(2);
str.vT = Params(3);
str.t  = Params(4);
str.r  = Params(5);
%**********************************************************************
%************************* SPCNN FUNCTIONS ****************************
function Y = spcnn(Data,Param)
% Initialize images
I = Data.Image.Norm;
Y = Data.Image.Zeros;
T = Data.Image.Ones;
% Matriz de sinapsis de interconexion con radio r
B = Param.B;
aT = Param.aT;
vT = Param.vT;
t = Param.t;
r = Param.r;
K = Data.Net.Synapsis{r};
% Proceso iterativo
for i = 1:t            
    F = I;
    L = imfilter(Y,K,'replicate');
    U = F.*(1+B*L);
    Y = im2double(U>T);             
    T = exp(-aT)*T + vT*Y;
end
Y = Y+1;
%**********************************************************************
function K = synapsis(r)
nk = 2*r + 1;
xi = r + 1;
yi = r + 1;
x = 1:nk;
y = (1:nk)';
K = meshgrid(abs(xi -x),abs(yi -y));
K = K + K' -1;
K = 1./2.^K;
K(r+1,r+1) = 0;