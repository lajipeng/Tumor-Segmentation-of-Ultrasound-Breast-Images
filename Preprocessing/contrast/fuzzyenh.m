% FUZZYENH Fuzzy enhancement.
%   J = FUZZYENH(I) performs the contrast enhancenment based on fuzzy logic 
%   to enhance the contrast of a breast ultrasound image I.
%
%   Example:
%   -------
%   I = imread('BUSreal.tif');
%   J = fuzzyenh(I);
%   figure; 
%   subplot 121; imshow(I); title('Original Image');
%   subplot 122; imshow(J); title('Contrasted Image');
%
%   See also CLAHE HISTEQU SACE SIGMOIDFILT  
%
%
%   Reference:
%   ----------
%   Y. Guo, H. Cheng, J. Huang, J. Tian, W. Zhao, L. Sun, Y. Su, 
%   "Breast ultrasound image enhancement using fuzzy logic," 
%   Ultrasound in Medicine and Biology, vol. 32 no. 2, pp. 237-247, 2006.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   FUZZYENH Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function g_out = fuzzyenh(go)
go = double(go);
%---------------------------------------------------------------------
% Paso 1. Normalizacion de la imagen de entrada
go_min = min(go(:));
go_max = max(go(:));
g_min = 0;
g_max = 255;
g = g_min + (((g_max-g_min)*(go-go_min))/(go_max-go_min)); % Ecuacion 1
g = fix(g);

%---------------------------------------------------------------------
% Paso 2. Fuzzificacion
% Paso 2.1. Calculo del histograma
n = 255;
hg = accumarray(fix(g(:))+1,ones(numel(g),1),[n+1 1],@sum,0);  % Histograma de la imagen

% Paso 2.2. Calculo del punto medio "y" de la funcion de membresia
% utilizando la maxima entropia con metodo de Kapur
pk = hg./sum(hg);   % Distribucion de probabilidad
Ps = cumsum(pk)+eps;                % Ecuacion 5
Hs = -cumsum(pk.*log(pk+eps));
Ps1 = 1-Ps; Ps1(Ps1<=0) = eps;
Ha = log(Ps)+(Hs./Ps);              % Ecuacion 3
Hb = log(Ps1)+((Hs(end)-Hs)./Ps1);  % Ecuacion 4
[~,ind] = max(Ha+Hb);
y = ind-1;

% Paso 2.3. Calculo de los valores extremos de la funcion de membresia "x"
% "z"
ind1 = [1,1:n];     % Desplazamiento hacia adelante
ind2 = [2:n+1,n+1]; % Desplazamiento hacia atras
peaks = (hg>hg(ind1))&(hg>hg(ind2));  % Calcula picos
x = find(peaks,1,'first')-1; % Valor inicial
z = find(peaks,1,'last')-1;  % Valor final

% Paso 2.4. Crea la funcion de membresia S-shape
% Ecuacion 2
S = zeros(size(hg));  % Inicializa funcion de membresia
gf = linspace(0,n,n+1);
% Rangos de la funcion
ind1 = gf<=x;
ind2 = (x<gf)&(gf<=y);
ind3 = (y<gf)&(gf<=z);
ind4 = z<gf;
% Funcion de membresia S(g;x,y,z)
S(ind1) = 0;
S(ind2) = ((gf(ind2)-x).^2)./((y-x)*(z-x));
S(ind3) = 1-((gf(ind3)-z).^2)./((z-y)*(z-x));
S(ind4) = 1;
% Fuzzifica la imagen
mu_ = S(g+1);

%---------------------------------------------------------------------
% Paso 3. Extrae informacion de bordes
% Paso 3.1. Calcula magnitud del gradiente con Sobel
hx = fspecial('sobel');
hy = hx';
Gx = conv2(mu_,hx,'same');
Gy = conv2(mu_,hy,'same');
delta_mu = sqrt(Gx.*Gx + Gy.*Gy);

% Paso 3.2. Normaliza informacion de bordes
delta_mu_min = min(delta_mu(:));
delta_mu_max = max(delta_mu(:));
e_mu = (delta_mu-delta_mu_min)./(delta_mu_max-delta_mu_min); % Ecuacion 10

%---------------------------------------------------------------------
% Paso 4. Extrae informacion de textura
% Paso 4.1. Crea mascaras de convolucion
L5 = [1,4,6,4,1];
E5 = [-1,-2,0,2,1];
S5 = [-1,0,2,0,-1];
M1 = L5'*E5;
M2 = L5'*S5;
M3 = E5'*L5;
M4 = S5'*L5;

% Paso 4.2. Calcula informacion de textura usando la imagen fuzzificada
T1 = abs(conv2(mu_,M1,'same'));
T2 = abs(conv2(mu_,M2,'same'));
T3 = abs(conv2(mu_,M3,'same'));
T4 = abs(conv2(mu_,M4,'same'));

% Paso 4.3. Textura conjunta
f_mu = (T1./max(T1(:))).*(T2./max(T2(:))).*(T3./max(T3(:))).*(T4./max(T4(:))); % Ecuacion 11

%---------------------------------------------------------------------
% Paso 5. Calcula el contraste
w = 5;   % Tamano de la ventana
% Paso 5.1. Calcula media local de tamano wxw
hw = ones(w);
numer = conv2(mu_.*f_mu.*e_mu,hw,'same');
denom = conv2(f_mu.*e_mu,hw,'same');
mu_w = numer./(denom+eps);    % Ecuacion 14

% Paso 5.2. Calcula contraste
C_mu = abs(mu_-mu_w)./(abs(mu_+mu_w)+eps); % Ecuacion 13

%---------------------------------------------------------------------
% Paso 6. Calcula valor de amplificacion k(i,j) y contraste modificado
% Paso 6.1. Calcula el contraste de la imagen original
s = 3;  % tamano de la ventana
hs = ones(s)/s^2;
g_s = conv2(g,hs,'same');  % Ecuacion 23 (media local)
C = abs(g-g_s)./(abs(g+g_s)+eps);    % Ecuacion 22

% Paso 6.2. Determina el valor medio del contraste
R = 0.01; % Weber ratio
G_R = C < R;
C_R = mean(C(G_R));  % Ecuacion 24
% Paso 6.3. Determina las constante de amplificacion minima y maxima
k_max = log(R)/log(C_R);    % Ecuacion 25 
k_min = log(R)/log(0.001);  % Ecuacion 26

% Paso 6.4. Calcula la entropia difusa local
E_mu = mu_.*f_mu.*e_mu;  % Ecuacion 20
phi_w = E_mu./(conv2(E_mu,hw,'same')+eps); % Ecuacion 19
E_n = -(1/log10(w*w))*conv2(phi_w.*log10(phi_w+eps),hw,'same'); % Ecuacion 18

% Paso 6.5. Calcula el exponente de amplificacion
E_n_min = min(E_n(:));  % Ecuacion 27a
E_n_max = max(E_n(:));  % Ecuacion 27b
k = k_min + (((E_n-E_n_min)*(k_max-k_min))./(E_n_max-E_n_min));  % Ecuacion 27

% Paso 6.6. Modifica el valor del contraste
C_mu_2 = C_mu.^k;  % Ecuacion 15

% Paso 6.7. Imagen fuzzificada modificada con Ecuacion 16
mu_2 = zeros(size(mu_));
ind1 = mu_>=mu_w;
ind2 = ~ind1;
mu_2(ind1) = mu_w(ind1).*(1+C_mu_2(ind1))./(1-C_mu_2(ind1));
mu_2(ind2) = mu_w(ind2).*(1-C_mu_2(ind2))./(1+C_mu_2(ind2));
mu_2(mu_2>1) = 1;
mu_2(mu_2<0) = 0;

%---------------------------------------------------------------------
% Paso 7. Defuzzificacion
S_inv = zeros(size(mu_2));  % Inicializa funcion de membresia
% Rangos de la funcion de membresia inversa
ind1 = (0<=mu_2)&(mu_2<((y-x)/(z-x)));
ind2 = (((y-x)/(z-x))<=mu_2)&(mu_2<=1);
% Funcion de membresia inversa Ecuacion 17
S_inv(ind1) = sqrt(mu_2(ind1)*(y-x)*(z-x));
S_inv(ind2) = (z-x-sqrt((1-mu_2(ind2))*(z-y)*(z-x)));
% Constante
gmin = min(S_inv(:));
gmax = max(S_inv(:));
K = gmin +((gmax-gmin)/(z-x));
% Defuzzifica image
g_out = uint8(K*S_inv);