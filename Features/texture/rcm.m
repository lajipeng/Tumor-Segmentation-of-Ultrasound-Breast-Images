function [x,feats] = rcm(I,BW)
% GLCM with ranklets
[~,J] = multiranklet(I,[2 4 8]);
[y,x] = find(BW);
xmn = min(x);
xmx = max(x);
ymn = min(y);
ymx = max(y);
J2 = J(ymn:ymx,xmn:xmx,:);
BW2 = BW(ymn:ymx,xmn:xmx);
% Distance 1
[x11,f11] = glcm_feats(J2(:,:,1),BW2,1,'R2');
[x21,f21] = glcm_feats(J2(:,:,2),BW2,1,'R4');
[x31,f31] = glcm_feats(J2(:,:,3),BW2,1,'R8');
% Distance 4
[x14,f14] = glcm_feats(J2(:,:,1),BW2,4,'R2');
[x24,f24] = glcm_feats(J2(:,:,2),BW2,4,'R4');
[x34,f34] = glcm_feats(J2(:,:,3),BW2,4,'R8');
% Features
x = [x11 x21 x31 x14 x24 x34];
feats = [f11,f21,f31,f14,f24,f34];
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [x,feats] = glcm_feats(I,B,D,str)
warning off;
I(~B) = NaN;
Offset = D*[0 1; -1 1; -1 0; -1 -1]; % Distance 1
CMs = graycomatrix(I,'Of',Offset,'NumLevels',64,'GrayLimits',[-1 1]);
[size_CM_1,size_CM_2,~] = size(CMs);
[i,j] = meshgrid(1:size_CM_1,1:size_CM_2);
idx1 = (i+j)-1;
idx2 = abs(i-j)+1;
entro = zeros(1,4);
corrm = zeros(1,4);
dissi = zeros(1,4);
homom = zeros(1,4);
for l = 1:4 % number of CMs
    CMaux = CMs(:,:,l);
    % Normalize GLCM
    CM_sum = sum(CMaux(:));
    Pij = CMaux./CM_sum;  % Normalize each CM
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
    % Dissimilarity
    dissi(l) = sum(sum(abs(i-j).*Pij));
    % Homogeneity Matlab
    homom(l) = sum(sum(Pij./(1+abs(i-j)))); 
    % Entropy
    entro(l) = -sum(sum(Pij.*log(Pij+eps)));
    %
    s_x = sum(sum(Pij.*((i-u_x).^2)))^0.5;
    s_y = sum(sum(Pij.*((j-u_y).^2)))^0.5;
    % Correlation Matlab;
    corm = sum(sum(Pij.*(i-u_x).*(j-u_y)));
    corrm(l) = corm/((s_x*s_y)+0.0001);    
end
xm = [mean(entro) mean(corrm)...
      mean(dissi) mean(homom)]; % Average
xd = [mean(entro-xm(1)) mean(corrm-xm(2)) ...
      mean(dissi-xm(3)) mean(homom-xm(4))]; % Absolute deviation
xr = [max(entro)-min(entro) max(corrm)-min(corrm) ...
      max(dissi)-min(dissi) max(homom)-min(homom)]; % Range
str = ['D' num2str(D) '_' str];
featsm = {['entro_mean_' str],['corr_mean_' str],...
          ['dissi_mean_' str],['homo_mean_' str]};
featsd = {['entro_mad_' str],['corr_mad_' str],...
          ['dissi_mad_' str],['homo_mad_' str]};
featsr = {['entro_rng_' str],['corr_rng_' str],...
          ['dissi_rng_' str],['homo_rng_' str]};      
feats = [featsm featsr featsd];
x = [xm xd xr];