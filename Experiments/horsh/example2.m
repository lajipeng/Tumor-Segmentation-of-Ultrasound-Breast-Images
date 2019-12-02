%-- Author：Wang, Peng. Fudan University 
%-- Date：2019/11/30
%-- Fuction：Breast Tumor Segmentation of Ultra Images 
%-- Reference:
%   https://www.springerprofessional.de/en/busat-a-matlab-toolbox-for-breast-ultrasound-image-analysis/
%   https://www.tamps.cinvestav.mx/~wgomez/downloads.html
%% loading and initializating image
clear;clc;
load('C:\Github\医学成像技术\project2\Tumor-Segmentation-of-Ultrasound-Breast-Images\Input.mat')
load('C:\Github\医学成像技术\project2\Tumor-Segmentation-of-Ultrasound-Breast-Images\GT.mat')
%% segmentation, compute dice, show the comparasion between the input and ground truth 
nFrames = 95;
horsh_dice = [];
for i=1:nFrames   
    Smanual = logical(imfill(GT{i},'holes'));
    I = Input{i};
    S = horsch(I,Smanual);
    horsh_dice = [horsh_dice,dice(S,Smanual)];
    imshowpair(Smanual,S,'falsecolor');
    h = findobj(gcf,'type','image');
    result = get(h,'CData');
    filename =  sprintf('results/comparasion/horsch%d.jpg',i);
    imwrite(result,filename);
end
horsh_dice_mean = mean(horsh_dice);
save('results/dice/horsh_dice.mat','horsh_dice','horsh_dice_mean');
%% compute time consuming 
t1=clock;
for i=1:nFrames   
    Smanual = logical(imfill(GT{i},'holes'));
    I = Input{i};
    S = horsch(I,Smanual);
end
t2=clock;
time = etime(t2,t1)/nFrames;
disp(['etime程序总运行时间：',num2str(time)]);
save('results/time/time.mat','time');