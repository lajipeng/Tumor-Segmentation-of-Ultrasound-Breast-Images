%-- Author：Wang, Peng. Fudan University 
%-- Date：2019/11/30
%-- Fuction：Breast Tumor Segmentation of Ultra Images 
%-- Reference:
%   https://www.springerprofessional.de/en/busat-a-matlab-toolbox-for-breast-ultrasound-image-analysis/
%   https://www.tamps.cinvestav.mx/~wgomez/downloads.html
%% loading and initializating image
clear;clc;
load('I:\学习\医学成像技术\project2\US Toolbox Ver. 2.0\Input.mat')
load('I:\学习\医学成像技术\project2\US Toolbox Ver. 2.0\GT.mat')
%% segmentation, compute dice, show the comparasion between the input and ground truth 
nFrames = 95;
autosegment_dice = [];
for i=53:nFrames   
    Smanual = logical(imfill(GT{i},'holes'));
    Smanual = Smanual(7:215,7:295);
    I = Input{i};
    [I,S] = autosegment(I);
    autosegment_dice = [autosegment_dice,dice(S,Smanual)];
    imshowpair(Smanual,S,'falsecolor');
    h = findobj(gcf,'type','image');
    result = get(h,'CData');
    filename = sprintf('results/comparasion/autosegment%d.jpg',i);
    imwrite(result,filename);
end
autosegment_dice_mean = mean(autosegment_dice);
save('results/dice/autosegment _dice.mat','autosegment_dice','autosegment_dice_mean');
%% compute time consuming 
t1=clock;
for i=1:nFrames   
    I = Input{i};
    S = autosegment(I);
end
t2=clock;
time = etime(t2,t1)/nFrames;
disp(['etime程序总运行时间：',num2str(time)]);
save('results/time/time.mat','time');
