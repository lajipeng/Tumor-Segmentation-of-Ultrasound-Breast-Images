%-- Author：Wang, Peng. Fudan University 
%-- Date：2019/11/30
%-- Fuction：Breast Tumor Segmentation of Ultra Images 
%-- Reference:
%   https://github.com/FubuFabian/Breast-Tumor-Segmentation-Matlab
%   https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1194626&tag=1
%% loading and initializating image
clear;clc;
load('C:\Github\医学成像技术\project2\Tumor-Segmentation-of-Ultrasound-Breast-Images\Input.mat')
load('C:\Github\医学成像技术\project2\Tumor-Segmentation-of-Ultrasound-Breast-Images\GT.mat')
%% segmentation, compute dice, show the comparasion between the input and ground truth 
nFrames = 95;
thresholdsegment_dice = [];
for i=1:nFrames   
    Smanual = logical(imfill(GT{i},'holes'));
    I = Input{i};
    S = thresholdsegment(I);
    thresholdsegment_dice = [thresholdsegment_dice,dice(S,Smanual)];
    imshowpair(Smanual,S,'falsecolor');
    h = findobj(gcf,'type','image');
    result = get(h,'CData');
    filename = sprintf('results/comparasion/THsegment%d.jpg',i);
    imwrite(result,filename);
end
thresholdsegment_dice_mean = mean(thresholdsegment_dice);
save('results/dice/thresholdsegment _dice.mat','thresholdsegment_dice','thresholdsegment_dice_mean');
%% compute time consuming 
t1=clock;
for i=1:nFrames   
    I = Input{i};
    S = thresholdsegment(I);
end
t2=clock;
time = etime(t2,t1)/nFrames;
disp(['etime程序总运行时间：',num2str(time)]);
save('results/time/time.mat','time');
