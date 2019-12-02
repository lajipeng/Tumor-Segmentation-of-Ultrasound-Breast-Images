function ModelBWCompare(subset)
%-- Author：Wang, Peng. Fudan University 
%-- Date：2019/11/30
%-- Fuction：show the piture of segmentation result from diffirent model
%% loading and initializating image
load('C:/Github/医学成像技术/project2/Tumor-Segmentation-of-Ultrasound-Breast-Images/Input.mat')
%% imshow 
original_image = Input{subset};
path_autosegment = sprintf('C:/Github/医学成像技术/project2/US Toolbox Ver. 2.0/Experiments/autosegment/results/comparasion/autosegment%d.jpg',subset);
autosegment_image = imread(path_autosegment);
path_horsh = sprintf('C:/Github/医学成像技术/project2/US Toolbox Ver. 2.0/Experiments/horsh/results/comparasion/horsch%d.jpg',subset);
horsh_image = imread(path_horsh);
path_wtsdsegment = sprintf('C:/Github/医学成像技术/project2/US Toolbox Ver. 2.0/Experiments/wtsdsegment/results/comparasion/wtsdsegment%d.jpg',subset);
wtsdsegment_image = imread(path_wtsdsegment);
path_threshold = sprintf('C:/Github/医学成像技术/project2/US Toolbox Ver. 2.0/Experiments/threshold/results/comparasion/THsegment%d.jpg',subset);
threshold_image = imread(path_threshold);
figure;
subplot 151; imshow(original_image,'InitialMagnification', 'fit'); xlabel('(a)');
subplot 152; imshow(threshold_image,'InitialMagnification', 'fit'); xlabel('(b)');
subplot 153; imshow(horsh_image,'InitialMagnification', 'fit'); xlabel('(c)');
subplot 154; imshow(autosegment_image,'InitialMagnification', 'fit'); xlabel('(d)');
subplot 155; imshow(wtsdsegment_image,'InitialMagnification', 'fit'); xlabel('(e)');
end

