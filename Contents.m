% The Breast Ultrasound Analysis Toolbox (BUSAT) contains 70 functions (m-files)
% to perfom image analysis including: preprocessing, segmentation, morphological
% and texture features, and binary classification (for two classes). This 
% toolbox was developed with MATLAB R2014a.
%
% Developed by Wilfrido Gomez Flores (wgomez@tamps.cinvestav.mx)
%           Arturo Rodriguez Cristerna (arodriguez@tamps.cinvestav.mx)
%           Wagner Coelho de Albuquerque Pereira (wagner.coelho@ufrj.br)
%
% Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%
% Please, cite the following paper where this toolbox was introduced:
%
%   Arturo Rodriguez-Cristerna, Wilfrido Gomez-Flores, and Wagner Coelho de 
%   Albuquerque-Pereira, ''BUSAT: A MATLAB Toolbox for Breast Ultrasound Image 
%   Analysis'', In Proceedings: 9th Mexican Congress on Pattern Recognition 
%   (MCPR 2017), LNCS 10267, pp. 268-277, 2017.
%
%
% IMPORTANT: First run "RUN_ME_FIRST.m" file to add this toolbox to search path.
%
% Symbol "*" indicates that the function requires a compiled C code.
% Despite the compiled functions are provided for MACI64, PCWIN64, and
% GLNXA64, sometimes it is necessary to re-compile from the source code
% using the mex function. Please see the notes in the "help" of the
% respective m-functions to get more details about the source C codes.
%
% 
%   IMAGE PREPROCESSING
%   -------------------
%       CONTRAST ENHANCEMENT:
%           clahe       - Contrast-limited adaptive histogram equalization.
%           fuzzyenh    - Fuzzy enhancement.
%           histequ     - Histogram equalization.
%           sace*       - Adaptive contrast enhancement based on sigmoidal mapping function.
%           sigmoidfilt - Sigmoid filter.
% 
%       DESPECKLING:
%           adf         - Anisotropic diffusion filtering.
%           adlg        - Anisotropic diffusion guided by Log-Gabor filters.
%           chmf        - Circular hybrid median filter.
%           isf         - Interference-based speckle filter.
%           isfad       - Interference-based speckle filter followed by anisotropic diffusion.
%
%       DOMAIN TRANSFORMATION:
%           multiloggabor - Multichannel decomposition by Log-Gabor filters.
%           multiranklet* - Multichannel decomposition by ranklet transform.
%           phasecong     - Phase congruency of a gray-scale image.
%           quantization  - Gray-level quantization.
%
%   LESION SEGMENTATION
%   -------------------
%       autosegment     - Automatic segmentation of breast lesions using texture analysis.
%       horsch          - Semiautomatic segmentation of breast lesions with the Horsch's method.
%       wtsdsegment     - Semiautomatic segmentation of breast lesions using watershed transform.
%
%   FEATURE EXTRACTION
%   ------------------
%       MORPHOLOGICAL:
%           convhulldiff    - Morphological features based on convex hull.
%           equivellipse    - Morphological features based on equivalent ellipse.
%           fourierfactor   - Fourier descriptors of lesion contour.
%           fouriershape    - Fourier-based shape features of a lesion.
%           fractalcontour  - Fractal dimension of the contour of a lesion.
%           geometric       - Geometric features of the lesion.
%           margclass       - Margin class features.
%           nrl             - Morphological features based on the normalized radial lenght curve.
%           nspd_li         - Number of protuberances and depressions and lobulation index.
%           polymodel       - Polygonal model of lesion contour.
%           spiculation     - Spiculation features of the lesion.
%
%       TEXTURE:
%           autocorr        - Gray-level autocorrelation coefficient.
%           autocov         - Gray-level autocovariance coefficients.
%           avmass          - Gray-level average-based texture features.
%           bdip_bvlc       - BDIP and BVLC moments.
%           bgc             - Binary gradient contours.
%           clxcurve        - Complexity curve.
%           fractaltexture  - Fractal features of gray-level image.
%           glcm            - Gray-level co-occurrence matrix features.
%           histfeatures    - First order statistics from gray-level histogram.
%           lawsenergy      - Laws' texture energy measures.
%           lbpv            - Local binary pattern variance from phase congruency.
%           nrg             - Normalized radial gradient.
%           pab             - Posterior acoustic behavior.
%
%       BI-RADS:
%           birads_feats    - BI-RADS features of US lexicon.
%           bound_feats     - Boundary features.
%           echo_feats      - Echo pattern features.
%           margin_feats    - Margin features.
%           orient_feats    - Orientation features.
%           shape_feats     - Shape features.
%
%   CLASSIFICATION
%   --------------
%       CLASSIFIERS:
%           bootstrap632AUC - AUC estimation using the .632+ bootstrap method.
%           bootstrap632ERR - Classification error estimation using the .632+ bootstrap method.
%           classifyLDA     - Classify data with Fisher's linear discriminant analysis.
%           classifyRBFN    - Classify data with radial basis function network.
%           classifySVM*    - Classify data with non-linear support vector machine.
%           classifyLSVM*   - Classify data with linear support vector machine.
%           classifyRF      - Classify data with random forest classifier.
%           classperf       - Classification performance.
%           trainLDA        - Train Fisher's linear discriminant analysis classifier.
%           trainRBFN       - Train radial basis function network.
%           trainSVM*       - Train non-linear support vector machine classifier.
%           trainLSVM*      - Train linear support vector machine classifier.
%           trainRF         - Train random forest classifier.
%
%       DATA NORMALIZATION:
%           minmaxnorm      - Minimum-Maximum normalization.
%           softmaxnorm     - Sigmoidal normalization.
%           statnorm        - Statistical normalization.
%       
%       FEATURE RANKING:
%           featselect      - Feature selection based on MRMR criterion.
%           mrmr_corr       - Minimal-redundancy-maximal-relevance (MRMR) criterion based on correlation.
%           mrmr_mi         - Minimal-redundancy-maximal-relevance (MRMR) criterion based on mutual information.
%
%   MISCELLANEOUS
%   -------------
%       AUC             - Area under the ROC curve.
%       ROCAUC          - Receiver operating characteristic curve.
%       cropROI         - Crop a region of interest.
%       delintumor      - Overlay tumor countour.
%
%
%