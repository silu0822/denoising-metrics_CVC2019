function [score,psnr] = denoisingMetrics(noisyImgName, denoisedImgName, cleanImgName)

% This is an implementation of our paper:
% =========================================================================
% Si Lu, "No-reference Image Denoising Quality Assessment". Computer Vision
% Conference, CVC 2019, Las Vegas.
% =========================================================================
%
% License: Only for academic use. All rights reserved to the authors of the
% paper (Si Lu). If you have any questions, comments or suggestions please 
% contact Si Lu at lusi@pdx.edu or daniel08220822@gmail.com. Please cite 
% our paper if you use this piece of code.
%
% How to run: 
% 1. go into `impl` folder and run `mex srMex.cpp`.
% Please note that this  is a personal implementation of an existing work: 
% "Chen, Jia, Chi-Keung Tang, and Jue Wang. Noise brush: interactive high 
% quality image-noise separation. ACM Transactions on Graphics (TOG). ACM, 
% 2009. PLEASE MAKE SURE TO CITE THE ABOVE PAPER AS WELL
%
% 2. go back to the root folder, run demo.m to see an running example  

%% load images and compute ground truth PSNR
dColorImage = imread(denoisedImgName);
nColorImge  = imread(noisyImgName);
cColorImge  = imread(cleanImgName);

imwrite(nColorImge(:,:,1),'nImg_R.png');
imwrite(nColorImge(:,:,2),'nImg_G.png');
imwrite(nColorImge(:,:,3),'nImg_B.png');

imwrite(dColorImage(:,:,1),'dImg_R.png');
imwrite(dColorImage(:,:,2),'dImg_G.png');
imwrite(dColorImage(:,:,3),'dImg_B.png');

psnr = PSNR(double(cColorImge), double(dColorImage));

%% R channel 
featuresR = computeFeatures('nImg_R.png', 'dImg_R.png'); 
load('RFR_PSNR_model_Rch_Ours_03.mat'); 
scoreR = regRF_predict(featuresR,model);

%% G channel
featuresG = computeFeatures('nImg_G.png', 'dImg_G.png'); 
load('RFR_PSNR_model_Gch_Ours_03.mat'); 
scoreG = regRF_predict(featuresG,model);

%% B channel
featuresB = computeFeatures('nImg_B.png', 'dImg_B.png'); 
load('RFR_PSNR_model_Bch_Ours_03.mat'); 
scoreB = regRF_predict(featuresB,model);

%% fusion R,G & B by mean
score = mean([scoreR, scoreG, scoreB]);

delete 'nImg_R.png' 'nImg_G.png' 'nImg_B.png' 'dImg_R.png' 'dImg_G.png' 'dImg_B.png'