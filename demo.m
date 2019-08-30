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
% 1. go into `core` folder and run `mex srMex.cpp`.
% Please note that this  is a personal implementation of an existing work: 
% "Chen, Jia, Chi-Keung Tang, and Jue Wang. Noise brush: interactive high 
% quality image-noise separation. ACM Transactions on Graphics (TOG). ACM, 
% 2009. PLEASE MAKE SURE TO CITE THE ABOVE PAPER AS WELL
%
% 2. go back to the root folder, run `demo.m`
%
% Please refer to our paper for more details

clear all;
close all hidden;

% Add path for all dependencies
folder = fileparts(which('demo.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

cleanImgName    = 'data/ex1_clean.png';
noisyImgName    = 'data/ex1_noisy.png';

denoisedImgName = 'data/ex1_denoised1.png';
[score,psnr]    = denoisingMetrics(noisyImgName, denoisedImgName, cleanImgName);
fprintf('psnr: %f, Our metrics: %f \n', psnr, score);

denoisedImgName = 'data/ex1_denoised2.png';
[score,psnr]    = denoisingMetrics(noisyImgName, denoisedImgName, cleanImgName);
fprintf('psnr: %f, Our metrics: %f \n', psnr, score);

denoisedImgName = 'data/ex1_denoised3.png';
[score,psnr]    = denoisingMetrics(noisyImgName, denoisedImgName, cleanImgName);
fprintf('psnr: %f, Our metrics: %f \n', psnr, score);

denoisedImgName = 'data/ex1_denoised4.png';
[score,psnr]    = denoisingMetrics(noisyImgName, denoisedImgName, cleanImgName);
fprintf('psnr: %f, Our metrics: %f \n', psnr, score);
