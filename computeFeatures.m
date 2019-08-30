function features = computeFeatures(nImgName, dImgName)

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
%
% Please refer to our paper for more details

features = zeros(1,19);
nImage   = imread(nImgName);
dImage   = imread(dImgName);

%% Compute feature SR (Please compile the SR binary file first using `mex srMex.cpp`!)
[~, SR1] = srMex(double(nImage), double(dImage), 1); features(1,1) = SR1;
[~, SR2] = srMex(double(nImage), double(dImage), 2); features(1,1) = SR2;
[~, SR3] = srMex(double(nImage), double(dImage), 3); features(1,1) = SR3;

%% Compute feature SC, not our work, please cite SC paper as well. 
%  See util/metric.m for more details
features(1,4) = SC(double(nImage), double(dImage), 6);
features(1,5) = SC(double(nImage), double(dImage), 8);
features(1,6) = SC(double(nImage), double(dImage), 10);

%% Compute feature VR
features(1,7)  = VR(double(nImage), double(dImage), 0.5, 1, 1);
features(1,8)  = VR(double(nImage), double(dImage), 1.0, 1, 1);
features(1,9)  = VR(double(nImage), double(dImage), 0.5, 2, 1);
features(1,10) = VR(double(nImage), double(dImage), 1.0, 2, 1);
features(1,11) = VR(double(nImage), double(dImage), 0.5, 2, 2);
features(1,12) = VR(double(nImage), double(dImage), 1.0, 2, 2);

%% Compute feature GH
nSigma = 20;
features(1,13) = GH(double(nImage), double(dImage), nSigma);

%% Compute feature SMG
features(1,14) = SGM(dImage, 0.40);
features(1,15) = SGM(dImage, 0.50);
features(1,16) = SGM(dImage, 0.60);

%% Compute feature SS 
features(1,17) = SS(dImage, 0.97);
features(1,18) = SS(dImage, 0.98);
features(1,19) = SS(dImage, 0.99);
