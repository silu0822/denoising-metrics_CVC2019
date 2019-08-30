function output = GH(noisyImg, denoisedImg, nSigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gradient Histogram Preservation - compute feature for Gradient Histogram 
% Preservation for image denoising quality assessment as defined in paper:
%
% =========================================================================
% Si Lu, "No-reference Image Denoising Quality Assessment". Computer Vision
% Conference, CVC 2019, Las Vegas.
% =========================================================================
% 
% Version 1.0 (08/10/2019)
%
% This function compute GH feature for image denoising quality assesment
%
% Inputs:
%
%   noisyImg    - uint8 noisy image to be denoised, only 1 channel allowed
%   denoisedImg - uint8 denoised image, only 1 channel allowed
%   nSigma      - sigma of the noisy image
%
% Outputs:
%
%   GH     - computed Gradient Histogram Preservation(GH)
% 
% See demo.m in this same code package for examples on how to use this
% function to compute features for denoising quality assessment.
%
% License: Only for academic use. All rights reserved to the authors of the
% paper (Si Lu). If you have any questions, comments or suggestions please 
% contact Si Lu at lusi@pdx.edu or daniel08220822@gmail.com. Please cite 
% our paper if you use this piece of code.
%
% =========================================================================
% Si Lu, "No-reference Image Denoising Quality Assessment". Computer Vision
% Conference, CVC 2019, Las Vegas.
% =========================================================================
%
%
% Please also refer to the paper below if you use our code as our feature 
% is derived from this work and uses their implementation at:
% https://github.com/tingfengainiaini/GHPBasedImageRestoration
% =========================================================================
% Zuo, W., Zhang, L., Song, C., Zhang, D.: Texture enhanced image denoising 
% via gradient histogram preservation. In: IEEE Conference on Computer 
% Vision and Pattern Recognition (CVPR), pp. 1203–1210 (2013)
% =========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noisyImg    = double(noisyImg);
denoisedImg = double(denoisedImg);
[h,w,c] = size(noisyImg);
[a, b]   =   HistEst(noisyImg, nSigma);
            
x = 0:1:255;

[Gx,Gy] = gradient(denoisedImg);
GM      = Gx.*Gx+Gy.*Gy;
GM      = GM.^0.5;
denoiseHist = hist(GM(:),x);
denoiseHist = denoiseHist/sum(denoiseHist(:));
histDif     = abs(denoiseHist-a);
output      = sum(histDif(:));