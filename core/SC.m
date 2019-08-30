function [ score ] = SC( image_noisy, image_denoised, patchSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Adaptive Structure Correlation - compute feature for Adaptive Structure
% Correlation for image denoising quality assessment as defined in paper:
%
% =========================================================================
% Si Lu, "No-reference Image Denoising Quality Assessment". Computer Vision
% Conference, CVC 2019, Las Vegas.
% =========================================================================
% 
% Version 1.0 (08/10/2019)
%
% This function compute SC feature for image denoising quality assesment
%
% Inputs:
%
%   image_noisy    - uint8 noisy image to be denoised, only 1 channel allowed
%   image_denoised - uint8 denoised image, only 1 channel allowed
%   patchSize      - patch size, default 8
%
% Outputs:
%
%   SC     - computed Adaptive Structure Correlation(SC)
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
% is derived from:
% =========================================================================
% Kong, X., Li, K., Yang, Q., Wenyin, L., Yang, M.H.: A new image quality
% metric for image auto-denoising. In: IEEE International Conference on 
% Computer Vision (ICCV), pp. 2888–2895 (2013)
% =========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


img_m = image_noisy - image_denoised;    %Compute the "Method Noise Image" (MNI)
if(image_noisy==image_denoised)
    score=1;
else
    [N] = SSIM(image_noisy,img_m, patchSize);  %Compute the map of noise reduction measurement N (see Eq. 3 in the paper)
    
    [P] = SSIM(image_noisy,image_denoised, patchSize);  %Compute the map of Structure preservation measurement P (see Eq. 4 in the paper)
    
    score = -corr(real(N(:)),real(P(:)),'type','Pearson'); %Compute linear correlation of N and P and normalize it as the quality score. final score ranges [-1,1]
end
end

