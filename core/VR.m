function output = VR(noisyImg, denoisedImg, lambda, ld, ls)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variational Denoising Residual - compute feature for Variational 
% Denoising Residual for image denoising quality assessment as defined in
% paper:
%
% =========================================================================
% Si Lu, "No-reference Image Denoising Quality Assessment". Computer Vision
% Conference, CVC 2019, Las Vegas.
% =========================================================================
% 
% Version 1.0 (08/10/2019)
%
% This function compute VR feature for image denoising quality assesment
%
% Inputs:
%
%   noisyImg    - uint8 noisy image to be denoised, only 1 channel allowed
%   denoisedImg - uint8 denoised image, only 1 channel allowed
%   lambda      - balancing factor between data term and smoothness term
%   ld          - 1 or 2, means L1 or L2 norm used in data term
%   ls          - 1 or 2, means L1 or L2 norm used in smoothness term
%
% Outputs:
%
%   VR     - computed Variational Denoising Residual(VR)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noisyImg    = double(noisyImg);
denoisedImg = double(denoisedImg);
[h,w,c] = size(noisyImg);
if (c~=1)
    fprintf('Input image must be uint8 grayscale image');
    output = 0;
    return;
end

if( (ld~=1&&ld~=2) || (ls~=1&&ls~=2) )
    fprintf('ld and ls should be either 1 or 2');
    output = 0;
    return;
end

if (lambda<=0)
    fprintf('lambda should be positive');
    output = 0;
    return;
end

if(ld==1)
    imgDif = abs(noisyImg - denoisedImg);
else
    imgDif = (noisyImg - denoisedImg).^2; 
end
dataTerm = sum(imgDif(:));

[Gmag,~] = imgradient(denoisedImg);
if(ls==1)
    Gmag_norm = abs(Gmag);
else
    Gmag_norm = Gmag.^2;
end
smoothnessTerm = sum(Gmag_norm(:));

output = (dataTerm + lambda * smoothnessTerm)/w/h;
