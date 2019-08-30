function output = SGM(img, m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Small Gradient Magnitude - compute feature for Small Gradient Magnitude 
% for image denoising quality assessment as defined in the paper:
%
% =========================================================================
% Si Lu, "No-reference Image Denoising Quality Assessment". Computer Vision
% Conference, CVC 2019, Las Vegas.
% =========================================================================
% 
% Version 1.0 (08/10/2019)
%
% This function compute SGM feature for image denoising quality assesment
%
% Inputs:
%
%   img - input uint8 image, only one channel allowed
%   m   - SMG is the the standard deviation of the m(%) smallest non-zero 
%   gradient magnitudes, m should be in range (0,1)
%
% Outputs:
%
%   sGM - computed Small Gradient Magnitude(SMG)
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
[height,width]=size(img);

[GX,GY] = gradient(double(img));
GM      = (GX.*GX+GY.*GY).^0.5;
GM      = reshape(GM,[height*width,1]);
index   = find(GM~=0);
GM2     = GM(index(:));

n1 = round(size(GM2,1)*m);
[smallestNElements, ~] = getNElements(GM2(:), n1);
output = std(smallestNElements);
end

function [smallestNElements, smallestNIdx] = getNElements(A, n)
     [ASorted, AIdx] = sort(A);
     smallestNElements = ASorted(1:n);
     smallestNIdx = AIdx(1:n);
end