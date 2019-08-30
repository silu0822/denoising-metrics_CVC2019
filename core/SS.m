function output = SS(img, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% localSelfSimilarity - compute feature for local self-similarity for image
% denoising quality assessment as defined in the paper:
%
% =========================================================================
% Si Lu, "No-reference Image Denoising Quality Assessment". Computer Vision
% Conference, CVC 2019, Las Vegas.
% =========================================================================
% 
% Version 1.0 (08/10/2019)
%
% This function compute SS feature for image denoising quality assesment
%
% Inputs:
%
%   img   - input uint8 image, only one channel allowed
%   alpha - recvoery percentage to estimate SS, from 0 to 1 [0,1)
%
% Outputs:
%
%   SS    - computed localSelfSimilarity(SS)
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

img     = double(img);
[h,w,c] = size(img);
if (c~=1)
    fprintf('Input image must be uint8 grayscale image');
    output = 0;
    return;
end

step        = 15;
blockCount  = 0;
blockMatrix = zeros(step*step, uint16(h/step)*uint16(w/step));
for i=1 : step : h - 2*step
    for j=1 : step : w - 2*step
        blockCount = blockCount+1;
        currBlock  = img(i+1 : i+step, j+1:j+step);
        currBlock  = reshape(currBlock, step*step, 1);
        blockMatrix(:,blockCount) = currBlock(:,1);
    end
end
S = svd(blockMatrix);  S=S/sum(S(:));

C=size(S,1);
sc=0;
for t=2:C
    S(t,1)=S(t-1,1)+S(t,1);
end
for t=1:C
    if(S(t,1)>alpha)
        sc=t;
        break;
    end
end
output=sc/C;