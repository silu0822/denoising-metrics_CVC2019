/** ------------------------------------------- **/
/** - Structural Residual (SR) implementation - **/
/** ------------------------------------------- **/
/**
 * Author: Si Lu, Portland State University
 *
 * Description: This is a personal implementation of an existing paper: "Chen, 
 *   Jia, Chi-Keung Tang, and Jue Wang. Noise brush: interactive high quality
 *   image-noise separation. ACM Transactions on Graphics (TOG). ACM, 2019"
 *
 * Reference: We use this implementation of to compute feature Structural
 *   Residual for image denoising quality assessment as defined our paper:
 *   "Si Lu, No-reference Image Denoising Quality Assessment. Computer Vision
 *   Conference, CVC 2019, Las Vegas."
 *
 * How to run:
 *   1. In MATLAB, run `mex srMex.cpp` to compile
 *   2. Usage in matlab: [srMap, SR] = srMex(nImage, dImage, index);
 *     Input:
 *       nImage: noisy input image in double format
 *       dImage: denoised image in double format
 *       index:  1, 2 or 3 indicating which set of parameter to choose
 *     Output:
 *       srMap:  a 2D image with the same size as input images, represent SR
 *               values at each pixel
 *       SR:     Final output SR value
 **/

#include <mex.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <vector>
#include <numeric>
#include <cmath>

double matMut(double a, double b, double D1, double D2, double D3, double D4)
{
    return a*a*D1 + a*b*(D2 + D3) + b*b*D4;
}

void symetrize(
    const std::vector<double> &img
,   std::vector<double> &img_sym
,   const unsigned width
,   const unsigned height
,   const unsigned chnls
,   const unsigned N
){
    //! Declaration
    const unsigned w = width + 2 * N;
    const unsigned h = height + 2 * N;

    if (img_sym.size() != w * h * chnls)
        img_sym.resize(w * h * chnls);

    for (unsigned c = 0; c < chnls; c++)
    {
        unsigned dc = c * width * height;
        unsigned dc_2 = c * w * h + N * w + N;

        //! Center of the image
        for (unsigned i = 0; i < height; i++)
            for (unsigned j = 0; j < width; j++, dc++)
                img_sym[dc_2 + i * w + j] = img[dc];

        //! Top and bottom
        dc_2 = c * w * h;
        for (unsigned j = 0; j < w; j++, dc_2++)
            for (unsigned i = 0; i < N; i++)
            {
                img_sym[dc_2 + i * w] = img_sym[dc_2 + (2 * N - i - 1 + 1) * w];
                img_sym[dc_2 + (h - i - 1) * w] = img_sym[dc_2 + (h - 2 * N + i -1) * w];
            }

        //! Right and left
        dc_2 = c * w * h;
        for (unsigned i = 0; i < h; i++)
        {
            const unsigned di = dc_2 + i * w;
            for (unsigned j = 0; j < N; j++)
            {
                img_sym[di + j] = img_sym[di + 2 * N - j - 1 + 1];
                img_sym[di + w - j - 1] = img_sym[di + w - 2 * N + j -1];
            }
        }
    }
    return;
}

void mexFunction(int nlhs, mxArray *plhs[],          // output
                 int nrhs, const mxArray *prhs[])    // input
{
    // Input: n_image (double), d_image (double), SR_ind
    // Output: SR_map, SR_score
	// Usage in MATLAB: [srMap, SR] = srMex(nImage, dImage, index); 
    const mwSize* dims;
    double *rawNoisyImage, *rawDenoisedImage;
    
    // Argument checking
    if (nlhs!=2 || nrhs != 3) {
        mexErrMsgTxt("Two ouutput arguments and three input arguments are required.") ;
        mexErrMsgTxt("Usage: [SR_map, SR] = srMex(noisyImage, denoisedImage, SR_ind).") ;
    }

    // load input
    int numElements    = (int)mxGetNumberOfElements(prhs[0]);
    dims               = mxGetDimensions(prhs[0]) ;
    rawNoisyImage      = (double*)mxGetData(prhs[0]) ;//mxGetData returns a void pointer, so cast it
    rawDenoisedImage   = (double*)mxGetData(prhs[1]) ;//mxGetData returns a void pointer, so cast it
    const int width   = dims[1]; 
    const int height  = dims[0]; //Note: first dimension provided is height and second is width
    const int imageSize = width * height;

    std::vector<double> noisyImage(imageSize);
    std::vector<double> denoisedImage(imageSize);
    std::vector<double> noise(imageSize);

    // select parameter settings
	const double SRInd = mxGetScalar(prhs[2]);
    double sigma_d, sigma_c;
    if (SRInd == 1)
    {
        sigma_d = 1.0;  sigma_c = 4.0;    
    }
    else if (SRInd == 2)
    {
        sigma_d = 4.0;  sigma_c = 10.0;   
    }else
    {
        sigma_d = 10.0; sigma_c = 30.0;   
    }

    // reading data from column-major MATLAB matrics to row-major C matrices
    // (i.e perform transpose)
    int x, y, ii;
    for(x = 0, ii = 0; x < width; x++)
    {
        for(y = 0; y < height; y++)
        {
            int i = (int)(y*width+x);
            noisyImage[i]    = rawNoisyImage[ii];
            denoisedImage[i] = rawDenoisedImage[ii];
            noise[i]         = noisyImage[i] - denoisedImage[i];
            ii++;
        }
    }

    // internal parameters
    const int templateWindowSize = 7;
    const int searchWindowSize   = 7;
    const int tr          = templateWindowSize >> 1;
    const int sr          = searchWindowSize >> 1;
    const int bb          = sr + tr;
    const int D           = searchWindowSize*searchWindowSize;
    const int height_b    = height + 2 * bb;
    const int width_b     = width  + 2 * bb;
    const int imageSize_b = height_b * width_b;

    // Create large size image for bounding box;
    std::vector<double> nImg(imageSize_b, 0);
    std::vector<double> dImg(imageSize_b, 0);
    std::vector<double> nse (imageSize_b, 0);

    symetrize(noisyImage,    nImg, width, height, 1, bb);
    symetrize(denoisedImage, dImg, width, height, 1, bb);
    symetrize(noise,         nse, width, height, 1, bb);

    // Calculate Gradient: Gx Gy
    std::vector<double> Gx(imageSize, 0.0);
	std::vector<double> Gy(imageSize, 0.0);
	
    for (int i = 0; i<height; i++)
    {
        for (int j = 0; j<width; j++)
        {
            int index_b = ( i + bb ) * width_b + ( j + bb );
            Gx [i*width+j] = (dImg[index_b+1]-dImg[index_b-1])/2.0;
            Gy [i*width+j] = (dImg[index_b+width_b]-dImg[index_b-width_b])/2.0;
        }
    }
    
    // Create large size gradient for bounding box;
    std::vector<double> Gxb(imageSize_b);
    std::vector<double> Gyb(imageSize_b);
    symetrize(Gx, Gxb, width, height, 1, bb);
    symetrize(Gy, Gyb, width, height, 1, bb);

    // Calculate Dp_xx, Dp_yy and Dp_xy
    std::vector<double> Dp_xx(imageSize, 0.0), Dp_yy(imageSize, 0.0), Dp_xy(imageSize, 0.0);

    int ptr = 0; 
    for (int i = 0; i<height; i++)
    {
        for (int j = 0; j<width; j++)
        {    
            for (int I = -sr; I <= sr; I++)
            {
                for (int J = -sr; J <= sr; J++)
                {
                    double tmpGx = Gxb[(i+bb+I)*width_b + j+bb+J];
                    double tmpGy = Gyb[(i+bb+I)*width_b + j+bb+J];
                    Dp_xx[ptr]   = Dp_xx[ptr] + tmpGx*tmpGx;
                    Dp_yy[ptr]   = Dp_yy[ptr] + tmpGy*tmpGy;
                    Dp_xy[ptr]   = Dp_xy[ptr] + tmpGx*tmpGy;
                }
            }
            Dp_xx[ptr] = Dp_xx[ptr]/(double)D;
            Dp_yy[ptr] = Dp_yy[ptr]/(double)D;
            Dp_xy[ptr] = Dp_xy[ptr]/(double)D;
            ptr++;
        }
    }

    std::vector<double> Dp_xxB(imageSize_b, 0);
    std::vector<double> Dp_yyB(imageSize_b, 0);
    std::vector<double> Dp_xyB(imageSize_b, 0);

    symetrize(Dp_xx, Dp_xxB, width, height, 1, bb);
    symetrize(Dp_yy, Dp_yyB, width, height, 1, bb);
    symetrize(Dp_xy, Dp_xyB, width, height, 1, bb);

    // Spatial distance weight
	std::vector<double> sw((sr*sr * 2) + 1, 0);
    double dis = 0, v = 0;
    for (int i = 0; i <= sr; i++)
    {
        for (int j = i; j <= sr; j++)
        {
            dis = pow(double(i), 2) + pow(double(j), 2);
            v = std::exp( -1.0*dis / sigma_d / sigma_d );
            sw[i*i + j*j] = v;
        }
    }

    // Color weight
	std::vector<double> w(256 * 256, 0);
    int emax = 0;
    for (int i = 0; i < 256 * 256 ; i++)
    {
        double v = std::exp(-1.0*i / sigma_c / sigma_c);
        w[i] = v;
        if (v<0.001)
        {
            emax = i;
            break;
        }
    }
    for (int i = emax; i < 256 * 256; i++) w[i] = 0.0;

    double ws = 0, wc = 0, wd = 0;
    int cdif = 0;

    const double sigma_sd = 10.0;
    double a, b, D1, D2, D3, D4;
    double c, d, D5, D6, D7, D8;
    
	double SR = 0;
	std::vector<double> SR_map(imageSize,0);
    for (int i = 0; i<height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            double cu_w = 0, cu_sum = 0;
			D1 = Dp_xxB[(i+bb)*width_b+j+bb];
			D2 = Dp_xyB[(i+bb)*width_b+j+bb];
			D3 = D2;
			D4 = Dp_yyB[(i+bb)*width_b+j+bb];
			
			//search current pixel's neighborhood
            for (int l = -sr; l <= sr; l++)
            {
                for (int k = -sr; k <= sr; k++)
                {
                    //calculate ws
                    ws = sw[l*l + k*k];
                    //calculate wc
					cdif = (int)(dImg[(i+bb)*width_b+j+bb] - dImg[(i+bb+l)*width_b+j+bb+k]);
                    wc = w[cdif*cdif];
                    //calculate wd
                    a = -k; b = -l;
                    double tmpwd = matMut(a, b, D1, D2, D3, D4);
                    wd = std::exp(-1 * tmpwd / sigma_sd / sigma_sd);

                    c = k; d = l;
					D5    = Dp_xxB[(i+bb+l)*width_b+j+bb+k];
					D6    = Dp_xyB[(i+bb+l)*width_b+j+bb+k];
					D7    = D6;
					D8    = Dp_yyB[(i+bb+l)*width_b+j+bb+k];
                    tmpwd = matMut(c, d, D5, D6, D7, D8);
                    wd    = wd + std::exp(-1 * tmpwd / sigma_sd / sigma_sd);
                    wd    = wd / 2;

                    cu_w   = cu_w + ws*wc*wd;
					cu_sum = cu_sum + ws*wc*wd*(nse[(i+bb+l)*width_b+j+bb+k]);
					int a = 0;
                }
            }
			
			if(cu_w==0)
			{
			    continue;
            }			
            double SR_ij = std::pow(cu_sum / (cu_w + 0.00001),2);
            SR += SR_ij;
            SR_map[(int)(i*width+j)] = SR_ij;
					
        }
		
    }
	
    SR = pow(SR / height / width, 0.5);
	
	// output
    double* outSRMap;
    plhs[0] = mxCreateNumericMatrix(height,width,mxDOUBLE_CLASS,mxREAL);
    outSRMap = (double*)mxGetData(plhs[0]);
	
	// copying data from row-major C matrix to column-major MATLAB matrix (i.e. perform transpose)
    for(x = 0, ii = 0; x < width; x++)
    {
        for(y = 0; y < height; y++)
        {
            int i0 = y*width+x;
            outSRMap[ii] = SR_map[i0];
            ii++;
        }
    }
    plhs[1] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
    double *SR_score = (double*)mxGetData(plhs[1]);//gives a void*, cast it to int*
    *SR_score = SR;

	return;
}