#ifndef _NORM_4_PT_H_
#define _NORM_4_PT_H_
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <climits>
#include <sstream>
#include <memory>

typedef struct tagABPoint2f
{
    double x;
    double y;
} ABPoint2f ;

//高斯消元法求解线性方程组
void Gauss(double A[][9], int equ, int var, double* coeffs);

void AB_getPerspectiveTransform(std::vector<ABPoint2f> src, std::vector<ABPoint2f>  dst, 
                                std::vector<std::vector<float> > &M);

void RGB2float(const unsigned char *RGBImage, int width, int height, float *dstImg);

void Float2RGB(const float *floatImage, int width, int height, unsigned char *dstImg);

void warpAffineLinear(const float *pfSrcImg, float *pfDstImg,
                      int nSrcWidth, int nSrcHeight, int nDstWidth,
                      int nDstHeight, int nChannels, std::vector<std::vector<float> > M);



#endif