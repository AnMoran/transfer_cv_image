#ifndef _COMMON_H_
#define _COMMON_H_

#include <iostream>
#include <vector>

void RGB2float(const unsigned char *RGBImage, int width, int height, float *dstImg);

void Float2RGB(const float *floatImage, int width, int height, unsigned char *dstImg);

void warpAffineLinear(const float *pfSrcImg, float *pfDstImg,
                      int nSrcWidth, int nSrcHeight, 
                      int nDstWidth,int nDstHeight, int nChannels, std::vector<std::vector<float> > M);
#endif