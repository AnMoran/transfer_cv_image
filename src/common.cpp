#include "common.h"

//rgbrgbrgb……图转换为rrrrrgggggbbbbb的float数组
void RGB2float(const unsigned char *RGBImage, int width, int height, float *dstImg)
{
    for (int i = 0; i < width * height; i++)
    {
        dstImg[i] = RGBImage[3 * i];
        dstImg[i + width * height] = RGBImage[3 * i + 1];
        dstImg[i + 2 * width * height] = RGBImage[3 * i + 2];
    }
}

//rrrrrgggggbbbbb的float数组转换为rgbrgbrgb……图
void Float2RGB(const float *floatImage, int width, int height, unsigned char *dstImg)
{
    for (int i = 0; i < width * height; i++)
    {
        dstImg[3 * i] = floatImage[i];
        dstImg[3 * i + 1] = floatImage[i + width * height];
        dstImg[3 * i + 2] = floatImage[i + 2 * width * height];
    }
}

//透视变换，输入多通道数据是rrrrrrggggggbbbbb，不是rgbrgbrgb
void warpAffineLinear(const float *pfSrcImg, float *pfDstImg,
                      int nSrcWidth, int nSrcHeight, 
                      int nDstWidth,int nDstHeight, int nChannels, std::vector<std::vector<float> > M)
{
    int i, j;
    float row = 0.0f;
    float col = 0.0f;
    float *pDstImg = pfDstImg;
    int m = (nSrcHeight - 1) * nSrcWidth;

    float D = M[0][0] * M[1][1] - M[0][1] * M[1][0];
    D = D != 0 ? 1.0f / D : 0.0f;
    float A11 = M[1][1] * D, A22 = M[0][0] * D;
    M[0][0] = A11;
    M[0][1] *= -D;
    M[1][0] *= -D;
    M[1][1] = A22;
    float b1 = -M[0][0] * M[0][2] - M[0][1] * M[1][2];
    float b2 = -M[1][0] * M[0][2] - M[1][1] * M[1][2];
    M[0][2] = b1;
    M[1][2] = b2;

    memset(pDstImg, 0, sizeof(float) * nDstWidth * nDstHeight * nChannels);

    int nDstWH = nDstWidth * nDstHeight;
    int nSrcWH = nSrcWidth * nSrcHeight;
    float *pDstImg2 = pfDstImg + nDstWH;
    float *pDstImg3 = pDstImg2 + nDstWH;
    for (i = 0; i < nDstHeight; i++)
    {
        for (j = 0; j < nDstWidth; j++)
        {
            col = M[0][0] * j + M[0][1] * i + M[0][2];
            row = M[1][0] * j + M[1][1] * i + M[1][2];
            if (row < 0 || col < 0 || row > nSrcHeight - 1 || col > nSrcWidth - 1)
            {
                pDstImg++;
                pDstImg2++;
                pDstImg3++;
                continue;
            }
            else if (row < nSrcHeight - 1)
            {
                float x = row - (int)row;
                float y = col - (int)col;
                float a = (1 - x) * (1 - y);
                float b = (1 - x) * y;
                float c = x * (1 - y);
                float d = x * y;
                int offset = (int)row * nSrcWidth + (int)col;
                float *pSrcImg = (float *)pfSrcImg + offset;
                float *pSrcImg2 = (float *)pfSrcImg + offset + nSrcWH;
                float *pSrcImg3 = (float *)pfSrcImg + offset + nSrcWH * 2;
                *pDstImg++ = (float)((float)pSrcImg[0] * a + (float)pSrcImg[1] * b + (float)pSrcImg[nSrcWidth] * c + (float)pSrcImg[nSrcWidth + 1] * d);
                *pDstImg2++ = (float)((float)pSrcImg2[0] * a + (float)pSrcImg2[1] * b + (float)pSrcImg2[nSrcWidth] * c + (float)pSrcImg2[nSrcWidth + 1] * d);
                *pDstImg3++ = (float)((float)pSrcImg3[0] * a + (float)pSrcImg3[1] * b + (float)pSrcImg3[nSrcWidth] * c + (float)pSrcImg3[nSrcWidth + 1] * d);
            }
            else
            {
                *pDstImg++ = (float)(pfSrcImg[m + (int)col]);
                *pDstImg2++ = (float)(pfSrcImg[m + (int)col + nSrcWH]);
                *pDstImg2++ = (float)(pfSrcImg[m + (int)col + nSrcWH * 2]);
            }
        }
    }
}
