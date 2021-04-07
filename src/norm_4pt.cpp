#include "norm_4pt.h"

//高斯消去法求解线性方程组
void Gauss(double A[][9], int equ, int var, double *coeffs)
{
    int row, col;
    for (row = 0, col = 0; col < var && row < equ; col++, row++)
    {
        int max_r = row;
        for (int i = row + 1; i < equ; i++)
        {
            if ((1e-12) < fabs(A[i][col]) - fabs(A[max_r][col]))
            {
                max_r = i;
            }
        }
        if (max_r != row)
            for (int j = 0; j < var + 1; j++)
                std::swap(A[row][j], A[max_r][j]);
        for (int i = row + 1; i < equ; i++)
        {
            if (fabs(A[i][col]) < (1e-12))
                continue;
            double tmp = -A[i][col] / A[row][col];
            for (int j = col; j < var + 1; j++)
            {
                A[i][j] += tmp * A[row][j];
            }
        }
    }
    for (int i = var - 1; i >= 0; i--)
    { //计算唯一解。
        double tmp = 0;
        for (int j = i + 1; j < var; j++)
        {
            tmp += A[i][j] * (*(coeffs + j));
        }
        coeffs[i] = (A[i][var] - tmp) / A[i][i];
    }
}

/* Calculates coefficients of perspective transformation
 * which maps <quad> into rectangle ((0,0), (w,0), (w,h), (h,0)):
 *
 *      c00*xi + c01*yi + c02
 * ui = ---------------------
 *      c20*xi + c21*yi + c22
 *
 *      c10*xi + c11*yi + c12
 * vi = ---------------------
 *      c20*xi + c21*yi + c22
 *
 * Coefficients are calculated by solving linear system:
 * / x0 y0  1  0  0  0 -x0*u0 -y0*u0 \ /c00\ /u0\
 * | x1 y1  1  0  0  0 -x1*u1 -y1*u1 | |c01| |u1|
 * | x2 y2  1  0  0  0 -x2*u2 -y2*u2 | |c02| |u2|
 * | x3 y3  1  0  0  0 -x3*u3 -y3*u3 |.|c10|=|u3|,
 * |  0  0  0 x0 y0  1 -x0*v0 -y0*v0 | |c11| |v0|
 * |  0  0  0 x1 y1  1 -x1*v1 -y1*v1 | |c12| |v1|
 * |  0  0  0 x2 y2  1 -x2*v2 -y2*v2 | |c20| |v2|
 * \  0  0  0 x3 y3  1 -x3*v3 -y3*v3 / \c21/ \v3/
 *
 * where:
 *   (xi, yi) = (quad[i][0], quad[i][1])
 *        cij - coeffs[i][j], coeffs[2][2] = 1
 *   (ui, vi) - rectangle vertices
 */
void AB_getPerspectiveTransform(std::vector<ABPoint2f> src, std::vector<ABPoint2f> dst,
                                std::vector<std::vector<float> > &M)
{
    double x0 = src[0].x, x1 = src[1].x, x2 = src[3].x, x3 = src[2].x;
    double y0 = src[0].y, y1 = src[1].y, y2 = src[3].y, y3 = src[2].y;
    double u0 = dst[0].x, u1 = dst[1].x, u2 = dst[3].x, u3 = dst[2].x;
    double v0 = dst[0].y, v1 = dst[1].y, v2 = dst[3].y, v3 = dst[2].y;
    double A[8][9] = {
        {x0, y0, 1, 0, 0, 0, -x0 * u0, -y0 * u0, u0},
        {x1, y1, 1, 0, 0, 0, -x1 * u1, -y1 * u1, u1},
        {x2, y2, 1, 0, 0, 0, -x2 * u2, -y2 * u2, u2},
        {x3, y3, 1, 0, 0, 0, -x3 * u3, -y3 * u3, u3},
        {0, 0, 0, x0, y0, 1, -x0 * v0, -y0 * v0, v0},
        {0, 0, 0, x1, y1, 1, -x1 * v1, -y1 * v1, v1},
        {0, 0, 0, x2, y2, 1, -x2 * v2, -y2 * v2, v2},
        {0, 0, 0, x3, y3, 1, -x3 * v3, -y3 * v3, v3},
    };
    double *coeffs = new double[8];
    Gauss(A, 8, 8, coeffs);
    // *(coeffs + 8) = 1;
    // std::cout<<*(coeffs)<<"\n";
    M[0][0] = (float)(*coeffs);
    M[0][1] = (float)coeffs[1];
    M[0][2] = (float)coeffs[2];
    M[1][0] = (float)coeffs[3];
    M[1][1] = (float)coeffs[4];
    M[1][2] = (float)coeffs[5];
    M[2][0] = (float)coeffs[6];
    M[2][1] = (float)coeffs[7];
    M[2][2] = 1;
    delete []coeffs;
}

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
