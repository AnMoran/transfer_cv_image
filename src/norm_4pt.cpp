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
void AB_getPerspectiveTransform(std::vector<std::vector<float>> src, std::vector<std::vector<float>> dst,
                                std::vector<std::vector<float> > &M)
{
    double x0 = src[0][0], x1 = src[1][0], x2 = src[3][0], x3 = src[2][0];
    double y0 = src[0][1], y1 = src[1][1], y2 = src[3][1], y3 = src[2][1];
    double u0 = dst[0][0], u1 = dst[1][0], u2 = dst[3][0], u3 = dst[2][0];
    double v0 = dst[0][1], v1 = dst[1][1], v2 = dst[3][1], v3 = dst[2][1];
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



