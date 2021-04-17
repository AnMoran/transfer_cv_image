#ifndef _NORM_5_PT_H
#define _NORM_5_PT_H
#include<memory.h>
#include<iostream>
#include<vector>
#include<stdlib.h>
#include<math.h>
#include<iomanip>

#ifndef __NORM_WIDTH
#define NORM_WIDTH 112
#endif

#ifndef __NORM_HEIGHT
#define NORM_HEIGHT 112
#endif

int NormImg_5_points(const float* pImg, int nImgWidth, int nImgHeight,int channels,
	std::vector<std::vector<float>> src,std::vector<std::vector<float>> dst, float *pfNormImg);

bool center_and_normalize_points(std::vector<std::vector<float>> points, std::vector<std::vector<float>> &points_new, 
	std::vector<std::vector<float>> &point_matrix);

bool CalMatDot(std::vector<std::vector<float>> A, std::vector<std::vector<float>> B, std::vector<std::vector<float>>&C);

float Creat_M(float *p, int m, int n, int k);

bool Gauss(std::vector<std::vector<float>>A, std::vector<std::vector<float> > &B);

float get_norm(float *x, int n);

float MatDet(float *p, int n);

float normalize(float *x, int n);

void orth(float *a, float *b, int n);

inline float product(float*a, float *b, int n);

bool svd(std::vector<std::vector<float>> A, int K, std::vector<std::vector<float>> &U, std::vector<float> &S, 
	std::vector<std::vector<float>> &V);

void print(std::vector<std::vector<float>> &A);

#endif
