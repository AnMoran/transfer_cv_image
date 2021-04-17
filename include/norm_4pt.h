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

//高斯消元法求解线性方程组
void Gauss(double A[][9], int equ, int var, double* coeffs);

void AB_getPerspectiveTransform(std::vector<std::vector<float>> src, std::vector<std::vector<float>>  dst, 
                                std::vector<std::vector<float> > &M);


#endif