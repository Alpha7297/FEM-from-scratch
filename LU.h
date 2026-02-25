#ifndef LU_H
#define LU_H
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>
extern int np;//LU分解分块数，np=16是反复测试的最好结果，32核CPU上分解32767^2稠密矩阵速度可以达到270s
void my_solver(int n, double *A, double *b);//替代法解方程，不需要分配多余内存，求解Ax=b方程，最终x保存在b
#endif