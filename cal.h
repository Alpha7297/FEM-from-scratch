#ifndef CAL_H
#define CAL_H
#include "geometry.hpp"
#include "LU.h"
#include<cstdio>
#include<cstdlib>
#include<algorithm>
#include<omp.h>
extern int np;
typedef struct PHI{
    std::vector<double> phi;
}PHI;
PHI dirc_cal_phi(double (*g)(point),double (*rho)(point),std::vector<point> p,std::vector<triangle> tri);//迪利克雷边界求解
//输入分别为 边界电势 负的内部电荷(nabla^2 phi=rho) 所有点集 三角形集
PHI neum_cal_phi(double (*g)(point),double (*rho)(point),std::vector<point> p,std::vector<triangle> tri);//诺伊曼边界求解
//输入分别为 边界法向电场 负的内部电荷(nabla^2 phi=rho) 所有点集 三角形集
//没有检查通量和电荷总和关系，不一致会奇异
#endif