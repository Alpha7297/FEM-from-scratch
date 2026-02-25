#include "geometry.hpp"
#include "LU.h"
#include<cstdio>
#include<cstdlib>
#include<algorithm>
#include<functional>
#include<omp.h>
extern int np;
using std::vector;
//储存结果
typedef struct PHI{
    vector<double> phi;
}PHI;
double triangle_integral(std::function<double(point)> rho,triangle a){
    point p1=(a.a*1/6+a.b*1/6+a.c*2/3);
    point p2=(a.a*1/6+a.b*2/3+a.c*1/6);
    point p3=(a.a*2/3+a.b*1/6+a.c*1/6);
    double res=fabs(a.d)/6*(rho(p1)+rho(p2)+rho(p3));
    return res;
}
PHI dirc_cal_phi(double (*g)(point),double (*rho)(point),vector<point> p,vector<triangle> tri){
    //读取
    for(int i=0;i<tri.size();i++){
        triangle curr=tri[i];
        p[curr.ia].tri_idx.push_back(i);
        p[curr.ib].tri_idx.push_back(i);
        p[curr.ic].tri_idx.push_back(i);
    }
    int N=p.size();
    double* K=(double*)calloc((N+5)*(N+5),sizeof(double));//K phi=f
    double* f=(double*)calloc(N+5,sizeof(double));
    int len=0;
    for(;len<p.size();len++){
        if(!p[len].edge){
            break;
        }
    }
    //并行对每一个内部点的电势求偏导，得到矩阵方程
    //E=\sum_i S_i(p_i1^2+p_i2^2)
    #pragma omp parallel for
    for(int i=len;i<p.size();i++){
        for(const auto&idx:p[i].tri_idx){
            triangle curr=tri[idx];
            if(i==curr.ia){
                auto temp=[rho,curr](point a)->double{
                    return rho(a)*(curr.p1[0]*a.x+curr.p2[0]*a.y+curr.p3[0]);
                };
                f[i-len]-=triangle_integral(temp,curr);
                K[(i-len)*(N-len)+i-len]+=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[0]+curr.p2[0]*curr.p2[0]);
                if(p[curr.ib].edge){
                    double phi0=g(p[curr.ib]);
                    f[i-len]-=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[1]*phi0+curr.p2[0]*curr.p2[1]*phi0);
                }
                else{
                    K[(i-len)*(N-len)+curr.ib-len]+=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[1]+curr.p2[0]*curr.p2[1]);
                }
                if(p[curr.ic].edge){
                    double phi0=g(p[curr.ic]);
                    f[i-len]-=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[2]*phi0+curr.p2[0]*curr.p2[2]*phi0);
                }
                else{
                    K[(i-len)*(N-len)+curr.ic-len]+=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[2]+curr.p2[0]*curr.p2[2]);
                }
            }
            else if(i==curr.ib){
                auto temp=[rho,curr](point a)->double{
                    return rho(a)*(curr.p1[1]*a.x+curr.p2[1]*a.y+curr.p3[1]);
                };
                f[i-len]-=triangle_integral(temp,curr);
                K[(i-len)*(N-len)+i-len]+=fabs(curr.d)/2.0*(curr.p1[1]*curr.p1[1]+curr.p2[1]*curr.p2[1]);
                if(p[curr.ia].edge){
                    double phi0=g(p[curr.ia]);
                    f[i-len]-=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[1]*phi0+curr.p2[0]*curr.p2[1]*phi0);
                }
                else{
                    K[(i-len)*(N-len)+curr.ia-len]+=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[1]+curr.p2[0]*curr.p2[1]);
                }
                if(p[curr.ic].edge){
                    double phi0=g(p[curr.ic]);
                    f[i-len]-=fabs(curr.d)/2.0*(curr.p1[1]*curr.p1[2]*phi0+curr.p2[1]*curr.p2[2]*phi0);
                }
                else{
                    K[(i-len)*(N-len)+curr.ic-len]+=fabs(curr.d)/2.0*(curr.p1[1]*curr.p1[2]+curr.p2[1]*curr.p2[2]);
                }
            }
            else{
                auto temp=[rho,curr](point a)->double{
                    return rho(a)*(curr.p1[2]*a.x+curr.p2[2]*a.y+curr.p3[2]);
                };
                f[i-len]-=triangle_integral(temp,curr);
                K[(i-len)*(N-len)+i-len]+=fabs(curr.d)/2.0*(curr.p1[2]*curr.p1[2]+curr.p2[2]*curr.p2[2]);
                if(p[curr.ib].edge){
                    double phi0=g(p[curr.ib]);
                    f[i-len]-=fabs(curr.d)/2.0*(curr.p1[2]*curr.p1[1]*phi0+curr.p2[2]*curr.p2[1]*phi0);
                }
                else{
                    K[(i-len)*(N-len)+curr.ib-len]+=fabs(curr.d)/2.0*(curr.p1[2]*curr.p1[1]+curr.p2[2]*curr.p2[1]);
                }
                if(p[curr.ia].edge){
                    double phi0=g(p[curr.ia]);
                    f[i-len]-=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[2]*phi0+curr.p2[0]*curr.p2[2]*phi0);
                }
                else{
                    K[(i-len)*(N-len)+curr.ia-len]+=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[2]+curr.p2[0]*curr.p2[2]);
                }
            }
        }
    }
    //LU分解
    my_solver(N-len,K,f);
    PHI res;
    for(int i=0;i<len;i++){
        //加入边界点电势
        res.phi.push_back(g(p[i]));
    }
    for(int i=0;i<N-len;i++){
        //加入内部点电势
        res.phi.push_back(f[i]);
    }
    free(K);
    free(f);
    return res;
}
double integral(std::function<double(point)> g,point a,point b) {
    double L=distance(a,b);
    point mid=(a+b)/2;
    point dir=(b-a)/L;
    double t1=-1.0/sqrt(3);
    double t2=1.0/sqrt(3);
    point p1=mid+dir*(L/2*t1);
    point p2=mid+dir*(L/2*t2);
    return L/2*(g(p1)+g(p2));
}
PHI neum_cal_phi(double (*g)(point),double (*rho)(point),vector<point> p,vector<triangle> tri){
    for(int i=0;i<tri.size();i++){
        triangle curr=tri[i];
        p[curr.ia].tri_idx.push_back(i);
        p[curr.ib].tri_idx.push_back(i);
        p[curr.ic].tri_idx.push_back(i);
    }
    int N=p.size();
    double* K=(double*)calloc((N+5)*(N+5),sizeof(double));//K phi=f
    double* f=(double*)calloc(N+5,sizeof(double));
    //并行对每一个内部点的电势求偏导，得到矩阵方程
    //固定最后一个节点电势为0
    #pragma omp parallel for
    for(int i=0;i<N-1;i++){
        for(const auto&idx:p[i].tri_idx){
            triangle curr=tri[idx];
            if(i==curr.ia){
                auto tempp=[rho,curr](point a)->double{
                    return rho(a)*(curr.p1[0]*a.x+curr.p2[0]*a.y+curr.p3[0]);
                };
                f[i]-=triangle_integral(tempp,curr);
                K[i*(N-1)+i]+=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[0]+curr.p2[0]*curr.p2[0]);
                K[i*(N-1)+curr.ib]+=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[1]+curr.p2[0]*curr.p2[1]);
                if(curr.ic!=N-1){
                    K[i*(N-1)+curr.ic]+=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[2]+curr.p2[0]*curr.p2[2]);
                }
                if(p[i].edge){
                    auto temp=[g,curr](point pt)->double{
                        return g(pt)*(curr.p1[0]*pt.x+curr.p2[0]*pt.y+curr.p3[0]);
                    };
                    if(p[curr.ib].edge){
                        f[i]+=integral(temp,p[i],p[curr.ib]);
                    }
                    if(p[curr.ic].edge){
                        f[i]+=integral(temp,p[i],p[curr.ic]);
                    }
                }
            }
            else if(i==curr.ib){
                auto tempp=[rho,curr](point a)->double{
                    return rho(a)*(curr.p1[1]*a.x+curr.p2[1]*a.y+curr.p3[1]);
                };
                f[i]-=triangle_integral(tempp,curr);
                K[i*(N-1)+i]+=fabs(curr.d)/2.0*(curr.p1[1]*curr.p1[1]+curr.p2[1]*curr.p2[1]);
                K[i*(N-1)+curr.ia]+=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[1]+curr.p2[0]*curr.p2[1]);
                if(curr.ic!=N-1){
                    K[i*(N-1)+curr.ic]+=fabs(curr.d)/2.0*(curr.p1[1]*curr.p1[2]+curr.p2[1]*curr.p2[2]);
                }
                if(p[i].edge){
                    auto temp=[g,curr](point pt)->double{
                        return g(pt)*(curr.p1[1]*pt.x+curr.p2[1]*pt.y+curr.p3[1]);
                    };
                    if(p[curr.ia].edge){
                        f[i]+=integral(temp,p[curr.ia],p[curr.ib]);
                    }
                    if(p[curr.ic].edge){
                        f[i]+=integral(temp,p[i],p[curr.ic]);
                    }
                }
            }
            else{
                auto tempp=[rho,curr](point a)->double{
                    return rho(a)*(curr.p1[2]*a.x+curr.p2[2]*a.y+curr.p3[2]);
                };
                f[i]-=triangle_integral(tempp,curr);
                K[i*(N-1)+i]+=fabs(curr.d)/2.0*(curr.p1[2]*curr.p1[2]+curr.p2[2]*curr.p2[2]);
                K[i*(N-1)+curr.ia]+=fabs(curr.d)/2.0*(curr.p1[2]*curr.p1[0]+curr.p2[2]*curr.p2[0]);
                K[i*(N-1)+curr.ib]+=fabs(curr.d)/2.0*(curr.p1[2]*curr.p1[1]+curr.p2[2]*curr.p2[1]);
                if(p[i].edge){
                    auto temp=[g,curr](point pt)->double{
                        return g(pt)*(curr.p1[2]*pt.x+curr.p2[2]*pt.y+curr.p3[2]);
                    };
                    if(p[curr.ia].edge){
                        f[i]+=integral(temp,p[curr.ia],p[curr.ic]);
                    }
                    if(p[curr.ib].edge){
                        f[i]+=integral(temp,p[curr.ib],p[curr.ic]);
                    }
                }
            }
        }
    }
    //LU分解
    my_solver(N-1,K,f);
    PHI res;
    for(int i=0;i<N-1;i++){
        //加入内部点电势
        res.phi.push_back(f[i]);
    }
    res.phi.push_back(0);
    free(K);
    free(f);
    return res;
}