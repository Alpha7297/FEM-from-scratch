#include "geometry.hpp"
#include<cstdio>
#include<cstdlib>
#include<algorithm>
#include<functional>
#include<omp.h>
double alpha=1;
using std::vector;
//储存结果
typedef struct PHI{
    vector<double> phi;
}PHI;
vector<double> phi;
vector<double> new_phi;
double triangle_integral(std::function<double(point)> rho,triangle a){
    point p1=(a.a*1/6+a.b*1/6+a.c*2/3);
    point p2=(a.a*1/6+a.b*2/3+a.c*1/6);
    point p3=(a.a*2/3+a.b*1/6+a.c*1/6);
    double res=fabs(a.d)/6*(rho(p1)+rho(p2)+rho(p3));
    return res;
}
PHI dirc_cal_phi(double (*g)(point),double (*rho)(point),vector<point> p,vector<triangle> tri){
    phi.clear();
    new_phi.clear();
    //读取
    for(int i=0;i<tri.size();i++){
        triangle curr=tri[i];
        p[curr.ia].tri_idx.push_back(i);
        p[curr.ib].tri_idx.push_back(i);
        p[curr.ic].tri_idx.push_back(i);
    }
    int N=p.size();
    int len=0;
    for(;len<p.size();len++){
        if(!p[len].edge){
            break;
        }
    }
    for(int i=0;i<N;i++){
        phi.push_back(g(p[i]));
        new_phi.push_back(phi.back());
    }
    //并行对每一个内部点的电势求偏导，得到矩阵方程
    //E=\sum_i S_i(p_i1^2+p_i2^2)
    for(int epoch=0;;epoch++){
        double err=0;
        #pragma omp parallel for reduction(max:err)
        for(int i=len;i<N;i++){
            double K=0;
            double f=0;
            for(const auto&idx:p[i].tri_idx){
                triangle curr=tri[idx];
                if(i==curr.ia){
                    auto temp=[rho,curr](point a)->double{
                        return rho(a)*(curr.p1[0]*a.x+curr.p2[0]*a.y+curr.p3[0]);
                    };
                    f-=triangle_integral(temp,curr);
                    K+=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[0]+curr.p2[0]*curr.p2[0]);
                    f-=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[1]+curr.p2[0]*curr.p2[1])*phi[curr.ib];
                    f-=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[2]+curr.p2[0]*curr.p2[2])*phi[curr.ic];
                }
                else if(i==curr.ib){
                    auto temp=[rho,curr](point a)->double{
                        return rho(a)*(curr.p1[1]*a.x+curr.p2[1]*a.y+curr.p3[1]);
                    };
                    f-=triangle_integral(temp,curr);
                    K+=fabs(curr.d)/2.0*(curr.p1[1]*curr.p1[1]+curr.p2[1]*curr.p2[1]);
                    f-=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[1]+curr.p2[0]*curr.p2[1])*phi[curr.ia];
                    f-=fabs(curr.d)/2.0*(curr.p1[1]*curr.p1[2]+curr.p2[1]*curr.p2[2])*phi[curr.ic];
                }
                else{
                    auto temp=[rho,curr](point a)->double{
                        return rho(a)*(curr.p1[2]*a.x+curr.p2[2]*a.y+curr.p3[2]);
                    };
                    f-=triangle_integral(temp,curr);
                    K+=fabs(curr.d)/2.0*(curr.p1[2]*curr.p1[2]+curr.p2[2]*curr.p2[2]);
                    f-=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[2]+curr.p2[0]*curr.p2[2])*phi[curr.ia];
                    f-=fabs(curr.d)/2.0*(curr.p1[1]*curr.p1[2]+curr.p2[1]*curr.p2[2])*phi[curr.ib];
                }
            }
            new_phi[i]=f/K;
            if(err<fabs(phi[i]-new_phi[i])){
                err=fabs(phi[i]-new_phi[i]);
            }
        }
        #pragma omp parallel for
        for(int i=len;i<N;i++){
            phi[i]=(1-alpha)*phi[i]+alpha*new_phi[i];
        }
        if(err<EPS){
            break;
        }
    }
    PHI res;
    for(int i=0;i<N;i++){
        res.phi.push_back(phi[i]);
    }
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
    phi.clear();
    new_phi.clear();
    for(int i=0;i<tri.size();i++){
        triangle curr=tri[i];
        p[curr.ia].tri_idx.push_back(i);
        p[curr.ib].tri_idx.push_back(i);
        p[curr.ic].tri_idx.push_back(i);
    }
    int N=p.size();
    for(int i=0;i<N;i++){
        phi.push_back(0);
        new_phi.push_back(0);
    }
    //并行对每一个内部点的电势求偏导，得到矩阵方程
    //固定最后一个节点电势为0
    for(int epoch=0;;epoch++){
        double err=0;
        #pragma omp parallel for reduction(max:err)
        for(int i=0;i<N-1;i++){
            double K=0,f=0;
            for(const auto&idx:p[i].tri_idx){
                triangle curr=tri[idx];
                if(i==curr.ia){
                    auto tempp=[rho,curr](point a)->double{
                        return rho(a)*(curr.p1[0]*a.x+curr.p2[0]*a.y+curr.p3[0]);
                    };
                    f-=triangle_integral(tempp,curr);
                    K+=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[0]+curr.p2[0]*curr.p2[0]);
                    f-=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[1]+curr.p2[0]*curr.p2[1])*phi[curr.ib];
                    f-=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[2]+curr.p2[0]*curr.p2[2])*phi[curr.ic];
                    if(p[i].edge){
                        auto temp=[g,curr](point pt)->double{
                            return g(pt)*(curr.p1[0]*pt.x+curr.p2[0]*pt.y+curr.p3[0]);
                        };
                        if(p[curr.ib].edge){
                            f+=integral(temp,p[i],p[curr.ib]);
                        }
                        if(p[curr.ic].edge){
                            f+=integral(temp,p[i],p[curr.ic]);
                        }
                    }
                }
                else if(i==curr.ib){
                    auto tempp=[rho,curr](point a)->double{
                        return rho(a)*(curr.p1[1]*a.x+curr.p2[1]*a.y+curr.p3[1]);
                    };
                    f-=triangle_integral(tempp,curr);
                    K+=fabs(curr.d)/2.0*(curr.p1[1]*curr.p1[1]+curr.p2[1]*curr.p2[1]);
                    f-=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[1]+curr.p2[0]*curr.p2[1])*phi[curr.ia];
                    f-=fabs(curr.d)/2.0*(curr.p1[1]*curr.p1[2]+curr.p2[1]*curr.p2[2])*phi[curr.ic];
                    if(p[i].edge){
                        auto temp=[g,curr](point pt)->double{
                            return g(pt)*(curr.p1[1]*pt.x+curr.p2[1]*pt.y+curr.p3[1]);
                        };
                        if(p[curr.ia].edge){
                            f+=integral(temp,p[i],p[curr.ia]);
                        }
                        if(p[curr.ic].edge){
                            f+=integral(temp,p[i],p[curr.ic]);
                        }
                    }
                }
                else{
                    auto tempp=[rho,curr](point a)->double{
                        return rho(a)*(curr.p1[2]*a.x+curr.p2[2]*a.y+curr.p3[2]);
                    };
                    f-=triangle_integral(tempp,curr);
                    K+=fabs(curr.d)/2.0*(curr.p1[2]*curr.p1[2]+curr.p2[2]*curr.p2[2]);
                    f-=fabs(curr.d)/2.0*(curr.p1[0]*curr.p1[2]+curr.p2[0]*curr.p2[2])*phi[curr.ia];
                    f-=fabs(curr.d)/2.0*(curr.p1[1]*curr.p1[2]+curr.p2[1]*curr.p2[2])*phi[curr.ib];
                    if(p[i].edge){
                        auto temp=[g,curr](point pt)->double{
                            return g(pt)*(curr.p1[2]*pt.x+curr.p2[2]*pt.y+curr.p3[2]);
                        };
                        if(p[curr.ia].edge){
                            f+=integral(temp,p[i],p[curr.ia]);
                        }
                        if(p[curr.ib].edge){
                            f+=integral(temp,p[i],p[curr.ib]);
                        }
                    }
                }
            }
            new_phi[i]=f/K;
            if(err<fabs(phi[i]-new_phi[i])){
                err=fabs(phi[i]-new_phi[i]);
            }
        }
        #pragma omp parallel for
        for(int i=0;i<N-1;i++){
            phi[i]=(1-alpha)*phi[i]+alpha*new_phi[i];
        }
        if(err<1e-5){
            break;
        }
    }
    PHI res;
    for(int i=0;i<N;i++){
        res.phi.push_back(phi[i]);
    }
    return res;
}