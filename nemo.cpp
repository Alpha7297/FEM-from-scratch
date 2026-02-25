#include "gen.h"
#include "cal.h"
#include<chrono>
//注意所有点到原点距离要小于100，如果想更大在gen.cpp里修改大三角形和分桶上限
using std::vector;
using namespace std::chrono;
extern int np;
extern double r_avg;
extern int k_poss;
double g(point a){//诺伊曼边界条件，E_n(x,y)=g({x,y})
    return 1/(4*PI);
}
double f(double theta){//r(theta)
    return 2;
}
double rho(point a){//\nabla^2 phi=rho
    if(fabs(a.x-1)<0.5&&fabs(a.y)<0.5){
        return 1;
    }
    return 0;
}
int main(void){
    np=16;//LU分解分块数，np=16是反复测试的最好结果，32核CPU上分解32767^2稠密矩阵速度可以达到270s
    r_avg=0.03;//点平均间隔，太小会MLE，耗时指数增长
    k_poss=30;//泊松圆盘取点，每个活跃点可生成随机点个数，设置越大速度越慢，但可能更均匀
    auto start_time=steady_clock::now();

    vector<point> edges=gen_edge(f);//第一种生成边界方式，给定r(theta)，目前不支持自交图形，如r=cosnθ的玫瑰线
    //vector<point> origin_edges={//第二种生成边界方式，直接给定边界上的点
    //    {0,0},{2,0},{2,2},{0,2}
    //};
    //vector<point> edges=better_edge(origin_edges);//由于原来点间隔可能过大，通过gen.cpp中的better_edge函数自动转换成更好的网格
    RES ans=triangle_gen(edges);//生成三角网格

    auto end_time=steady_clock::now();
    duration<double> elapse_seconds=end_time-start_time;
    printf("网格生成用时%.5lfs,生成%zu个点,%zu个三角形\n",elapse_seconds,ans.pr.size(),ans.trir.size());

    start_time=steady_clock::now();
    
    PHI phi_data=neum_cal_phi(g,rho,ans.pr,ans.trir);//求解迪利克雷边界下的电势分布

    end_time=steady_clock::now();
    elapse_seconds=end_time-start_time;
    printf("求解矩阵用时%.5lfs\n",elapse_seconds);

    FILE* pf=fopen("points.txt","w");
    for(int i=0;i<ans.pr.size();i++){
        fprintf(pf,"%lf %lf\n",ans.pr[i].x,ans.pr[i].y);
    }
    FILE* trif=fopen("triangles.txt","w");
    for(int i=0;i<ans.trir.size();i++){
        fprintf(trif,"%d %d %d\n",ans.trir[i].ia,ans.trir[i].ib,ans.trir[i].ic);
    }
    FILE* phif=fopen("phi.txt","w");
    for(int i=0;i<ans.pr.size();i++){
        fprintf(phif,"%lf\n",phi_data.phi[i]);
    }
    fclose(pf);
    fclose(trif);
    fclose(phif);
    return 0;
}