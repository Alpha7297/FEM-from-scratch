#include "geometry.hpp"
#include "gen.h"
#include<iostream>
using std::vector;
vector<point> edges;
polygon* poly;
double f(double theta){
    return 1+cos(theta);
}
void test(point a){
    if(poly->in(a)){
        printf("%lf %lf in\n",a.x,a.y);
    }
    else{
        printf("%lf %lf not in\n",a.x,a.y);
    }
}
int main(void){
    edges=gen_edge(f);
    polygon polyy=polygon(edges);
    FILE* pf=fopen("points.txt","w");
    for(const auto& po:polyy.edges){
        fprintf(pf,"%lf %lf\n",po.x,po.y);
    }
    fclose(pf);
    poly=&polyy;
    vector<point> test_point={
        {0.5,1},{0.5,1.5},{2.5,1.5},{1,1},{2,1.5}
    };
    for(int i=0;i<test_point.size();i++){
        test(test_point[i]);
    }
    return 0;
}