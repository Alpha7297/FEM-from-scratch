#ifndef GEN_H
#define GEN_H
#include "geometry.hpp"
#include<iostream>
#include<algorithm>
#include<cstdio>
#include<cstdlib>
#include<map>
extern double r_avg;
extern int k_poss;
typedef struct RES{
    std::vector<point> pr;
    std::vector<triangle> trir;
}RES;
RES triangle_gen(std::vector<point> edges);
std::vector<point> gen_edge(double (*f)(double));
std::vector<point> better_edge(std::vector<point> origin);
#endif