#include "geometry.hpp"
#include<iostream>
#include<algorithm>
#include<cstdio>
#include<cstdlib>
#include<map>
#include<omp.h>
double r_avg=0.05;//平均点间距
int k_poss=30;//每个活跃点生成可能点
using std::vector;
typedef struct r_point{
    double r;
    double theta;
    bool operator<(const r_point& other){
        return theta<other.theta;
    }
}r_point;
typedef struct RES{
    vector<point> pr;
    vector<triangle> trir;
}RES;
typedef struct bucket{//按x距离分桶加速搜索，空间换时间
    double xend;
    vector<point> ps;
}bucket;
typedef struct edge{
    int id1,id2;
    bool operator<(const edge& other)const{//重载运算符用于map
        if(id1!=other.id1){
            return id1<other.id1;
        }
        return id2<other.id2;
    }
}edge;
class buckets{
public:
    vector<bucket> b;
    int size;
    buckets(vector<point> edges){
        if(edges.empty()){
            size=0;
            return;
        }
        double minx=edges[0].x,maxx=edges[0].x;
        for(int i=1;i<edges.size();i++){//找最大最小值
            if(minx>edges[i].x){
                minx=edges[i].x;
            }
            if(maxx<edges[i].x){
                maxx=edges[i].x;
            }
        }
        double dx;
        if((maxx-minx)/(2*r_avg)>1000){//防止分过多的桶直接MLE
            dx=(maxx-minx)/1000.0;
        }
        else{
            dx=2*r_avg;
        }
        int idx=0;
        size=(maxx-minx+dx-1)/dx+2;
        for(int i=0;i<size;i++){
            bucket temp;
            temp.xend=minx+i*dx;
            b.push_back(temp);
        }
        b.back().xend=100;//一个很大的数保证所有点都在范围内
        for(int i=0;i<edges.size();i++){
            int left=0,right=size-1;
            int res=right;
            while(left<=right){
                int mid=(left+right)/2;
                if(b[mid].xend>=edges[i].x){
                    res=mid;
                    right=mid-1;
                }
                else{
                    left=mid+1;
                }
            }
            b[res].ps.push_back(edges[i]);
        }
    }
};
vector<int> gen_active;
vector<point> gen_p;
vector<triangle> gen_tri;
buckets* gen_bs=nullptr;
polygon* gen_poly=nullptr;
int find(double x){//二分查找桶编号
    int left=0;
    int right=gen_bs->size-1;
    int res=right;
    while(left<=right){
        int mid=(right+left)/2;
        if(gen_bs->b[mid].xend>=x){
            res=mid;
            right=mid-1;
        }
        else{
            left=mid+1;
        }
    }
    return res;
}
void addpoint(point neww){//加入新的点
    neww.edge=0;
    gen_p.push_back(neww);
    int idx=find(neww.x);
    gen_bs->b[idx].ps.push_back(neww);
    vector<int> to_delete;//寻找外接圆包含这个点的三角形，进行修改
    std::map<edge,int> pairs;
    for(int i=0;i<gen_tri.size();i++){
        point o=gen_tri[i].o;
        if(distance(o,gen_p.back())<gen_tri[i].r){
            to_delete.push_back(i);
            edge temp={gen_tri[i].ia,gen_tri[i].ib};
            pairs[temp]++;
            temp={gen_tri[i].ia,gen_tri[i].ic};
            pairs[temp]++;
            temp={gen_tri[i].ib,gen_tri[i].ic};
            pairs[temp]++;
        }
    }
    vector<triangle> new_tri;
    int id2=0;
    for(int i=0;i<gen_tri.size();i++){
        if(id2<to_delete.size()&&to_delete[id2]==i){
            id2++;
            continue;
        }
        new_tri.push_back(gen_tri[i]);
    }
    gen_tri.swap(new_tri);
    int new_idx=gen_p.size()-1;
    for(const auto& pair:pairs){
        if(pair.second==1){
            const edge& e=pair.first;
            point a=gen_p[e.id1];
            point b=gen_p[e.id2];
            point c=gen_p[new_idx];
            gen_tri.push_back({a,b,c,e.id1,e.id2,new_idx});
        }
    }
    gen_active.push_back(new_idx);
}
int search(point curr,int id){
    vector<point> temp=gen_bs->b[id].ps;
    for(int i=0;i<temp.size();i++){
        if(distance(temp[i],curr)<r_avg){
            return 1;
        }
    } 
    return 0;
}
void delpoint(int idx){
    vector<int> to_delete;
    for(int i=0;i<gen_tri.size();i++){
        if(gen_tri[i].ia==idx||gen_tri[i].ib==idx||gen_tri[i].ic==idx){
            to_delete.push_back(i);
        }
    }
    vector<triangle> new_tri;
    int id2=0;
    for(int i=0;i<gen_tri.size();i++){
        if(id2<to_delete.size()&&to_delete[id2]==i){
            id2++;
            continue;
        }
        new_tri.push_back(gen_tri[i]);
    }
    gen_tri.swap(new_tri);
    for(int i=0;i<gen_tri.size();i++){
        if(gen_tri[i].ia>idx){
            gen_tri[i].ia--;
        }
        if(gen_tri[i].ib>idx){
            gen_tri[i].ib--;
        }
        if(gen_tri[i].ic>idx){
            gen_tri[i].ic--;
        }
    }
    gen_p.erase(gen_p.begin()+idx);
}
void delout(){//删除多边形外的点，使用重心在多边形外简单判断
    vector<int> to_delete;
    for(int i=0;i<gen_tri.size();i++){
        point a=gen_tri[i].a;
        point b=gen_tri[i].b;
        point c=gen_tri[i].c;
        point mid=(a+b+c)/3;
        if(!gen_poly->in(mid)){
            to_delete.push_back(i);
        }
    }
    vector<triangle> new_tri;
    int id2=0;
    for(int i=0;i<gen_tri.size();i++){
        if(id2<to_delete.size()&&to_delete[id2]==i){
            id2++;
            continue;
        }
        new_tri.push_back(gen_tri[i]);
    }
    gen_tri.swap(new_tri);
}
RES triangle_gen(vector<point> edges){
    for(int i=0;i<edges.size();i++){//进行微扰，防止共线或四点共圆导致退化
        edges[i].x+=(rand()%100)*1e-10;
        edges[i].y+=(rand()%100)*1e-10;
    }
    gen_p.clear();
    gen_tri.clear();
    gen_active.clear();
    buckets bss=buckets(edges);
    gen_bs=&bss;
    polygon polyy=polygon(edges);
    gen_poly=&polyy;
    gen_p.push_back({0,100});
    gen_p.push_back({-100,-100});
    gen_p.push_back({100,-100});//初始超级大三角形，之后删除
    gen_tri.push_back({gen_p[0],gen_p[1],gen_p[2],0,1,2});
    for(int i=0;i<edges.size();i++){//加入边界点
        addpoint(edges[i]);
        gen_p.back().edge=1;
    }
    while(!gen_active.empty()){
        point curr=gen_p[gen_active.back()];
        gen_active.pop_back();
        for(int i=0;i<k_poss;i++){
            double dis=(rand()%100)*0.01*r_avg+r_avg;
            double theta=(rand()%100)*0.01*2*PI;
            double x=curr.x+dis*cos(theta);
            double y=curr.y+dis*sin(theta);
            point next={x,y};
            if(!(gen_poly->in(next))){
                continue;
            }
            int idx=find(x);
            int flag=0;
            if(idx>0){
                flag=search(next,idx-1);
                if(flag)continue;
            }
            if(idx<gen_bs->size-1){
                flag=search(next,idx+1);
                if(flag)continue;
            }
            flag=search(next,idx);
            if(flag)continue;
            addpoint(next);
        }
    }
    delout();//删除多边形外三角形
    delpoint(2);
    delpoint(1);
    delpoint(0);
    return {gen_p,gen_tri};
}
vector<point> better_edge(vector<point> origin){
    gen_p.clear();
    point zero={0,0};
    int len=origin.size();
    for(int i=0;i<len;i++){
        point a=origin[i];
        point b=(i==len-1)?origin[0]:origin[i+1];
        if(distance(a,b)>2*r_avg){
            point delta=(b-a);
            point temp=delta/distance(zero,delta);
            delta=temp;
            point curr=a;
            while((curr.x>b.x)==(a.x>b.x)&&(curr.y>b.y)==(a.y>b.y)){
                gen_p.push_back(curr);
                curr=curr+delta*((rand()%100)*0.01+1)*r_avg;
            }
        }
    }
    return gen_p;
}
vector<point> gen_edge(double (*f)(double)){//输入r(theta)，使用类泊松圆盘取点
    gen_p.clear();
    gen_active.clear();
    gen_p.push_back({f(0),0});
    gen_active.push_back(0);
    while(!gen_active.empty()){
        point curr=gen_p[gen_active.back()];
        gen_active.pop_back();
        double theta0=atan4(curr);
        double rdot=(f(theta0+0.01)-f(theta0))/0.01;
        double r0=f(theta0);
        for(int i=0;i<k_poss;i++){
            double sign=(rand()%2==0)?1:-1;
            double thetan;
            if(sqrt(r0*r0+rdot*rdot)<EPS){
                thetan=theta0+sign*2*PI*(rand()%100)*0.01*0.01;//防止溢出
            }
            else{
                thetan=theta0+sign*((rand()%100)*0.01*r_avg+r_avg)/sqrt(r0*r0+rdot*rdot);//弧长在[r,2r]内
            }
            double x=f(thetan)*cos(thetan);
            double y=f(thetan)*sin(thetan);
            point next={x,y};
            int flag=0;
            for(auto& po:gen_p){
                if(distance(next,po)<r_avg){
                    flag=1;
                    break;
                }
            }
            if(!flag){
                gen_p.push_back(next);
                gen_active.push_back(gen_p.size()-1);
            }
        }
    }
    std::sort(gen_p.begin(),gen_p.end(),[](const point& a,const point& b){return atan4(a)<atan4(b);});
    return gen_p;
}