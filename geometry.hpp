#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP
#define EPS 1e-10
#define PI 3.141592653589
#include<cmath>
#include<vector>
class point{//point类
public:
    double x,y;
    int edge;//是否是边界上的点
    std::vector<int> tri_idx;//三角形邻边
    point(double xx=0,double yy=0):x(xx),y(yy){}
    double dist(point& other)const{
        return sqrt(pow(x-other.x,2)+pow(y-other.y,2));
    }
    point operator+(const point& other){
        return {x+other.x,y+other.y};
    }
    point operator-(const point& other){
        return {x-other.x,y-other.y};
    }
    point operator/(const double& other){
        return {x/other,y/other};
    }
    point operator*(const double& other){
        return {x*other,y*other};
    }
};
class polygon{
public:
    std::vector<point> edges;//顺时针或逆时针
    polygon(std::vector<point> edgess){
        edges=edgess;
    }
    int in(const point& other){//判断一个点在不在多边形内，射线法
        int count=0;
        for (int i=0;i<edges.size();i++) {
            point a=edges[i];
            point b=edges[(i+1)%edges.size()];
            if(fabs((b.y-a.y)*(other.x-a.x)-(b.x-a.x)*(other.y-a.y))<EPS){
                if(other.x>=fmin(a.x,b.x)-EPS&&other.x<=fmax(a.x,b.x)+EPS){
                    if(other.y>=fmin(a.y,b.y)-EPS&&other.y<=fmax(a.y,b.y)+EPS){
                        return 1;
                    }
                }
            }
            if(((a.y>other.y)!=(b.y>other.y))){
                double x_inter=a.x+(other.y-a.y)*(b.x-a.x)/(b.y-a.y);
                if(x_inter>other.x+EPS)count++;
            }
        }
        return count%2;
    }
};
class triangle{//三角形类
public:
    point a,b,c;
    int ia,ib,ic;
    double x1,x2,x3,y1,y2,y3;
    double p1[3],p2[3];//phi=(p1 /cdot phi)*x+(p2 /cdot phi)*y+p3
    double p3[3];
    double d;
    point o;//外接圆圆心
    double r;//外接圆半径
    triangle(point aa,point bb,point cc,int iaa=0,int ibb=0,int icc=0):a(aa),b(bb),c(cc),ia(iaa),ib(ibb),ic(icc){
        x1=a.x,x2=b.x,x3=c.x,y1=a.y,y2=b.y,y3=c.y;
        d=x1*(y2-y3)+y1*(x3-x2)+x2*y3-x3*y2;
        p1[0]=(y2-y3)/d,p1[1]=(y3-y1)/d,p1[2]=(y1-y2)/d;
        p2[0]=(x3-x2)/d,p2[1]=(x1-x3)/d,p2[2]=(x2-x1)/d;
        p3[0]=(x2*y3-x3*y2)/d,p3[1]=(x3*y1-y3*x1)/d,p3[2]=(x1*y2-x2*y1)/d;
        double temp=(a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y))*2;
        if(fabs(temp)<EPS){//共线，通过微扰避免了这种情况
            o.x=0;
            o.y=0;
            r=0;
        }
        double a1=a.x*a.x+a.y*a.y;
        double a2=b.x*b.x+b.y*b.y;
        double a3=c.x*c.x+c.y*c.y;
        o.x=(a1*(b.y-c.y)+a2*(c.y-a.y)+a3*(a.y-b.y))/temp;
        o.y=(a1*(c.x-b.x)+a2*(a.x-c.x)+a3*(b.x-a.x))/temp;
        r=sqrt((o.x-a.x)*(o.x-a.x)+(o.y-a.y)*(o.y-a.y));
    }
};
inline double distance(point& a,point& b){
    return sqrt(pow(a.x-b.x,2)+pow(a.y-b.y,2));
}
inline double atan4(point a){
    if(fabs(a.x)<EPS){
        if(a.y>0){
            return PI/2.0;
        }
        else{
            return 1.5*PI;
        }
    }
    if(fabs(a.y)<EPS){
        if(a.x>0){
            return 0;
        }
        else{
            return PI;
        }
    }
    if(a.x>0){
        if(a.y>0){
            return atan(a.y/a.x);
        }
        else{
            return 2*PI-atan(-a.y/a.x);
        }
    }
    else{
        return PI+atan(a.y/a.x);
    }
}
#endif