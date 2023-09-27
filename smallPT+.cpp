/*路径追踪的基本思想是从视点发出一条光线,光线与物体表面相交时根据表面的材质属性继续采样一个方向,发出另一条光线,如此
迭代,直到光线打到光源上(或逃逸出场景),然后用蒙特卡洛的方法,计算其贡献,作为像素的颜色值.
路径追踪会避开贡献小的路径,而在贡献大的路径附近做更多局部的探索.
路径追踪=光线追踪+蒙特卡洛方法*/
#define _CRT_SECURE_NO_WARNINGS//禁用fopen的警告
#define _SILENCE_AMP_DEPRECATION_WARNINGS
#define seeds_ 490
#include "stdafx.h"
#include <math.h>//smallpt,一个小型路径追踪器
#include <iostream>
#include <stdlib.h> 
#include <random>
#include <omp.h>
#include <Windows.h>
#include <GLFW/glfw3.h>
#include <amp.h>
#include <amp_math.h>
using namespace concurrency;
using namespace concurrency::precise_math;
using namespace std;
//#define double float//使用单精度浮点数以节省内存
#define float2 double
#define M_PI 3.1415926525//使用常量M_PI表示圆周率
//#define IMAGE_PATHNAME "image.ppm"
double erand48(unsigned short xsubi[3]) {//随机数
    return (double)rand() / (double)RAND_MAX;
}
double erand49(long seed) __GPU {
    // 设置参数
    long a = 1664525;
    long c = 1013904223;
    long m = (929496729);
    // 生成随机数
    long seed1 = (a * seed + c) % m;
    double random = static_cast<double>(seed1) / m;
    random < 0 ? random *= -1 : random = random;
    return random;
}
double erand47(long seed) __GPU {

    int a = 127;  // 随意选取的常数，可以自行修改
    int b = 4096; // 随意选取的常数，可以自行修改
    int p = 2001; // 选取一个大于2000的质数，可以自行修改

    // 使用算法生成随机数
    seed = (a * seed + b) % p;
    return seed;
}
/*    // 设置参数
    long a = 9753625;
    long c = 211309423;
    long m = (99940672);
    // 生成随机数
    long seed1 = (a * seed + c) % m;
    double random = static_cast<double>(seed1) / m;
    random < 0 ? random *= -1 : random = random; 
    return random;*/
struct Vec {//三维向量
    double x, y, z;//位置或颜色都能使用
    Vec(double x_ = 0, double y_ = 0, double z_ = 0) __GPU { x = x_; y = y_; z = z_; }//构造函数,x,y,z默认为零
    Vec operator-() const __GPU
    {
        return Vec(-x, -y, -z);
    }
    Vec operator+=(const Vec& a) __GPU
    {
        return Vec(x + a.x, y + a.y, z + a.z);
    }

    Vec operator-=(const Vec& a) __GPU
    {
        return Vec(x - a.x, y - a.y, z - a.z);
    }
    Vec operator+(const Vec& b) const __GPU { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec operator-(const Vec& b) const __GPU { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec operator*=(const double& a) __GPU
    {
        return Vec(x * a, y * a, z * a);
    }
    Vec operator/=(const double& a) __GPU
    {
        return Vec(x / a, y / a, z / a);
    }
    Vec operator*(double b) const __GPU { return Vec(x * b, y * b, z * b); }
    Vec divide(const Vec& a) __GPU
    {
        return Vec(x / a.x, y / a.y, z / a.z);
    };
    Vec mult(const Vec& a) __GPU
    {
        return Vec(x * a.x, y * a.y, z * a.z);
    };
    Vec operator/(double b) const __GPU { return Vec(x / b, y / b, z / b); }
    Vec operator%(Vec& b)__GPU { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
    double dot(const Vec& b) const __GPU { return x * b.x + y * b.y + z * b.z; }
    double length() __GPU {
        return sqrt(x * x + y * y + z * z);
    }
    Vec norm()__GPU
    {
        double len = length();
        if (len > 0)
        {
            x /= len;
            y /= len;
            z /= len;
        }
        return Vec(x, y, z);
    }
    void print() {
        std::cout << "<" << x << " " << y << " " << z << ">" << '\n';
    }
};
inline Vec VecRotate(Vec vec, double angle, Vec axis)__GPU {
    double c = cos(angle);
    double s = sin(angle);
    double t = 1 - c;
    double x = vec.x, y = vec.y, z = vec.z;
    double newX = (t * axis.x * axis.x + c) * x + (t * axis.x * axis.y - s * axis.z) * y + (t * axis.x * axis.z + s * axis.y) * z;
    double newY = (t * axis.x * axis.y + s * axis.z) * x + (t * axis.y * axis.y + c) * y + (t * axis.y * axis.z - s * axis.x) * z;
    double newZ = (t * axis.x * axis.z - s * axis.y) * x + (t * axis.y * axis.z + s * axis.x) * y + (t * axis.z * axis.z + c) * z;
    return Vec(newX, newY, newZ);
}
struct Ray { Vec o, d; Ray(Vec o_, Vec d_)__GPU : o(o_), d(d_) {} };//光线结构体
enum obj_type { S, P, T,BVH };//材质类型
struct Material {
    double ref, diff, spec, refr, refr_nt;//和光线交互的概率,漫反射概率,镜面反射概率,折射概率,折射率
    Vec e, c, i;//自发光颜色,本体颜色,杂质
    bool sssEnable;//次表面散射开关
    double sss, Distortion = 0, Power = 0, Scale = 0, WrapValue = 0;//发生概率,散射扭曲值,散射强度指数,散射强度缩放系数,包裹值(平滑过渡)
    bool sb;//穿透开关
    double sb1, sb2;//穿透发生概率,穿透光线的扭曲值
    bool MFR;//微表面反射开关
    double MFR1, MFR2;//微表面反射发生概率,扭曲值
    Material(double ref_=0, double diff_=0, double spec_=0, double refr_=0, double refr_nt_=0,Vec e_=Vec(), Vec c_ = Vec(1, 1, 1), Vec i_ = Vec(0, 0, 0),
        bool sssEnable_ = false, double sss_ = 0,
        double Distortion_=0, double Power_=0, double Scale_=0,double WrapValue_=0,bool sb_=false, double sb1_=0, double sb2_=0,
        bool MFR_=false, double MFR1_=0, double MFR2_=0)__GPU:
    ref(ref_),diff(diff_),spec(spec_),refr(refr_),refr_nt(refr_nt_),e(e_),c(c_),i(i_),sssEnable(sssEnable_),sss(sss_),Distortion(Distortion_),Power(Power_),
    Scale(Scale_),WrapValue(WrapValue_),sb(sb_),sb1(sb1_),sb2(sb2_),MFR(MFR_),MFR1(MFR1_),MFR2(MFR2_){}
};
struct Obj {
    double rad ,w, h, size;//半径,反射,漫反射,镜面反射,折射概率,折射率,宽,高
    Vec p,n, p1, p2, p3;//位置,自发光,颜色(这里的颜色并不是0-255,而是0-1),杂质,法线
    Material Material_;
    obj_type type;// 物体类型
    Obj(double rad_=0, Vec p_=Vec(), obj_type refl_=S, double w_=0, double h_=0, Vec n_=Vec(), double size_=0, Vec p1_=Vec(),
        Vec p2_=Vec(), Vec p3_=Vec(), Material Material__ = Material()
    )__GPU ://构造函数
        rad(rad_), p(p_), type(refl_), w(w_), h(h_), n(n_)
        , size(size_), p1(p1_), p2(p2_), p3(p3_),Material_(Material__) {}
    double intersect_sphere(const Ray& r) const  __GPU { //计算射线原点与球体之间的交点的距离,如果没有交点返回0
        Vec op = p - r.o; //光源指向球心的一条向量
        double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
        //eps是一个很小的量,代指0
        //t是射线与交点之间的距离,是方程t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0的解"."表示点乘
        //b是光源指向球心的向量和光源方向向量的夹角的余弦值
        //det没有解(<0)则没有相交^^^|op-t|=|r|         (op-t)^2-r^2=0
        if (det < 0) return 0; else det = sqrt(det);
        if (b - det < 0.00001) return 0;
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
        //返回t值,如果为0,则表示不相交
        //选择其中较小且大于0的解,如果det<0则无解,表示不相交,返回0
    };
    double intersect_plane(const Ray& r) const __GPU {
        double t = n.dot(p - r.o) / n.dot(r.d);
        double dn = n.dot(r.d);
        int a = 1;
        if (dn >= 0) a = -1;
        // 如果光线平行于或位于平面后面,则没有相交
        if (fabs(dn) < 1e-6 || dn >= 0) {
            return 0;
        }
        if (t < 0.00001) return 0;
        // 交点
        Vec intersection = r.o + (r.d * t);
        // 检查交点是否在平面范围内
        double half_width = w / 2.0;
        double half_height = h / 2.0;
        double z_distance = intersection.z - p.z; // 计算交点到平面中心的Z轴距离

        if (intersection.x < p.x - half_width || intersection.x > p.x + half_width ||
            intersection.y < p.y - half_height || intersection.y > p.y + half_height ||
            z_distance < 0 || z_distance > h) { // 检查交点在Z轴方向上是否在范围内
            return 0;
        }

        // 返回相交点与光线原点之间的距离
        return t;
    };
    double intersect_triangle(const Ray& r) const __GPU { //计算射线原点与球体之间的交点的距离,如果没有交点返回0
        /*Vec p2_ = p2 - p1;
        Vec p3_ = p3 - p1;
        Vec p1 = p;
        Vec p2 = p1 + p2_;
        Vec p3 = p1 + p3_;
        Vec nl = n;
        if (n.dot(r.d) > 0) nl=n * -1;
        if (fabs(nl.dot(r.d)) < -1e4) return 0;
        double t = (nl.dot(p1) - r.o.dot(nl)) / r.d.dot(nl);
        if (t < 0.0005f) return 0;
        Vec P = r.o +r.d*t;
        Vec a, b;
        a = p2 - p1;
        b = P - p1;
        Vec c1 = a%b;
        a = p3 - p2;
        b = P - p2;
        Vec c2 = a % b;
        a = p1 - p3;
        b = P - p3;
        Vec c3 = a % b;
        if (c1.dot(n) < 0 || c2.dot(n) < 0 || c3.dot(n) < 0) return 0;
        return t;*/
        const double EPSILON = 0.0000001; //要比较的小值
        Vec p1_ = p1 * size + p;//根据三角形位置计算三个顶点的坐标
        Vec p2_ = p2 * size + p;//根据三角形位置计算三个顶点的坐标
        Vec p3_ = p3 * size + p;//根据三角形位置计算三个顶点的坐标
        Vec p2 = p + p2_;//新的三点坐标
        Vec p3 = p + p3_;//新的三点坐标
        Vec p1 = p + p1_;//新的三点坐标
        Vec edge1 = p2 - p1;
        Vec edge2 = p3 - p1;
        Vec rd = r.d;
        Vec h = rd % edge2;
        double a = edge1.dot(h);
        Vec AB = p2 - p1, AC = p3 - p1, n0 = AB % AC;
        Vec n = n0.norm();
        double tdn = n.dot(r.d);
        if (a > -EPSILON && a < EPSILON || tdn >= 0)
            return 0.0; //射线平行于三角形

        double f = 1.0 / a;
        Vec s = r.o - p1;
        double u = f * s.dot(h);

        if (u < 0.0 || u > 1.0)
            return 0.0; //交点在三角形之外

        Vec q = s % edge1;
        double v = f * r.d.dot(q);

        if (v < 0.0 || u + v > 1.0)
            return 0.0; //交点在三角形之外

        double t = f * edge2.dot(q);

        if (t > EPSILON)
            return t; //找到交点

        return 0.0; //未找到交点
    };
};
struct Camera {
    Vec o, d;          // Camera position and direction
    Vec right, up;         // Basis vectors of the camera coordinate system
    double fov;            // Field of view

    Camera(Vec pos_ = Vec(0, 0, 0), Vec dir_ = Vec(0, 0, -1), double fov_ = 45.0)__GPU
        : o(pos_), d(dir_.norm()), fov(fov_) {
        right = Vec(1, 0, 0);
        up = Vec(0, 1, 0);
    }

    Ray getRay(double x, double y, double width, double height) const __GPU{
        Vec imagePlane = d + right * ((2 * x / width - 1) * tan(fov / 2 * M_PI / 180) * width / height) + up * (1 - 2 * y / height) * tan(fov / 2 * M_PI / 180);
        return Ray(o, imagePlane.norm());
    }
};
//这两个函数是亮度和颜色计算的辅助函数
inline double clamp(double x) __GPU { return x < 0 ? 0 : x>1 ? 1 : x; }//对于经过递归叠加后大于1的设为1,小于0的设为0
/*半径,反射概率,漫反射概率,镜面反射概率,折射概率,折射率,杂质,位置,自发光,颜色,类型,
宽,高,法线,大小,三点位置*/
Obj scenes1[] = {
  Obj(13,Vec(50,0,-500), P,
  1000,1000,Vec(0,1,0)
  ,1,Vec(1,1,1),Vec(2,2,2),Vec(3,3,3),Material(100,50,50,0,0,0)),//
    Obj(150,Vec(250,400,0), S,
  1000,1000,Vec(0,1,0)
  ,1,Vec(1,1,1),Vec(2,2,2),Vec(3,3,3),Material(100,100,0,0,0,Vec(12,12,12))),//
   Obj(30,Vec(50,80,0), S,
  1000,1000,Vec(0,1,0)
  ,1,Vec(1,1,1),Vec(2,2,2),Vec(3,3,3),Material(100,0,0,0,0,Vec(),Vec(1,.45,0),Vec(),false,0,0,0,0,0,0,0,0,1,100,.9)),//
   Obj(60,Vec(170,80,0), S,
  1000,1000,Vec(0,1,0)
  ,1,Vec(1,1,1),Vec(2,2,2),Vec(3,3,3),Material(100,0,100,0,0,Vec(),Vec(0,0,1))),//
   Obj(60,Vec(300,80,0), S,
  1000,1000,Vec(0,1,0)
  ,1,Vec(1,1,1),Vec(2,2,2),Vec(3,3,3),Material(100,0,0,100,1.5,Vec(),Vec(1,1,1))),//
  Obj(60,Vec(25,35,-125), T,
  1000,1000,Vec(0,1,0)
  ,1,Vec(0,-20,50),Vec(40,20,50),Vec(-40,20,50),Material(100,100,0,0,0,Vec(),Vec(1,0,1))),//
};
//inline double clamp(double x) { return x < 0 ? 0 : x>1 ? 1 : x; }//对于经过递归叠加后大于1的设为1,小于0的设为0
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }/*把0 - 1转换为rgb中的0 - 255, 设置了各1 / 2.2的调整值,
能让画面更亮*/
Vec debug;
inline bool intersect(const Ray& r_,double& t,double& num,array_view<Obj, 1>&scenes,int &id,double &d1)__GPU {
    double inf = t = 1e20;
    for (int i = int(num); i--;) {//求出最近交点
        if (scenes[i].type == S) {
            if ((d1 = scenes[i].intersect_sphere(r_)) && d1 < t) {
                t = d1;
                id = i;
            }
        }
        else if (scenes[i].type == P) {
            if ((d1 = scenes[i].intersect_plane(r_)) && d1 < t) {
                t = d1;
                id = i;
            }
        }
        else if (scenes[i].type == T) {
            if ((d1 = scenes[i].intersect_triangle(r_)) && d1 < t) {
                t = d1;
                id = i;
            }
        }
    }
    return t < inf;
}
Vec radiance(const Ray& r, int depth, array_view<double,1> Xi, Camera cam, array_view<Obj,1>scenes, int slength,Vec am) __GPU {
    int ri = 0;
    double ncg = 1;//空气折射率
    double t;// 相交距离
    int id = 0;// 相交对象的ID
    double num = slength, d1, inf = t = 1e20;
    Ray r_ = r;
    int depth_ = depth;
    Vec cl(0, 0, 0);   // 累积颜色
    Vec cf(1, 1, 1);  // 累积反射率
    while (true) {
        if (!intersect(r_, t, num, scenes, id, d1)) return cl;
        const Obj& obj = scenes[id];//被击中的对象
        Vec x = r_.o + r_.d * t, f = obj.Material_.c;
        Vec n, pn = obj.n;
        Vec AB = obj.p2 - obj.p1, AC = obj.p3 - obj.p1, n0_t = AB % AC;
        if (obj.type == S) {
            n = (x - obj.p).norm();
        }
        else if (obj.type == P) {
            n = (pn).norm();
        }
        else if (obj.type == T) {
            n = n0_t.norm();
        }
        Vec nl = n.dot(r_.d) < 0 ? n : n * -1;
        double fanshe = erand49(Xi[2 + ri]) * 100;
        Xi[2] = erand47(Xi[2 + ri]) * 2000;
        double fanshe2 = erand49(Xi[3 + ri]) * 100;
        Xi[3] = erand47(Xi[3 + ri]) * 2000;
        double zheshe = erand49(Xi[4 + ri]) * 100;
        Xi[4] = erand47(Xi[4 + ri]) * 2000;
        /*x为交点,n为球体法向量,nl用于修正法向量,如果球体法向量和光线方向向量的点积小于零,则法线变为相反方向
        (此时光线从内部发出),(从外部来的光线,法线向外;从内部来的光线,法线向内),f为球体颜色*/
        if (depth_ > 6) return cl;//当递归深度大于6,返回黑色
        double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; /*获取RGB三个值里的最高值*/
        cl = cl + cf.mult(obj.Material_.e)+am;
        if (++depth_ > 5) {//递归达到5时有机会返回
            if (erand49(Xi[5 + ri]) < p) {
                Xi[5] = erand47(Xi[5 + ri]) * 2000;
                f = f * (1 / p);
            }
            else {
                return cl;
            }
        }
        cf = cf.mult(f);
        if (fanshe <= obj.Material_.ref) {
            if (obj.Material_.diff >= fanshe2) {// 漫反射(在半球当中随即找一个方向,然后进行递归)
                double r1 = 2 * M_PI * erand49(Xi[0 + ri]), r2 = erand49(Xi[1 + ri]), r2s = sqrt(r2);
                Xi[0] = erand47(Xi[0 + ri]) * 2000;
                Xi[1] = erand47(Xi[1 + ri]) * 2023;
                //r1为随机选取的角度,范围是 0 到 2π 之间,r2是随机选择了一个距离(0-1),r2s是距离开方的结果
                Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1, 0) : Vec(1, 0, 0)) % w).norm(), v = w % u;
                //fabs()求浮点数绝对值
                /*根据法线构造了一个正交基, w与法线同向(在n维空间当中,由n个互相正交(垂直)的向量组成一个正交基)
                在w.x的绝对值>0.1的时候,u垂直于(0,1,0)和w的单位向量,否则是垂直于(1,0,0)和w的单位向量
                这样做的目的是当w.x等于或接近0时,可能会出现线性相关的情况*/
                Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();//反射光方向向量
                //return obj.Material_.e + f.mult(radiance(Ray(x, d), depth_, Xi, cam));
                r_ = Ray(x, d);
                continue;
            }/*全局光照方程,使用蒙特卡罗方法求解*/
            /*镜面反射,
             公式为 入射光线方向向量-2*法向量*入射光方向向量dot法向量=反射光线方向向量*/
            else if (obj.Material_.spec >= fanshe2) {
                Vec spec_d = r_.d - n * 2 * n.dot(r_.d);
                r_ = Ray(x, spec_d);
                continue;
            }
            else if (obj.Material_.sb1 >= fanshe2 && obj.Material_.sb) {
                Vec d = r_.d;
                d.x += erand49(Xi[6 + ri]) < 0.5 ? erand49(Xi[7 + ri]) * -1 * obj.Material_.sb2 : erand49(Xi[8 + ri]) * obj.Material_.sb2;
                d.z += erand49(Xi[9 + ri]) < 0.5 ? erand49(Xi[10 + ri]) * -1 * obj.Material_.sb2 : erand49(Xi[11 + ri]) * obj.Material_.sb2;
                d.y += erand49(Xi[12 + ri]) < 0.5 ? erand49(Xi[13 + ri]) * -1 * obj.Material_.sb2 : erand49(Xi[14 + ri]) * obj.Material_.sb2;
                Xi[9] = erand47(Xi[6 + ri]) * 2000; Xi[7] = erand47(Xi[7 + ri]) * 2000; Xi[8] = erand47(Xi[8 + ri]) * 2000;
                Xi[6] = erand47(Xi[9 + ri]) * 2000; Xi[10] = erand47(Xi[10 + ri]) * 2000; Xi[11] = erand47(Xi[11 + ri]) * 2000;
                Xi[12] = erand47(Xi[12 + ri]) * 2000; Xi[13] = erand47(Xi[13 + ri]) * 2000; Xi[14] = erand47(Xi[14 + ri]) * 2000;
                Vec spec_t = r_.d - n * 2 * n.dot(r_.d);
                //if ((Xi) < 0.1) return obj.e + f.mult(radiance(Ray(x, spec_t), depth, Xi));
                //d.x += (Xi) * obj.sb2;
                //d.z += (Xi) * obj.sb2;
                //d.y += (Xi) * obj.sb2;
                //return obj.Material_.e + f.mult(radiance(Ray(x, d.norm()), depth_, Xi, cam));
                r_ = Ray(x, d.norm());
                continue;
            }
            else if (obj.Material_.MFR1 >= fanshe2 && obj.Material_.MFR) {
                Vec d = n;
                d.x += erand49(Xi[15 + ri]) < 0.5 ? erand49(Xi[15 + ri]) * -1 * obj.Material_.MFR2 : erand49(Xi[16 + ri]) * obj.Material_.MFR2;
                d.z += erand49(Xi[17 + ri]) < 0.5 ? erand49(Xi[18 + ri]) * -1 * obj.Material_.MFR2 : erand49(Xi[19 + ri]) * obj.Material_.MFR2;
                d.y += erand49(Xi[20 + ri]) < 0.5 ? erand49(Xi[21 + ri]) * -1 * obj.Material_.MFR2 : erand49(Xi[22 + ri]) * obj.Material_.MFR2;
                Xi[15] = erand47(Xi[15 + ri]) * 2000; Xi[15] = erand47(Xi[15 + ri]) * 2000; Xi[16] = erand47(Xi[16 + ri]) * 2000;
                Xi[17] = erand47(Xi[17 + ri]) * 2000; Xi[18] = erand47(Xi[18 + ri]) * 2000; Xi[19] = erand47(Xi[19 + ri]) * 2000;
                Xi[20] = erand47(Xi[20 + ri]) * 2000; Xi[21] = erand47(Xi[21 + ri]) * 2000; Xi[22] = erand47(Xi[22 + ri]) * 2000;
                Vec spec_t = r_.d - d * 2 * d.dot(r_.d);
                //return obj.Material_.e + f.mult(radiance(Ray(x, spec_t.norm()), depth_, Xi, cam));
                r_ = Ray(x, spec_t.norm());
                continue;
            }
            //以下为折射
            Ray reflRay(x, r_.d - n * 2 * n.dot(r_.d));// 反射光线
            bool into = n.dot(nl) > 0; /* 入射的光线是否从外面进来?如果n和nl同向则光线从外边进来,否则光线从内部发出*/
            double nc = ncg, nt = obj.Material_.refr_nt, nnt = into ? nc / nt : nt / nc, ddn = r_.d.dot(nl), cos2t;/*nc为空气的折射率,nt为玻璃的折射
            率,nnt是原介质和目标介质的比值,ddn是光线方向向量和法线nl的夹角余弦值,cos2t是cos(t)^2*/
            double r1 = 2 * M_PI * erand49(Xi[23 + ri]), r2 = erand49(Xi[24 + ri]), r2s = sqrt(r2);
            Xi[21] = erand47(Xi[23 + ri]) * 2000; Xi[22] = erand47(Xi[24 + ri]) * 2000;
            Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1, 0) : Vec(1, 0, 0)) % w).norm(), v = w % u;
            Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
            if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) {   /*全反射,当光线从较高 折射率 的 介质 进入到较低折射率的介质
                时,如果入射角大于某一临界角θc(光线远离 法线 )时,折射角将变得足够大,折射后的光线将不会离开介质*/
                //return obj.Material_.e + f.mult(radiance(reflRay, depth_, Xi, cam));
                r_ = reflRay;
                continue;
            }
            else if (obj.Material_.refr >= fanshe2) {
                //return obj.Material_.e + f.mult(radiance(Ray(x, tdir), depth_, Xi, cam)) + obj.Material_.i;
                cf = cf + obj.Material_.i;
                Vec tdir = (r_.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();//折射光线的方向
                r_ = Ray(x, tdir);
            }
            else if (zheshe <= 50) {
                //return obj.Material_.e + radiance(reflRay, depth_, Xi, cam) * (obj.Material_.spec) * 0.9;
                cf = cf * (obj.Material_.spec) * 0.9;
                r_ = reflRay;
            }
            else {
                //return obj.Material_.e + f.mult(radiance(Ray(x, d), depth_, Xi, cam)) * (obj.Material_.diff) * 0.9;
                cf = cf * (obj.Material_.diff) * 0.9;
                r_ = Ray(x, d);
            }
            continue;
        }
        else {
            return cl;
        }
    }
}
#define jz 0.01
Vec DeNoisy(int x, int y, int i, int w, int h, int samps, array_view<Vec,1>c, int i1=0) __GPU{
    Vec yansehuancun;
    double p = c[i].x > c[i].y && c[i].x > c[i].z ? c[i].x : c[i].y > c[i].z ? c[i].y : c[i].z;
    double a = 0;
    for (int j = 1; j <= i1; j++) {
        if (c[i].x <= jz && c[i].y <= jz && c[i].z <= jz && x > 0 && y > 0 && samps <= 100) {//降噪
            if (y > 0 + j && i - w - j <= w * h && i - w - j >= 0 && !(&(c[i - w - j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i - w - j];
                a++;
            }
            if (y < h - 1 - j && i + w + j <= w * h && !(&(c[i + w + j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i + w + j];
                a++;
            }
            if (x > 0 + j && i - 1 - j <= w * h && !(&(c[i - 1 - j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i - 1 - j];
                a++;
            }
            if (x < w - (1 + j) && i + 1 + j <= w * h && !(&(c[i + 1 + j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i + 1 + j];
                //a++;
            }
        }
    }
    if (yansehuancun.x > 0 && yansehuancun.y > 0 && yansehuancun.z > 0)
        return c[i] = yansehuancun * (((1 - p) / 2.5 + p) / a * 1.5);
}
Vec DeNoisy2(int x, int y, int i, int w, int h, int samps, Vec c[], int i1=0) {
    Vec yansehuancun;
    double p = c[i].x > c[i].y && c[i].x > c[i].z ? c[i].x : c[i].y > c[i].z ? c[i].y : c[i].z;
    double a = 0;
    for (int j = 1; j <= i1; j++) {
        if (c[i].x <= jz && c[i].y <= jz && c[i].z <= jz && x > 0 && y > 0 && samps <= 64) {//降噪
            if (y > 0 + j + 1 && i - w - j - 1 <= w * h && i - w - j >= 0 && !(&(c[i - w - j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i - w - j - 1];
                a++;
            }
            if (y < h - 1 - j - 1 && i + w + j + 1 <= w * h && !(&(c[i + w + j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i + w + j + 1];
                a++;
            }
            if (x > 0 + j + w && i - 1 - j - w <= w * h && !(&(c[i - 1 - j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i - 1 - j - w];
                a++;
            }
            if (x < w - (1 + j - 1) && i + 1 + j + w <= w * h && !(&(c[i + 1 + j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i + 1 + j + w];
                //a++;
            }
        }
    }
    if (yansehuancun.x > 0 && yansehuancun.y > 0 && yansehuancun.z > 0)
        return c[i] = yansehuancun * (((1 - p) / 2.5 + p) / a * 1.5);
}
bool js( Vec a,Vec b) {
    if (fabs(a.x - b.x) > 0.02 && fabs(a.y - b.x) > 0.02 && fabs(a.z - b.z) > 0.02) {
        return true;
    }
    else {
        return false;
    }
}
void setwindow(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, width, height, 0, -1, 1); // 坐标系设置为像素坐标系
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}
void drawpoint(int x, int y, Vec c) {
    glColor3f(c.x, c.y, c.z);
    glVertex2f(static_cast<float>(x), static_cast<float>(y));
}
//全屏显示
void MaxScreen() {
    HWND Hwnd = GetForegroundWindow();
    ShowWindow(Hwnd, SW_MAXIMIZE);
}
Vec Vecfabs(Vec v) {
    return Vec(fabs(v.x), fabs(v.y), fabs(v.z));
}
#define angle 0.05
#define bc 15
#define DETECT 0
int main(int argc, char* argv[]) {
    int w = 600, h = 600; // 设置图像大小
    //RenderWindow window(VideoMode(800, 600), "绘制一个方形的点");
    if (!glfwInit())
        return -1;
    GLFWwindow* window = glfwCreateWindow(w, h, "Hello World", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }
    // 设置窗口相关回调函数
    glfwSetFramebufferSizeCallback(window, setwindow);
    // 使用OpenGL进行绘图
    glfwMakeContextCurrent(window);
    glMatrixMode(GL_PROJECTION);
    glOrtho(0, w, h, 0, -1, 1); // 坐标系设置为像素坐标系
    glMatrixMode(GL_MODELVIEW);
    MaxScreen();
    float fov = 45; // 视角
    Camera cam(Vec(50, 42, 395.6), (Vec(0, 0, -30) - Vec(50, 42, 395.6)).norm(), fov);
    Vec cx = Vec(w * fov / h, 0, 0); // 计算cx向量
    Vec cy = (cx % cam.d).norm() * fov; // 计算cy向量
    Vec* c = new Vec[w * h]; // 存储每个像素的颜色
    array_view<Vec, 1> cv(w*h,c);
    double seeds[seeds_];
    array_view<double, 1> seedsv(seeds_, seeds);
    int samps = 1;
    int slength1 = sizeof(scenes1) / sizeof(Obj);
    array_view<Obj, 1>scenesv(slength1, scenes1);
    while (!glfwWindowShouldClose(window)) {
        for (int iii = 0; iii < slength1; iii++) {
            cv[iii] = Vec();
        }
        std::random_device rd;  
        std::mt19937 gen(rd()); 
        std::uniform_real_distribution<double> dis(0.0, 2000);  // 创建一个0~1之间的均匀分布
        for (int i = 0; i < seeds_; i++) {
            seeds[i] = dis(gen);
        }
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // 清空颜色缓冲区
        glClear(GL_COLOR_BUFFER_BIT);  // 清空颜色缓冲区
        glBegin(GL_POINTS);//几何图形的起始点和类型
        Vec* c = new Vec[w * h]; // 存储每个像素的颜色
        parallel_for_each(cv.extent, [=](index<1> idx) restrict(amp) {
            int samps1 = samps;
            Vec r,am=Vec(.02,.02,.02); // 存储每个像素的颜色
            for (int sy = 0; sy < 2; sy++) //2x2子像素行
                for (int sx = 0; sx < 2; sx++, r = Vec()) { //2x2子像素列,r用来记录本次获得颜色值
                    for (int s = 0; s < samps1; s++) { //采样
                        //像素发出光线的方向(↓)
                        int x = idx[0] % h;
                        int y = idx[0] / w;
                        Vec d = cam.getRay(x, y, w, h).d;
                        d.y *= -1;
                        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, seedsv, cam,scenesv,slength1,am) * (1. / samps);
                    }
                    cv[idx] = cv[idx] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25; //因为是2*2的子像素,每个结果只占1/4,所以乘0.25
                    int x = idx[0] % h;
                    int y = idx[0] / w;
                    cv[i] = DeNoisy(x, y, i, w, h, samps, cv,2);
                    /*if (cv[idx - 1].x > 0.95 && cv[idx - 1].y > 0.95 && cv[idx - 1].z > 0.95) { //重要性采样
                        samps1 = 0;
                        cv[idx] = Vec(1, 1, 1);
                    }
                    else if (cv[idx - 1].x < 0.05 && cv[idx - 1].y < 0.05 && cv[idx - 1].z < 0.05) {
                        samps1 = 0;
                        cv[idx] = Vec();
                    }
                    else {
                        samps1 = samps;
                    }*/
                    if (cv[idx - 1].x > 0.95 && cv[idx - 1].y > 0.95 && cv[idx - 1].z > 0.95) { //重要性采样
                        samps1 = 1;
                    }
                    else if (cv[idx - 1].x < 0.05 && cv[idx - 1].y < 0.05 && cv[idx - 1].z < 0.05) {
                        samps1 = 1;
                    }
                    else {
                        samps1 = samps;
                    }
                }
            });
        scenesv.synchronize();
        glBegin(GL_POINTS); // 开始绘制点
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                int i = (h - y - 1) * w + x; // 计算像素在数组中的索引
                drawpoint(x, y, cv[i]);
                cv[i] = Vec();
            }
        }
        glEnd(); // 结束绘制点
        glfwSwapBuffers(window);/* 交换前后缓冲区，显示渲染结果 */
        glfwPollEvents();/* 处理窗口的事件 */
        if (GetAsyncKeyState(VK_ESCAPE) & 0x8000) {
            return 0;
        }
        else if (GetAsyncKeyState('A') & 0x8000) {
            samps = 1;
            cam.o = cam.o - cam.right * bc;
        }
        else if (GetAsyncKeyState('D') & 0x8000) {
            samps = 1;
            cam.o = cam.o + cam.right * bc;
        }
        else if (GetAsyncKeyState('W') & 0x8000) {
            samps = 1;
            cam.o = cam.o + cam.d * bc;
        }
        else if (GetAsyncKeyState('S') & 0x8000) {
            samps = 1;
            cam.o = cam.o - cam.d * bc;
        }
        else if (GetAsyncKeyState('E') & 0x8000) {
            samps = 1;
            cam.o = cam.o - cam.up * bc;
        }
        else if (GetAsyncKeyState('Q') & 0x8000) {
            samps = 1;
            cam.o = cam.o + cam.up * bc;
        }
        else if (GetAsyncKeyState(VK_RIGHT) & 0x8000) {
            samps = 1;
            cam.d = VecRotate(cam.d, -angle, cam.up);
            cam.right = VecRotate(cam.right, -angle, cam.up);
            cam.up = VecRotate(cam.up, -angle, cam.up);
        }
        else if (GetAsyncKeyState(VK_LEFT) & 0x8000) {
            samps = 1;
            cam.d = VecRotate(cam.d, angle, cam.up);
            cam.right = VecRotate(cam.right, angle, cam.up);
            cam.up = VecRotate(cam.up, angle, cam.up);
        }
        else if (GetAsyncKeyState(VK_UP) & 0x8000) {
            samps = 1;
            cam.d = VecRotate(cam.d, -angle, cam.right);
            cam.right = VecRotate(cam.right, -angle, cam.right);
            cam.up = VecRotate(cam.up, -angle, cam.right);
        }
        else if (GetAsyncKeyState(VK_DOWN) & 0x8000) {
            samps = 1;
            cam.d = VecRotate(cam.d, angle, cam.right);
            cam.right = VecRotate(cam.right, angle, cam.right);
            cam.up = VecRotate(cam.up, angle, cam.right);
        }
        samps ++;
    }
    debug.print();
    glEnd();//绘制完成
    glfwDestroyWindow(window);
    glfwTerminate();
}