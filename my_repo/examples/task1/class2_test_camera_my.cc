/*
    2020/7/25
    created by: Zhang Yichi
      
    已知：空间三维点坐标、相机的内参矩阵和外参矩阵
    计算：投影后的像素坐标
*/
#include <iostream>
#include <vector>
#include <cmath>
#include "Eigen/Dense"
//#include "math/matrix.h"
//#include "math/vector.h"

//用齐次形式表达二维和三维点坐标，最后一项均为1
typedef Eigen::Matrix<double, 4, 1> point3d;
typedef Eigen::Matrix<double, 3, 1> point2d;

class camera{
private:
    //焦距
    double f;
    //径向畸变系数k1、k2
    double k1, k2;
    //像素坐标系偏移u0、v0
    double u0, v0;
    //内参矩阵
    Eigen::Matrix3d inner_parameter_mat;

    //旋转矩阵
    Eigen::Matrix3d R;
    //平移向量
    Eigen::Matrix<double, 3, 1> t;
    //外参矩阵
    Eigen::Matrix<double, 3, 4> transform_mat;

public:
    camera
    (double f,  double k1, double k2, double u0, double v0, Eigen::Matrix3d R, Eigen::Matrix<double, 3, 1> t){
        inner_parameter_mat.setZero();
        inner_parameter_mat(0,0) = inner_parameter_mat(1,1) = f; 
        inner_parameter_mat(0,2) = u0;
        inner_parameter_mat(1,2) = v0;
        inner_parameter_mat(2,2) = 1;

        transform_mat.block(0,0,3,3).transpose() = R;
        transform_mat.block(0,3,3,1) = t;

        this->k1 = k1;
        this->k2 = k2;
    }

    void print_2dpoint(point2d p){
        std::cout<<"("<<p[0]<<", "<<p[1]<<", "<<p[2]<<")"<< std::endl;
    }

    void print_3dpoint(point3d p){
        std::cout<<"("<<p[0]<<", "<<p[1]<<", "<<p[2]<<", "<<p[3]<<")"<< std::endl;
    }

    void print_MxN_matrix(size_t m, size_t n, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat){
        if(mat.rows() != m || mat.cols() != n){
            std::cout<<"wrong matrix size!"<<std::endl;
            return;
        }
        std::cout<<"[";
        for(size_t i = 0; i < m ; ++i){ 
            for(size_t j = 0; j < n; ++j){   
                std::cout<<mat(i, j)<<"\t";
            }
            std::cout<<"\n";    
        }
        std::cout<<"]";
    }

    point2d projection_3Dto2D(point3d wor_point){

        std::cout<<"wor_point: ";
        print_3dpoint(wor_point);
        std::cout<<std::endl;

        //世界坐标系到相机坐标系
        //3x1 = 3x4 * 4x1
        //cam_point不用齐次表达，是一个三维坐标
        point2d cam_point = transform_mat * wor_point;
        std::cout<<"cam_point: ";
        print_2dpoint(cam_point);
        std::cout<<std::endl;

        //内参与外参矩阵
        std::cout<<"inner_parameter_mat: "<<std::endl;
        print_MxN_matrix(3, 3, inner_parameter_mat);
        std::cout<<"transform_mat: "<<std::endl;
        print_MxN_matrix(3, 4, transform_mat);

        //相机坐标系到像平面坐标系，未纠正径向畸变
        //3x1 = 3x3 * 3x1
        point2d ori_img_point = inner_parameter_mat * cam_point;
        std::cout<<"ori_img_point: ";
        print_2dpoint(ori_img_point);
        std::cout<<std::endl;

        //纠正径向畸变
        int r = pow(ori_img_point[0] , 2) + pow(ori_img_point[1] , 2);
        double distort_parameter = 1 + k1 * r + k2 * pow(r, 2); 
        std::cout<<"distort_parameter: "<<distort_parameter<<std::endl;

        point2d img_point = distort_parameter * ori_img_point;
        return img_point;
    }

};

int main(){
    
    //定义焦距
    double s_f = 1.5;

    //定义径向畸变系数
    double s_k1 = 0.005;
    double s_k2 = 0.005;

    //定义像平面坐标偏移量
    double s_u0 = 1;
    double s_v0 = 1;

    //定义旋转矩阵（正交矩阵）
    Eigen::Matrix3d s_R;
    s_R(0,0) = 2; s_R(0,1) = 0;   s_R(0,2) = 0;
    s_R(1,0) = 0; s_R(1,1) = 1.5; s_R(1,2) = 0; 
    s_R(2,0) = 0; s_R(2,1) = 0;   s_R(2,2) = -1; 

    //平移向量
    Eigen::Matrix<double, 3, 1> s_t;
    s_t[0] = 1; s_t[1] = -2; s_t[2] = 3;

    //初始化赋值
    camera my_cam(s_f, s_k1, s_k2, s_u0, s_v0, s_R, s_t);

    //输入三维点
    point3d input_3dpoint;
    input_3dpoint[3] = 1.0;

    input_3dpoint[0] = 0.5;
    input_3dpoint[1] = -1;
    input_3dpoint[2] = 1;

    point2d output_img_point = my_cam.projection_3Dto2D(input_3dpoint);
    
    my_cam.print_2dpoint(output_img_point);
    return 0;
}

