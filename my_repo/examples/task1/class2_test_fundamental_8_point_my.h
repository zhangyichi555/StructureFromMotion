/*
    2020/7/27
    created by: Zhang Yichi

    已知两对极约束的相机像素的n对特征点配对
    采用直接线性变换法计算基础矩阵

*/
#pragma once
#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/SVD"

struct point_2d{
    double x;
    double y;
    double z = 1.0;
};

typedef std::vector <point_2d> point_sets_vec;


Eigen::Matrix3d cal_funda_mat(point_sets_vec &point0, point_sets_vec &point1)
{

    int point_size = point0.size();

    //按照对极约束方程变换：Af = 0, f是基础矩阵9个元素组成的向量
    std::vector<double> f(9);
    Eigen::Matrix<double, Eigen::Dynamic, 9> A (point_size, 9);
    Eigen::Matrix3d fundamental_mat;

    //矩阵A的构造:
    //point0[0] (u0, v0, 1) ------match----- point1[0] (u1, v1, 1)
    //point0[1]             ------match----- point1[1]
    //......

    //A = 
    //rows(0)[u0*u0, u0*v1, u0, v0*u1, v0*v1, v0, u1, v1, 1]
    //rows(1)
    //......
    for(size_t i = 0; i < point_size; ++i){
        A(i, 0) = point0[i].x * point0[i].x;
        A(i, 1) = point0[i].x * point1[i].y;
        A(i, 2) = point0[i].x;
        A(i, 3) = point0[i].y * point1[i].x;
        A(i, 4) = point0[i].y * point1[i].y;
        A(i, 5) = point0[i].y;
        A(i, 6) = point1[i].x;
        A(i, 7) = point1[i].y;
        A(i, 8) = 1;
    }
     
    //JacobiSVD 分解中，奇异值是按照降序进行排序
    //V的列向量是A的右奇异向量，同时也是A^T * A的特征向量
    //当有8个特征点对，唯一解，A的最小奇异值S(8, 8)对应的奇异向量为最优解(V的第8列向量)
    //当特征点对大于8个，最小二乘法，A^T * A的最小奇异值对应的奇异向量为最优解
    Eigen::JacobiSVD < Eigen::MatrixXd > svd1(A, Eigen::ComputeThinU | Eigen::ComputeThinV ); 
    for(int i = 0; i <= 8; ++i){
        f[i] = svd1.matrixV()(i, 7);
    }
   
    //还原基础矩阵
    fundamental_mat(0,0) = f[0]; fundamental_mat(0,1) = f[1]; fundamental_mat(0,2) = f[2]; 
    fundamental_mat(1,0) = f[3]; fundamental_mat(1,1) = f[4]; fundamental_mat(1,2) = f[5]; 
    fundamental_mat(2,0) = f[6]; fundamental_mat(2,1) = f[7]; fundamental_mat(2,2) = f[8]; 
    
    //奇异值约束
    //JacobiSVD 分解中，奇异值默认以Eigen::Vector3d列向量输出
    Eigen::JacobiSVD<Eigen::MatrixXd> svd2(fundamental_mat, Eigen::ComputeThinU | Eigen::ComputeThinV ); 
    Eigen::Vector3d value_vector = svd2.singularValues();
    
    //重构奇异值矩阵，令第三个奇异值为0
    Eigen::Matrix3d constrainted_value_mat;
    constrainted_value_mat.setZero();
    constrainted_value_mat(0,0) = value_vector[0];
    constrainted_value_mat(1,1) = value_vector[1];

    //重构基础矩阵
    fundamental_mat = svd2.matrixU() * constrainted_value_mat * svd2.matrixV().transpose();
    
    return fundamental_mat;
}

 void print_MxN_matrix(size_t m, size_t n, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat){
        if(mat.rows() != m || mat.cols() != n){
            std::cout<<"input wrong matrix size when printing!"<<std::endl;
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


// int main(int argc, char* argv[]){

//     //输入的特征点对个数
//     size_t point_pairs_num = 8;

//     point_sets_vec pset1(point_pairs_num);
//     pset1[0].x = 0.180123 ; pset1[0].y= -0.156584; pset1[0].z=1.0;
//     pset1[1].x = 0.291429 ; pset1[1].y= 0.137662 ; pset1[1].z=1.0;
//     pset1[2].x = -0.170373; pset1[2].y= 0.0779329; pset1[2].z=1.0;
//     pset1[3].x = 0.235952 ; pset1[3].y= -0.164956; pset1[3].z=1.0;
//     pset1[4].x = 0.142122 ; pset1[4].y= -0.216048; pset1[4].z=1.0;
//     pset1[5].x = -0.463158; pset1[5].y= -0.132632; pset1[5].z=1.0;
//     pset1[6].x = 0.0801864; pset1[6].y= 0.0236417; pset1[6].z=1.0;
//     pset1[7].x = -0.179068; pset1[7].y= 0.0837119; pset1[7].z=1.0;

//     point_sets_vec pset2(point_pairs_num);
//     pset2[0].x = 0.208264 ; pset2[0].y= -0.035405 ; pset2[0].z = 1.0;
//     pset2[1].x = 0.314848 ; pset2[1].y=  0.267849 ; pset2[1].z = 1.0;
//     pset2[2].x = -0.144499; pset2[2].y= 0.190208  ; pset2[2].z = 1.0;
//     pset2[3].x = 0.264461 ; pset2[3].y= -0.0404422; pset2[3].z = 1.0;
//     pset2[4].x = 0.171033 ; pset2[4].y= -0.0961747; pset2[4].z = 1.0;
//     pset2[5].x = -0.427861; pset2[5].y= 0.00896567; pset2[5].z = 1.0;
//     pset2[6].x = 0.105406 ; pset2[6].y= 0.140966  ; pset2[6].z = 1.0;
//     pset2[7].x = -0.15257 ; pset2[7].y= 0.19645   ; pset2[7].z = 1.0;
 
//     if(point_pairs_num < 8){
//         std::cout<<"at least 8 point pairs to calculate fundamental matrix!"<<std::endl;
//         return -1;
//     }
    
//     Eigen::Matrix3d funda_mat = cal_funda_mat(pset1, pset2);
//     std::cout<<"svd success"<<std::endl;
//     //print_MxN_matrix(3, 3, funda_mat);

//     return 0;
// }
