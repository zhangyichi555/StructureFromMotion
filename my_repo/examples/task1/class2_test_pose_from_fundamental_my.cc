/*
    2020/7/31
    created by: Zhang Yichi

    已知：已知两对极约束的相机像素的n对特征点配对
    通过基础矩阵计算本征矩阵
    通过本征矩阵计算合适的两相机的姿态（R，t）
    相机姿态是将相机O1定为世界坐标系，计算O2相对于O1的外参矩阵

    计算相机姿态时，需要使用4个备选姿态（R，t）和一对匹配点进行三角测量得到空间三维点
    这里三角量测计算三维点采用直接线性变换法

*/

#include <iostream>
#include <vector>
#include <assert.h>
#include "Eigen/Dense"
#include "Eigen/SVD"

 
struct point_2d{
    double x;
    double y;
    double z = 1.0;
};

typedef std::vector <point_2d> point_sets_vec;

struct point_sets_pair_vec
{
    point_sets_vec set_a;
    point_sets_vec set_b;
};

//旋转矩阵和平移向量
typedef std::pair<Eigen::Matrix3d, Eigen::Matrix<double,3,1> > R_t_pair;
typedef std::vector<R_t_pair> vec_R_t_pair;



/*通过对极约束的基础矩阵，计算本征矩阵E*/
/*
    funda_mat 基础矩阵
    k1 相机1的内参矩阵
    k2 相机2的内参矩阵
*/
Eigen::Matrix3d cal_essential_mat  (Eigen::Matrix3d funda_mat, 
                                    Eigen::Matrix3d k1,
                                    Eigen::Matrix3d k2
                                    )
{
    Eigen::Matrix3d essential_mat = k2.transpose() * funda_mat * k1; 

    //奇异值约束：（（a1+a2）/2 , （a1+a2）/2, 0）
    Eigen::JacobiSVD < Eigen::MatrixXd > svd(essential_mat, Eigen::ComputeThinU | Eigen::ComputeThinV ); 
    Eigen::Vector3d value_vector = svd.singularValues();

    //重构奇异值矩阵
    Eigen::Matrix3d constrained_value_mat;
    constrained_value_mat.setZero();
    constrained_value_mat(0,0) = constrained_value_mat(1,1) = (value_vector[0]+value_vector[1]) / 2;

    //重构本征矩阵
    Eigen::Matrix3d constrained_essential_mat = svd.matrixU() * constrained_value_mat * svd.matrixV().transpose();

    return constrained_essential_mat;
}


/*通过三角量测的方法，使用直接线性变换法计算两相机像素特征点对的空间三维点*/
/*
    point_pair_vec  相机像素特征点对，这里只需要一对像素点
    k1              相机1的内参矩阵
    k2              相机2的内参矩阵
    R_t_1           相机1的外参矩阵
    R_t_2           相机2的外参矩阵
*/
Eigen::Matrix<double, 3, 1> triangulation  (Eigen::Matrix3d k1,
                                            Eigen::Matrix3d k2,
                                            point_sets_pair_vec point_pair_vec,
                                            Eigen::Matrix<double,3,4> R_t_1,
                                            Eigen::Matrix<double,3,4> R_t_2
                                            )
{
    //这里一对相机只有一对特征匹配点，参数矩阵A为4x4
    Eigen::Matrix4d A;

    //参考task3文档中的推导过程,填充参数矩阵A
    A.row(0) = point_pair_vec.set_a[0].x * R_t_1.row(2) - R_t_1.row(0);
    A.row(1) = point_pair_vec.set_a[0].y * R_t_1.row(2) - R_t_1.row(1);
    A.row(2) = point_pair_vec.set_b[0].x * R_t_2.row(2) - R_t_2.row(0);
    A.row(3) = point_pair_vec.set_b[0].y * R_t_2.row(2) - R_t_2.row(1);

    //最小二乘法求三维点坐标
    Eigen::JacobiSVD <Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV); 

    return svd.matrixV().col(3);      
}


/*通过本征矩阵得到的4个相机外参矩阵，分别计算其空间三维点，选择三维点在两相机中深度同时为正的外参矩阵*/
/*
    constrained_essential_mat 本征矩阵
    point_pair      相机像素特征点对
    k1              相机1的内参矩阵
    k2              相机2的内参矩阵
*/
R_t_pair select_proper_transfrom_mat(Eigen::Matrix3d constrained_essential_mat,
                                    point_sets_pair_vec point_pair_vec,
                                    Eigen::Matrix3d k1,
                                    Eigen::Matrix3d k2
                                    )
{
    //存放4个外参矩阵
    vec_R_t_pair transform_mat_container(4);

    //参数矩阵
    Eigen::Matrix3d s1;
    s1.setZero();
    s1(0,1) = -1; s1(1,0) = 1; s1(2,2) = 1;

    //三维点
    Eigen::Matrix<double,3,1> point_3d;

    //计算三维点时，将相机O1定为世界坐标系
    Eigen::Matrix<double,3,4> R_t_1;
    R_t_1.setZero();
    R_t_1(0,0) = R_t_1(1,1) = R_t_1(2,2) = 1;
    Eigen::Matrix<double,3,4> R_t_2;

    //两相机在世界坐标系中的朝向
    Eigen::Matrix<double,3,1> d1;
    d1[0] = d1[1] = 0;  d1[2] = 1;
    Eigen::Matrix<double,3,1> d2;   //旋转矩阵的第三行

    //正确的相机位姿
    R_t_pair correct_R_t_pair;
    
    //按照task2文档，计算四种外参矩阵,为相机O2的外参矩阵
    Eigen::JacobiSVD < Eigen::MatrixXd > svd(constrained_essential_mat, Eigen::ComputeThinU | Eigen::ComputeThinV ); 

    transform_mat_container[0].first = svd.matrixU() * s1 * svd.matrixV().transpose();
    transform_mat_container[0].second = svd.matrixU().col(2);

    transform_mat_container[1].first = svd.matrixU() * s1 * svd.matrixV().transpose();
    transform_mat_container[1].second = -svd.matrixU().col(2);

    transform_mat_container[2].first = svd.matrixU() * s1.transpose() * svd.matrixV().transpose();
    transform_mat_container[2].second = svd.matrixU().col(2);

    transform_mat_container[3].first = svd.matrixU() * s1.transpose() * svd.matrixV().transpose();
    transform_mat_container[3].second = -svd.matrixU().col(2);

    int count = 0;
    for(size_t i = 0; i < 4; ++i)
    {
        //pair类型转为matrix
        R_t_2.block(0,0,3,3) = transform_mat_container[i].first;
        R_t_2.block(0,3,3,1) = transform_mat_container[i].second;

        //计算三维点
        point_3d = triangulation(k1, k2, point_pair_vec, R_t_1, R_t_2);

        //旋转矩阵的第三行
        d2 = transform_mat_container[i].first.row(2).transpose();

       //按照task2文档，判断两向量点集是否同时为正
       if(
            point_3d.transpose() * d1 > 0 &&
            ( point_3d + 
              transform_mat_container[i].first.transpose() * 
              transform_mat_container[i].second
            ).transpose() * d2 > 0
        )
       {
           correct_R_t_pair.first = transform_mat_container[i].first;
           correct_R_t_pair.second = transform_mat_container[i].second;
           count++;
       }  
    }
    assert(count == 1);
    return correct_R_t_pair;
}
