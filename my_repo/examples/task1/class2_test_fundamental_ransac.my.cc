/*
    2020/7/29
    created by: Zhang Yichi

    已知两对极约束的相机像素的n对特征点配对
    采用ransac计算基础矩阵

    NOTE:
    txt文件中添加了离群点 

*/
#include <iostream>
#include <cmath>     //for pow(); ceil()
#include <algorithm> //for sort();
#include <assert.h>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include "Eigen/Dense"
#include "class2_test_fundamental_8_point_my.h"

/* 
    struct point_2d{
        double x;
        double y;
        double z = 1.0;
    };

    typedef std::vector <point_2d> point_sets_vec;
*/
struct point_sets_pair_vec
{
    point_sets_vec set_a;
    point_sets_vec set_b;
};

/*
    存放一个基础矩阵及其对应的内点个数
*/
struct inlinerNum_FundaMatrix{
    int inlinerNum;
    Eigen::Matrix3d FundaMatrix;
};

/*
    vector存放所有基础矩阵及其相应的内点个数
    便于得到内点个数最多的基础矩阵
*/
struct inlinerNum_FundaMatrix_vec{
    std::vector<int> inlinerNum_vec;
    std::vector<Eigen::Matrix3d> FundaMatrix_vec;
};


/*加载归一化后的匹配点对*/
/*
    file_dir 文件路径
*/
point_sets_pair_vec read_pointsets_file(std::string file_dir)
{
    //点对容器
    point_sets_pair_vec corr_pairs;

    std::ifstream ifs;
    ifs.open(file_dir, std::ios::in);
    assert(ifs.is_open());
    std::cout<<"successfully read file"<<std::endl;
    
    std::string data;    //存放读取的数据
    int line = 0;        //文件行数
    int point_pairs_num = 0; //点对个数

    while(getline(ifs, data))
    {
        std::stringstream ss_data(data);
        if(line == 0)
        {
            //txt文件中第一行是特征点对的个数,输出个数确定点对容器大小
            ss_data >> point_pairs_num;
            corr_pairs.set_a.resize(point_pairs_num);
            corr_pairs.set_b.resize(point_pairs_num);
            line++;
            continue;
        }
        if(line > 0)
        {
            //输出每行的4个数据赋给点对容器
            ss_data >>  corr_pairs.set_a[line -1].x
                    >>  corr_pairs.set_a[line -1].y
                    >>  corr_pairs.set_b[line -1].x
                    >>  corr_pairs.set_b[line -1].y;
        }
        line++;   
    }
    std::cout<<"文件中共有： "<< corr_pairs.set_a.size() <<"组点对"<<std::endl;
    return corr_pairs;
}

/*用于RANSAC采样成功所需要的采样次数*/
/*
    p -- 估计的内点的概率
    K --拟合模型需要的样本个数，对应基础矩阵num_samples=8
    z  -- 预期的采样成功的概率
                            log(1-z)
        需要的采样次数 M = -----------
                            log(1-p^K)
    Example: For p = 50%, z = 99%, n = 8: M = log(0.001) / log(0.99609) = 1176.
    需要采样1176次从而保证RANSAC的成功率不低于0.99.
*/
int cal_ransac_iterations_num(double p, int k, double z)
{
    double m = std::log(1 -z) / std::log(1 - (pow(p, k)));
    int num = static_cast<int> (std::ceil(m));
    std::cout<<"ransac需要的采样次数："<<num<<std::endl;
    return num;    
}


/*sampson距离计算当前基础矩阵模型的内点个数*/
/*
    sampson_threshold            判断阈值
    input_corr_pairs             读取的特征点对集
    input_rand_8_pairs_funda_mat 输入基础矩阵

    Sampson distance = (x2'Fx1)^2 / ( (Fx1)_1^2 + (Fx1)_2^2 + (x2'F)_1^2 + (x2'F)_2^2 )
*/
int cal_inliner_num_with_sampson(double sampson_threshold, 
                                point_sets_pair_vec const& input_corr_pairs,
                                Eigen::Matrix3d const& input_rand_8_pairs_funda_mat)
{
    //将input_corr_pairs转为矩阵类型便于计算
    Eigen::Matrix<double, 3, 1> x1;
    Eigen::Matrix<double, 1, 3> x2;

    //计数
    int count = 0;

    for(int i = 0; i < input_corr_pairs.set_a.size(); ++i){
        x1[0] = input_corr_pairs.set_a[i].x;
        x1[1] = input_corr_pairs.set_a[i].y;
        x1[2] = input_corr_pairs.set_a[i].z;  // = 1.0

        x2[0] = input_corr_pairs.set_b[i].x;
        x2[1] = input_corr_pairs.set_b[i].y;
        x2[2] = input_corr_pairs.set_b[i].z;  // = 1.0
        
        //分子
        double numerator = pow((x2 * input_rand_8_pairs_funda_mat * x1),2);
        //分母
        double denominator = (
            pow((input_rand_8_pairs_funda_mat * x1)[0], 2) +
            pow((input_rand_8_pairs_funda_mat * x1)[1], 2) +
            pow((x2 * input_rand_8_pairs_funda_mat)[0], 2) +
            pow((x2 * input_rand_8_pairs_funda_mat)[1], 2) );  

        //阈值判断
        if((numerator / denominator) < sampson_threshold)
        {
            count++;
        }
    }
    return count;
}


/*ransac迭代*/
/*
    sampson_threshold 判断阈值
    iter_num          ransac最大采样次数
    input_corr_pairs  读取的特征点对集 

    return：（具有最大内点数的基础矩阵 ，最大内点数）  
*/
inlinerNum_FundaMatrix ransac_iteration_process(double sampson_threshold, 
                                        int iter_num, 
                                        point_sets_pair_vec const& input_corr_pairs)
{
    //读取的特征点对集个数
    int input_corr_pairs_num = input_corr_pairs.set_a.size();

    //8组点对
    point_sets_pair_vec random_8_pairs;
    random_8_pairs.set_a.resize(8);
    random_8_pairs.set_b.resize(8);
    std::set<int> random_8_idx;

    //基础矩阵
    Eigen::Matrix3d rand_8_pairs_funda_mat;

    //内点个数_基础矩阵容器
    inlinerNum_FundaMatrix_vec num_mat_vec;

    for(int i = 0; i < iter_num; ++i)
    {
        //利用set容器 随机生成8个不重复的、位于（0～input_corr_pairs_num-1）的整数
        while(random_8_idx.size() < 8)
        {  
            int rand_num = rand() % input_corr_pairs_num;
            random_8_idx.insert(rand_num); 
        }
        
        //相应的8组点对
        int k = 0; 
        for(std::set<int>::iterator it = random_8_idx.begin(); it != random_8_idx.end(); ++it)
        {
            random_8_pairs.set_a[k].x = input_corr_pairs.set_a[*it].x;
            random_8_pairs.set_a[k].y = input_corr_pairs.set_a[*it].y;
            random_8_pairs.set_b[k].x = input_corr_pairs.set_b[*it].x;
            random_8_pairs.set_b[k].y = input_corr_pairs.set_b[*it].y; 
            ++k;  
        } 
        random_8_idx.clear();

        //8点法估计基础矩阵
        rand_8_pairs_funda_mat = cal_funda_mat(random_8_pairs.set_a, random_8_pairs.set_b);

        //计算sampson距离，计算内点个数
        int inliner_num = cal_inliner_num_with_sampson(sampson_threshold, input_corr_pairs, rand_8_pairs_funda_mat);
        std::cout<<"inliner: "<<inliner_num<<std::endl;

        //存放
        num_mat_vec.inlinerNum_vec.push_back(inliner_num);
        num_mat_vec.FundaMatrix_vec.push_back(rand_8_pairs_funda_mat);
    }

    //inlinerNum_vec最大值索引
    std::vector<int>::iterator max_iter = std::max_element(std::begin(num_mat_vec.inlinerNum_vec), std::end(num_mat_vec.inlinerNum_vec));
    int max_idx = std::distance(std::begin(num_mat_vec.inlinerNum_vec), max_iter);
    
    std::cout<<"max_idx: "<< max_idx <<std::endl;
    std::cout<<"*max_iter: "<< *max_iter <<std::endl;

    Eigen::Matrix3d max_inliner_funda_mat = num_mat_vec.FundaMatrix_vec[max_idx];
    inlinerNum_FundaMatrix max_inlinerNum_FundaMatrix;
    max_inlinerNum_FundaMatrix.inlinerNum = *max_iter;
    max_inlinerNum_FundaMatrix.FundaMatrix = max_inliner_funda_mat;

    return max_inlinerNum_FundaMatrix;
}

/*对读取的特征点集采样，提取最大内点时所有的inliner*/
/*
   double sampson_threshold        判断阈值
   input_largest_inliner_funda_mat 具有最大内点数的基础矩阵 
   int max_inlinerNum              最大内点个数
   input_corr_pairs                读取的特征点对集 
   
   返回：内点数最大时所有的内点点对
*/
point_sets_pair_vec sample_corr_pairs_to_MaxInliner(double sampson_threshold, 
                                                Eigen::Matrix3d input_largest_inliner_funda_mat,
                                                int max_inlinerNum,
                                                point_sets_pair_vec const& input_corr_pairs)
{
     //将input_corr_pairs转为矩阵类型便于计算
    Eigen::Matrix<double, 3, 1> x1;
    Eigen::Matrix<double, 1, 3> x2;

    //读取的特征点对集个数
    size_t input_corr_pairs_num = input_corr_pairs.set_a.size();

    //存放inliner point_pair容器
    point_sets_pair_vec inliner_point_pair;
    std::vector<double> vec_set_a_x, vec_set_a_y, vec_set_b_x, vec_set_b_y;

    for(int i = 0; i < input_corr_pairs.set_a.size(); ++i){
        x1[0] = input_corr_pairs.set_a[i].x;
        x1[1] = input_corr_pairs.set_a[i].y;
        x1[2] = input_corr_pairs.set_a[i].z;  // = 1.0

        x2[0] = input_corr_pairs.set_b[i].x;
        x2[1] = input_corr_pairs.set_b[i].y;
        x2[2] = input_corr_pairs.set_b[i].z;  // = 1.0
        
        //分子
        double numerator = x2 * input_largest_inliner_funda_mat * x1;
        //分母
        double denominator = (
            pow((input_largest_inliner_funda_mat * x1)[0], 2) +
            pow((input_largest_inliner_funda_mat * x1)[1], 2) +
            pow((x2 * input_largest_inliner_funda_mat)[0], 2) +
            pow((x2 * input_largest_inliner_funda_mat)[1], 2) );  

        //阈值判断
        if((numerator / denominator) < sampson_threshold)
        {
           vec_set_a_x.push_back(input_corr_pairs.set_a[i].x);
           vec_set_a_y.push_back(input_corr_pairs.set_a[i].y);
           vec_set_b_x.push_back(input_corr_pairs.set_b[i].x);
           vec_set_b_y.push_back(input_corr_pairs.set_b[i].y);   
        }    
    }
    for(int j = 0; j < max_inlinerNum; ++j)
    {
        inliner_point_pair.set_a[j].x = vec_set_a_x[j];
        inliner_point_pair.set_a[j].y = vec_set_a_y[j];
        inliner_point_pair.set_b[j].x = vec_set_b_x[j];
        inliner_point_pair.set_b[j].y = vec_set_b_y[j];
    }
    vec_set_a_x.clear(), vec_set_a_y.clear(), vec_set_b_x.clear(), vec_set_b_y.clear();

    return inliner_point_pair;
}
    
int main(int argc, char* argv[])
{   
    //读取文件
    std::string file_dir = 
    "/home/zhangyichi/VSCode/ImageBasedModelling/ImageBasedModelling/examples/task2/correspondences.txt";
    
    point_sets_pair_vec  corr_pairs = read_pointsets_file(file_dir);
    std::cout<<"文件读取完成"<<std::endl;

    //计算ransac最大采样次数
    /*
        p -- 估计的内点的概率
        K --拟合模型需要的样本个数，对应基础矩阵num_samples=8
        z  -- 预期的采样成功的概率
                                log(1-z)
            需要的采样次数 M = -----------
                                log(1-p^K)
        Example: For p = 50%, z = 99%, n = 8: M = log(0.001) / log(0.99609) = 1176.
        需要采样1176次从而保证RANSAC的成功率不低于0.99.
    */
    const double p = 0.5;
    const int k = 8;
    double z = 0.99;

    int m = cal_ransac_iterations_num(p, k, z);
    std::cout<<"最大采样次数计算完成"<<std::endl;

    //更新随机数
    time_t t;
    srand(time(&t));

    //ransac迭代
    const double sampson_threshold = 0.0015;
    inlinerNum_FundaMatrix num_mat = ransac_iteration_process(sampson_threshold, m, corr_pairs);
    std::cout<<"ransac迭代完成"<<std::endl;

    //具有最大内点数的基础矩阵
    Eigen::Matrix3d max_inliner_fundMat = num_mat.FundaMatrix;
    std::cout<<"基础矩阵："<<std::endl;
    print_MxN_matrix(3, 3, max_inliner_fundMat);
    std::cout<<"具有最大内点个数： "<<num_mat.inlinerNum<<std::endl;
}
