#include <vector>
#include <algorithm>
#include <math.h>
#include "opencv2/opencv.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <fstream>
#include <string>
#include <iostream>

#include <math.h>
#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>

#include "norm_4pt.h"
#include "norm_5pt.h"
#include "common.h"


int main()
{
    //读取车牌图像
    std::string img_path = "../img/tmp.jpg";
    cv::Mat img = cv::imread(img_path,1);
    std::cout<<"读取成功\n";
    //车牌图像的四个角,标的可能不太准
    std::vector<cv::Point2f> src = 
    {
        {77.0,113},
        {258,33},
        {267,102},
        {86,186}
    };
    //要归一化的4个目标点
    const float std_w = 440.0f, std_h = 140.0f;
    std::vector<cv::Point2f> points_dst = 
    {
        {.0f, .0f},
        {std_w, .0f},
        {std_w, std_h},
        {.0f, std_h}};

    //opencncv求取的透视变换矩阵
    cv::Mat M = cv::getPerspectiveTransform(src,points_dst);
 
    cv::Mat rectified;
    cv::warpPerspective(img, rectified, M, cv::Size(int(std_w), int(std_h)));
    std::string save_path = "../img/result_cv.jpg";
    cv::imwrite(save_path, rectified);
    std::cout<<M<<"\n";
    std::vector<std::vector<float>> _src = 
    {
        {77.0,113},
        {258,33},
        {267,102},
        {86,186}
    };
    std::vector<std::vector<float>> _dst = 
    {
        {.0f, .0f},
        {std_w, .0f},
        {std_w, std_h},
        {.0f, std_h}
    }; 
    std::vector<std::vector<float> > AB_M(3,std::vector<float>(3));
    // 本人实现的求取的 透视变换矩阵
    AB_getPerspectiveTransform(_src,_dst,AB_M);
    // 验证透视变换矩阵求取的是否正确
    for (int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
            std::cout<<AB_M[i][j]<<"\t";
        std::cout<<"\n";
    }
    int img_w = img.cols;
    int img_h = img.rows;
    //根据透视变换矩阵计算归一化后的图像
    float *img_f = new float[img_h*img_w*3];
    float *img_dst = new float[440*140*3];
    cv::Mat rectified_wpy(140,440,CV_8UC3);
    RGB2float(img.data,img_w,img_h,img_f);
    warpAffineLinear(img_f,img_dst,img_w,img_h,440,140,3,AB_M);
    Float2RGB(img_dst,440,140,rectified_wpy.data);
    delete []img_f;
    delete []img_dst;
    // 验证车牌是否由一般四边形归一化为矩形了
    std::string save_path_wpy = "../img/result_wpy.jpg";
    cv::imwrite(save_path_wpy, rectified_wpy);


    //验证五点归一化

    std::string jinsha_img_path = "../img/jinsha.jpg";
    cv::Mat jinsha_img = cv::imread(jinsha_img_path,1);
    if(jinsha_img.data!=NULL)
    {
        std::cout<<"load jinsha success!\n";
    }
    //原始人脸五点
    std::vector<std::vector<float>> face_src_points =
    {
        {181.f,123.f},
        {246.f,124.f},
        {211.f,163.f},
        {189.f,195.f},
        {238.f,193.f}
    } ;
    //目标人脸5点
    std::vector<std::vector<float>> face_dst_points =
    {
        {38.2946,51.6963},
        {73.5318,51.5014},
        {56.0252,71.7366},
        {41.5493,92.3655},
        {70.7299,92.2041}
    } ;
    float *jinsha_img_f = new float[jinsha_img.cols*jinsha_img.rows*3];
    RGB2float(jinsha_img.data,jinsha_img.cols,jinsha_img.rows,jinsha_img_f);

    float *face_img_f = new float[112*112*3];
    NormImg_5_points(jinsha_img_f,jinsha_img.cols,jinsha_img.rows,3,
                face_src_points,face_dst_points,face_img_f);
    cv::Mat face_img(112,112,CV_8UC3);
    Float2RGB(face_img_f,112,112,face_img.data);
    std::string face_save_path = "../img/jinsha_face.jpg";
    cv::imwrite(face_save_path,face_img);

    delete []jinsha_img_f;
    delete []face_img_f;

    return 0;
}