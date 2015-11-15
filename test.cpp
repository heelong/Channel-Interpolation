#include "kdtree.h"
#include <iostream>
#include <string>
#include <fstream>
#include <strstream>
#include <sstream>
#include <iomanip>
#include <time.h>

/*反距离加权函数,参数为指向邻域内点的容器的迭代器，和各点到查询点之间距离的容器的迭代器*/
double IIDV(std::vector<KDTree::point>::iterator XYZ_begin,std::vector<KDTree::point>::iterator XYZ_end,
	std::vector<double>::iterator d_begin,std::vector<double>::iterator d_end)
{
	double _sumx=0,_sumy=0;
	std::vector<double>::iterator id=d_begin;
	for(std::vector<KDTree::point>::iterator ix=XYZ_begin;ix!=XYZ_end&&id!=d_end;++ix)
	{
		double a=*id++;
		_sumx+=1 /(a*a);//计算权
		_sumy+=ix->z/(a*a);
	}
	return _sumy/_sumx;//返回加权距离
}
// 矩阵求逆
//矩阵a，n行n列，求逆后将逆矩阵存入a中
void inv(double *a,int n)
{ 
	int i,j,k;
	for(k=0;k<n;k++)
	{
		for(i=0;i<n;i++)
		{
			if(i!=k)
				*(a+i*n+k)=-*(a+i*n+k)/(*(a+k*n+k));
		}
		*(a+k*n+k)=1/(*(a+k*n+k));
		for(i=0;i<n;i++)
		{
			if(i!=k)
			{
				for(j=0;j<n;j++)
				{
					if(j!=k)
						*(a+i*n+j)+=*(a+k*n+j)* *(a+i*n+k);
				}
			}
		}
		for(j=0;j<n;j++)
		{
			if(j!=k)
				*(a+k*n+j)*=*(a+k*n+k);
		}
	}
}
// 矩阵相乘
//矩阵m1，i_1行j_12列，
//矩阵m2,j_12行j_2列
//结果矩阵result,i_1行j_2列
void mult(double *m1,double *m2,double *result,int i_1,int j_12,int j_2)
{ 
	int i,j,k;
	for(i=0;i<i_1;i++)
		for(j=0;j<j_2;j++)
		{
			result[i*j_2+j]=0.0;
			for(k=0;k<j_12;k++)
				result[i*j_2+j]+=m1[i*j_12+k]*m2[j+k*j_2];
		}
		return;
}

/*曲面拟合内插,参数分别为求取二次曲面的六个点的XYZ坐标值得迭代器，和待求高程点的XY坐标*/
double SFI(std::vector<KDTree::point>::iterator XYZ_begin,std::vector<KDTree::point>::iterator XYZ_end,std::pair<double,double> coordinate)
{
	double *_eigenX=new double[6*6];//六行六列的矩阵
	double *_eigenY=new double[6];					//六行一列的矩阵
	double *_eigenFac=new double[6];
	int i=0;
	for(std::vector<KDTree::point>::iterator ix=XYZ_begin;ix!=XYZ_end;++ix)
	{
		_eigenX[i*6+0]=ix->x*ix->x;
		_eigenX[i*6+1]=ix->x*ix->y;
		_eigenX[i*6+2]=ix->y*ix->y;
		_eigenX[i*6+3]=ix->x;
		_eigenX[i*6+4]=ix->y;
		_eigenX[i*6+5]=1;
		_eigenY[i]=ix->z;
		++i;
	}
	inv(_eigenX,6);
	mult(_eigenX,_eigenY,_eigenFac,6,6,1);

	double x=coordinate.first;
	double y=coordinate.second;

	return _eigenFac[0]*x*x+_eigenFac[1]*x*y+_eigenFac[2]*y*y+_eigenFac[3]*x+_eigenFac[4]*y+_eigenFac[5];
	delete[]_eigenX;delete[]_eigenY;delete[]_eigenFac;
}
int main(int argc, char *argv[])
{
	double dur;
    clock_t start,end;
    start = clock();

	KDTree::ExamplarSet exm_set;
	std::vector<KDTree::point> _cloudOriginal;
	std::vector<KDTree::point> _cloudInterpolation;
	std::vector<KDTree::point> _cloud;
	std::vector<KDTree::point> _cloud2;
	std::ifstream infile1("E:\\PatternRecognition\\data\\ChannelData\\OriginalSegmentation2.txt");
	std::ifstream infile2("E:\\PatternRecognition\\data\\ChannelData\\InterpolationSegmentation2.txt");
	std::ofstream outfile1("E:\\PatternRecognition\\data\\ChannelData\\1.txt");
	std::ofstream outfile2("E:\\PatternRecognition\\data\\ChannelData\\2.txt");
	if(!infile1)
	{
		std::cout<<"数据文件打开失败！"<<std::endl;
		return 0;
	}
	if(!infile2)
	{
		std::cout<<"数据文件打开失败！"<<std::endl;
		return 0;
	}
	std::string line;
	float a;
	while(getline(infile1,line))
	{	
		double a1,a2,a3;
		std::istringstream l_stream1(line);
		KDTree::point _xyz;
		l_stream1>>a>>a>>a1>>a2>>a3;
		_xyz.x=a1-603000;
		_xyz.y=a2-3350000;
		_xyz.z=a3;
//		outfile1<<_xyz.x<<" "<<_xyz.y<<" "<<_xyz.z<<std::endl;
		_cloudOriginal.push_back(_xyz);
	}
	while(getline(infile2,line))
	{	
		double a1,a2,a3;
		std::istringstream l_stream1(line);
		KDTree::point _xyz;
		l_stream1>>a>>a>>a1>>a2>>a3;
		_xyz.x=a1-603000;
		_xyz.y=a2-3350000;
		_xyz.z=a3;
//		outfile2<<_xyz.x<<" "<<_xyz.y<<" "<<_xyz.z<<std::endl;
		_cloudInterpolation.push_back(_xyz);
	}
	exm_set.createfromcloud(_cloudOriginal.begin(),_cloudOriginal.end());
	KDTree::KDTree kdtree;
	kdtree.create(exm_set);

	/*反距离加权*/
	for(std::vector<KDTree::point>::iterator ix=_cloudInterpolation.begin();ix!=_cloudInterpolation.end();++ix)
	{
		KDTree::_Examplar exm(3);//待插值的点
		KDTree::point _axyz;//用来存储插值后的点云
		_axyz.x=ix->x;_axyz.y=ix->y;
		exm[0]=ix->x;exm[1]=ix->y;exm[2]=ix->z;


		std::vector<KDTree::point> pointSearch;
		std::vector<double> pointRadiusSquaredDistance;
		double radius=10;

		std::vector<std::pair<KDTree::_Examplar, double>> res_result;
		int num = kdtree.findNearest(exm, 10, res_result);//R半径搜索
		if(num==0)
			break;
		for(int i=0;i<num;++i)
		{
			KDTree::point _xyz;
			_xyz.x=res_result[i].first[0];
			_xyz.y=res_result[i].first[1];
			_xyz.z=res_result[i].first[2];
			pointSearch.push_back(_xyz);
			pointRadiusSquaredDistance.push_back(res_result[i].second);
		}
		_axyz.z=IIDV(pointSearch.begin(),pointSearch.end(),pointRadiusSquaredDistance.begin(),pointRadiusSquaredDistance.end());
		_cloud.push_back(_axyz);

		res_result.clear();
		pointSearch.clear();
		pointRadiusSquaredDistance.clear();
	}
	for(std::vector<KDTree::point>::iterator ix=_cloud.begin();ix!=_cloud.end();++ix)
		outfile1<<std::setprecision(10)<<ix->x+603000<<" "<<ix->y+3350000<<" "<<ix->z<<std::endl;

	/********************************二次曲面拟合**************************************/
	for(std::vector<KDTree::point>::iterator ix=_cloudInterpolation.begin();ix!=_cloudInterpolation.end();++ix)
	{
		KDTree::_Examplar exm(3);//待插值的点
		KDTree::point _axyz;//用来存储插值后的点云
		_axyz.x=ix->x;_axyz.y=ix->y;
		exm[0]=ix->x;exm[1]=ix->y;exm[2]=ix->z;
		std::vector<KDTree::point> pointSearch;
		double k=6;
		std::vector<std::pair<KDTree::_Examplar, double>> res_result;
		if(kdtree.find_K_Nearest(exm,k,res_result)==k)
		{
			for(int i=0;i<k;++i)
			{
				KDTree::point _xyz;
				_xyz.x=res_result[i].first[0];
				_xyz.y=res_result[i].first[1];
				_xyz.z=res_result[i].first[2];
				pointSearch.push_back(_xyz);
			}
			std::pair<double,double> coordinate(exm[0],exm[1]);
			_axyz.z=SFI(pointSearch.begin(),pointSearch.end(),coordinate);
			_cloud2.push_back(_axyz);
			res_result.clear();
			pointSearch.clear();
		}
	}
	for(std::vector<KDTree::point>::iterator ix=_cloud2.begin();ix!=_cloud2.end();++ix)
		outfile2<<std::setprecision(10)<<ix->x+603000<<" "<<ix->y+3350000<<" "<<ix->z<<std::endl;

	_cloudOriginal.clear();
	_cloudInterpolation.clear();
	_cloud.clear();
	_cloud2.clear();
	infile1.close();
	infile2.close();
	outfile1.close();
	outfile2.close();
	end = clock();
    dur = (double)(end - start);
	std::cout<<"总耗时："<<dur/CLOCKS_PER_SEC<<std::endl;
	return 0;
}
