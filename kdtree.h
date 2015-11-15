#pragma once
#include <memory.h>
#include <stdio.h>
#include "assert.h"
#include <vector>

class TrainData
{

};

namespace KDTree
{
	struct point
	{
		float x;
		float y;
		float z;
	};
	struct _Examplar
	{
	public:

		_Examplar():dom_dims(0){}		//数据维度初始化为0
		_Examplar(const std::vector<float> elt, int dims)	//带有完整的两个参数的构造函数
		{													//这里const是为了保护原数据不被修改
			if(dims > 0)
			{
				dom_elt = elt;
				dom_dims = dims;
			}
			else
			{
				dom_dims = 0;
			}
		}
		_Examplar(int dims)	//只含有维度信息的constructor
		{
			if(dims > 0)
			{
				dom_elt.resize(dims);
				dom_dims = dims;
			}
			else
			{
				dom_dims = 0;
			}
		}
		_Examplar(const _Examplar& rhs)	//拷贝构造函数
		{
			if(rhs.dom_dims > 0)
			{
				dom_elt = rhs.dom_elt;
				dom_dims = rhs.dom_dims;
			}
			else
			{
				dom_dims = 0;
			}
		}
		_Examplar& operator=(const _Examplar& rhs)	//重载"="运算符
		{
			if(this == &rhs) 
				return *this;

			releaseExamplarMem();//清除现有数据向量，为读入数据准备

			if(rhs.dom_dims > 0)
			{
				dom_elt = rhs.dom_elt;
				dom_dims = rhs.dom_dims;
			}

			return *this;
		}
		~_Examplar()
		{
		}
		float& dataAt(int dim)	//定义访问控制函数
		{
			assert(dim < dom_dims);
			return dom_elt[dim];
		}
		float& operator[](int dim)	//重载"[]"运算符，实现下标访问
		{
			return dataAt(dim);//返回指针，从而可以改变数组内部的值
		}
		const float& dataAt(int dim) const	//定义只读访问函数
		{
			assert(dim < dom_dims);
			return dom_elt[dim];
		}
		const float& operator[](int dim) const	//重载"[]"运算符，实现下标只读访问
		{
			return dataAt(dim);
		}
		void create(int dims)	//创建数据向量
		{
			releaseExamplarMem();
			if(dims > 0)
			{
				dom_elt.resize(dims);	//控制数据向量维度
				dom_dims = dims;
			}
		}
		int getDomDims() const 	//获得数据向量维度信息
		{
			return dom_dims;
		}
		void set(std::vector<float>val)
		{
			if(dom_dims > 0)
			{
				for(int i=0;i<dom_dims;i++)
				{
					dom_elt[i] = val[i];
				}
			}			
		}
		void setTo(float val)	//数据向量初始化设置
		{
			if(dom_dims > 0)
			{
				for(int i=0;i<dom_dims;i++)
				{
					dom_elt[i] = val;
				}
			}
		}
	private:
		void releaseExamplarMem()	//清除现有数据向量
		{
			dom_elt.clear();
			dom_dims = 0;
		}
	private:
		std::vector<float> dom_elt;	//每个数据定义为一个float类型的向量
		int dom_dims;					//数据向量的维度
	};

	float Distance_exm(const _Examplar &x, const _Examplar &y);	//定义的距离函数

	class ExamplarCompare	//定义数据向量比较类，产生的对象用于sort的comp
	{
	public:
		ExamplarCompare(const int dim) : _dim(dim){}	//这里的dim是指待比较的方向
		bool
			operator()(const _Examplar &x, const _Examplar &y) const
		{
			return x[_dim] < y[_dim];
		}
	private:
		int _dim;	// don't make this const so that an assignment operator can be auto-generated
	};


	struct _HyperRectangle	//定义表示数据范围的超矩形结构
	{
		_Examplar min;		//统计数据集中所有数据向量每个维度上最小值组成的一个数据向量
		_Examplar max;		//统计数据集中所有数据向量每个维度上最大值组成的一个数据向量
		_HyperRectangle() {}
		_HyperRectangle(_Examplar mx, _Examplar mn)
		{
			assert (mx.getDomDims() == mn.getDomDims());
			min = mn;
			max = mx;
		}
		_HyperRectangle(const _HyperRectangle& rhs)
		{
			min = rhs.min;
			max = rhs.max;
		}
		_HyperRectangle& operator= (const _HyperRectangle& rhs)
		{
			if(this == &rhs)
				return *this;
			min = rhs.min;
			max = rhs.max;
			return *this;
		}
		void create(_Examplar mx, _Examplar mn)
		{
			assert (mx.getDomDims() == mn.getDomDims());
			min = mn;
			max = mx;
		}
	};

	class ExamplarSet : public TrainData	//整个数据集类，由一个抽象类TrainData派生
	{
	private:
		//_Examplar *_ex_set;
		std::vector<_Examplar> _ex_set;		//定义含有若干个_Examplar类数据向量的数据集
		size_t _size;							//数据集大小
		int _dims;							//数据集中每个数据向量的维度
	public:
		ExamplarSet():_size(0), _dims(0){}
		ExamplarSet(std::vector<_Examplar> ex_set, size_t size, int dims);
		ExamplarSet(size_t size, int dims);
		ExamplarSet(const ExamplarSet& rhs);
		ExamplarSet& operator=(const ExamplarSet& rhs);
		~ExamplarSet(){}

		_Examplar& examplarAt(int idx)
		{ 
			assert(idx < _size);
			return _ex_set[idx]; 
		}
		_Examplar& operator[](int idx)
		{
			return examplarAt(idx);
		}
		const _Examplar& examplarAt(int idx) const
		{
			assert(idx < _size);
			return _ex_set[idx];
		}
		void create(size_t size, int dims);
		void createfromcloud(std::vector<point>::iterator,std::vector<point>::iterator);
		int getDims() const { return _dims;}
		int getSize() const { return _size;}
		_HyperRectangle calculateRange();
		bool empty() const
		{
			return (_size == 0);
		}

		void sortByDim(int dim);	//按某个方向维的排序函数
		bool remove(int idx);		//去除数据集中排序后指定位置的数据向量
		void push_back(const _Examplar& ex)	//添加某个数据向量至数据集末尾
		{
			_ex_set.push_back(ex);
			_size++;
		}

		int readData(char *strFilePath);	//从文件读取数据集
	private:
		void releaseExamplarSetMem()		//清除现有数据集
		{
			_ex_set.clear();
			_size = 0;
		}
	};

	class KDTreeNode
	{
	private:
		int _split_dim;			 //该节点的最大区分度方向维
		_Examplar _dom_elt;		//该节点的数据向量
		_HyperRectangle _range_hr; //表示数据范围的超矩形结构,它表示的就是这一节点所代表的空间范围Range
	public:
		KDTreeNode *_left_child, *_right_child, *_parent;	//该节点的左右子树和父节点
	public:
		KDTreeNode():_left_child(0), _right_child(0), _parent(0), 
			_split_dim(0){}//构造函数

		KDTreeNode(KDTreeNode *left_child, KDTreeNode *right_child, 
			KDTreeNode *parent, int split_dim, _Examplar dom_elt, _HyperRectangle range_hr):
		_left_child(left_child), _right_child(right_child), _parent(parent),
			_split_dim(split_dim), _dom_elt(dom_elt), _range_hr(range_hr){}//带参数的构造函数

		KDTreeNode(const KDTreeNode &rhs);//拷贝构造函数

		KDTreeNode& operator=(const KDTreeNode &rhs);

		_Examplar& getDomElt() { return _dom_elt; }

		_HyperRectangle& getHyperRectangle(){ return _range_hr; }

		int& splitDim(){ return _split_dim; }

		void create(KDTreeNode *left_child, KDTreeNode *right_child, 
			KDTreeNode *parent, int split_dim, _Examplar dom_elt,  _HyperRectangle range_hr);
	};

	class KDTree	//k-d tree结构定义
	{
	public:
		KDTreeNode *_root;		//k-d tree的根节点
	public:
		KDTree():_root(NULL){}
		void create(const ExamplarSet &exm_set);		//创建k-d tree，实际上调用createKDTree
		void destroy();									//销毁k-d tree，实际上调用destroyKDTree
		~KDTree(){ destroyKDTree(_root); }
		std::pair<_Examplar, float> findNearest(_Examplar target);	//查找最近邻点函数，返回值是pair类型
		//实际是调用findNearest_i
		//查找距离在range范围内的近邻点，返回这样近邻点的个数，实际是调用findNearest_range
		int findNearest(_Examplar target, float range, std::vector<std::pair<_Examplar, float>> &res_nearest);

		//查找target最临近的k个点，返回查找到的k个点以及距离,实际调用findKNearest
		int find_K_Nearest(_Examplar target,int k,std::vector<std::pair<_Examplar, float>> &res_nearest);

	private:
		KDTreeNode* createKDTree(const ExamplarSet &exm_set);
		void destroyKDTree(KDTreeNode *root);
		std::pair<_Examplar, float> findNearest_i(KDTreeNode *root, _Examplar target);

		int findKNearest(KDTreeNode *root,_Examplar target,
			int k,std::vector<std::pair<_Examplar, float>> &res_nearest);

		int findNearest_range(KDTreeNode *root, _Examplar target, float range, 
			std::vector<std::pair<_Examplar, float>> &res_nearest);
	};

	class KDTreeCompare	//定义最近邻点向量比较类，产生的对象用于sort的comp
	{
	public:
		//		KDTreeCompare();	//这里的dim是指待比较的方向
		bool
			operator()(const std::pair<_Examplar, float> &x, const std::pair<_Examplar, float> &y) const
		{
			return x.second < y.second;
		}
		//private:
		//	int _dim;	// don't make this const so that an assignment operator can be auto-generated
	};
}
