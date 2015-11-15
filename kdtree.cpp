#include <algorithm>
#include <fstream>
#include "kdtree.h"
#include <limits> 

KDTree::ExamplarSet::ExamplarSet(std::vector<_Examplar> ex_set, size_t size, int dims)
{
	if(size > 0)
		_size = size;
	else
		_size = 0;
	if(dims > 0)
		_dims = dims;
	else
		_dims = 0;

	_ex_set = ex_set;
}


KDTree::ExamplarSet::ExamplarSet( size_t size, int dims )
{
	if(size > 0)
		_size = size;
	else
		_size = 0;
	if(dims > 0)
		_dims = dims;
	else
		_dims = 0;
}

KDTree::ExamplarSet::ExamplarSet( const ExamplarSet& rhs )
{
	if(rhs._size > 0)
		_size = rhs._size;
	else
		_size = 0;
	if(rhs._dims > 0)
		_dims = rhs._dims;
	else
		_dims = 0;

	_ex_set = rhs._ex_set;
}

KDTree::ExamplarSet& KDTree::ExamplarSet::operator=( const ExamplarSet& rhs )
{

	if(this == &rhs) 
		return *this;

	releaseExamplarSetMem();

	if(rhs._size > 0)
		_size = rhs._size;
	else
		_size = 0;
	if(rhs._dims > 0)
		_dims = rhs._dims;
	else
		_dims = 0;

	_ex_set = rhs._ex_set;

	return *this;
}

void KDTree::ExamplarSet::create(size_t size, int dims )
{
	releaseExamplarSetMem();
	if(size > 0 && dims > 0)
	{
		_ex_set.resize(size);
		_size = size;
		_dims = dims;
		for(int i=0;i<_size;i++)
		{
			_ex_set[i].create(_dims);
		}
	}
}
void KDTree::ExamplarSet::createfromcloud(std::vector<point>::iterator begin,std::vector<point>::iterator end)
{
	size_t n=end-begin;
	create(n,3);
	size_t i=0;
	for(std::vector<point>::iterator iter=begin;iter!=end;++iter)
	{
		std::vector<float>val;
		val.push_back(iter->x);val.push_back(iter->y);val.push_back(iter->z);
		_ex_set[i++].set(val);
	}
}

KDTree::_HyperRectangle KDTree::ExamplarSet::calculateRange()
{
	assert(_size > 0);
	assert(_dims > 0);
	_Examplar mn(_dims);
	_Examplar mx(_dims);

	for(int j=0;j<_dims;j++)
	{
		mn.dataAt(j) = (*this)[0][j];	//初始化最小范围向量
		mx.dataAt(j) = (*this)[0][j];	//初始化最大范围向量
	}

	for(int i=1;i<_size;i++)
	{
		for(int j=0;j<_dims;j++)
		{
			if( (*this)[i][j] < mn[j] )
				mn[j] = (*this)[i][j];
			if( (*this)[i][j] > mx[j] )
				mx[j] = (*this)[i][j];
		}
	}
	_HyperRectangle hr(mx, mn);

	return hr;
}

void KDTree::ExamplarSet::sortByDim( int dim )
{
	ExamplarCompare cmp(dim);
	std::sort(_ex_set.begin(), _ex_set.end(), cmp);
}

bool KDTree::ExamplarSet::remove( int idx )
{
	if(idx >=0 && idx < _size)
	{
		_ex_set.erase(_ex_set.begin() + idx);
		_size --;
		return true;
	}
	else
		return false;
}

int KDTree::ExamplarSet::readData( char *strFilePath )
{
	std::ifstream fin(strFilePath);
	assert(fin != NULL);

	float temp;

	int row_id = 0, column_id = 0;;
	while(!fin.eof())
	{
		//! 获得一个数据样本

		for(column_id = 0;column_id < _dims;column_id++)
		{
			fin >> temp;
			(*this)[row_id][column_id] = temp;
		}

		row_id ++;
		if(row_id == _size)
			break;
	}

	fin.close();
	return 0;
}

float KDTree::Distance_exm( const _Examplar &x, const _Examplar &y )//计算多为空间向量的欧式空间距离
{
	float dis;
	if(x.getDomDims() == y.getDomDims() && x.getDomDims() > 0)
	{
		dis = 0.0;
		for(int i=0;i<x.getDomDims();i++)
		{
			dis += (x[i] - y[i]) * (x[i] - y[i]);
		}
		dis = sqrt(dis);
	}
	else
		dis = -1.0;
	return dis;
}

void KDTree::KDTreeNode::create( KDTreeNode *left_child, KDTreeNode *right_child,
	KDTreeNode *parent, int split_dim, _Examplar dom_elt,  _HyperRectangle range_hr)
{
	this->_left_child = left_child;
	this->_right_child = right_child;
	this->_parent = parent;
	this->_split_dim = split_dim;
	this->_dom_elt = dom_elt;
	this->_range_hr = range_hr;
}

KDTree::KDTreeNode::KDTreeNode( const KDTreeNode &rhs )
{
	this->_left_child = rhs._left_child;
	this->_right_child = rhs._right_child;
	this->_parent = rhs._parent;
	this->_split_dim = rhs._split_dim;
	this->_dom_elt = rhs._dom_elt;
	this->_range_hr = rhs._range_hr;
}

KDTree::KDTreeNode& KDTree::KDTreeNode::operator=( const KDTreeNode &rhs )
{
	if(this == &rhs) 
		return *this;
	this->_left_child = rhs._left_child;
	this->_right_child = rhs._right_child;
	this->_parent = rhs._parent;
	this->_split_dim = rhs._split_dim;
	this->_dom_elt = rhs._dom_elt;
	this->_range_hr = rhs._range_hr;

	return *this;
}


void KDTree::KDTree::create( const ExamplarSet &exm_set )
{
	_root = createKDTree(exm_set);
}

KDTree::KDTreeNode* KDTree::KDTree::createKDTree( const ExamplarSet &exm_set )
{
	if(exm_set.empty())
		return NULL;
	//当前创建kd-tree的数据集
	ExamplarSet exm_set_copy(exm_set);


	int dims = exm_set_copy.getDims();
	int size = exm_set_copy.getSize();

	//! 计算每个维的方差，选出方差值最大的维
	float var_max = -0.1f; //最大方差
	float avg, var;
	int dim_max_var = -1;//方差最大的维
	for(int i=0;i<dims;i++)
	{
		avg = 0;
		var = 0;
		//! 求某一维的总和
		for(int j=0;j<size;j++)
		{
			avg += exm_set_copy[j][i];
		}
		//! 求平均
		avg /= size;
		//! 求方差
		for(int j=0;j<size;j++)
		{
			var += ( exm_set_copy[j][i] - avg ) * 
				( exm_set_copy[j][i] - avg );
		}
		var /= size;
		if(var > var_max)
		{
			var_max = var;
			dim_max_var = i;
		}
	}

	//! 确定节点的数据矢量

	_HyperRectangle hr = exm_set_copy.calculateRange();  //统计节点空间范围，获得最大包含范围
	exm_set_copy.sortByDim(dim_max_var);				 //将所有数据向量按最大区分度方向排序调用标准库函数sort()
	int mid = size / 2;
	_Examplar exm_split = exm_set_copy.examplarAt(mid);//取出排序结果的中间节点
	exm_set_copy.remove(mid);							//将中间节点作为父（根）节点，所以将其从数据集中去除

	//! 确定左右节点

	ExamplarSet exm_set_left(0, exm_set_copy.getDims());
	ExamplarSet exm_set_right(0, exm_set_copy.getDims());
	exm_set_right.remove(0);//调用vector容器删除指定元素函数

	int size_new = exm_set_copy.getSize();
	for(int i=0;i<size_new;i++)
	{
		_Examplar temp = exm_set_copy[i];
		if( temp.dataAt(dim_max_var) < 
			exm_split.dataAt(dim_max_var) )
			exm_set_left.push_back(temp);
		else
			exm_set_right.push_back(temp);
	}

	KDTreeNode *pNewNode = new KDTreeNode(0, 0, 0, dim_max_var, exm_split, hr);//创建该节点的最大最小向量的范围

	pNewNode->_left_child = createKDTree(exm_set_left);
	if(pNewNode->_left_child != NULL)
		pNewNode->_left_child->_parent = pNewNode;//指向自己

	pNewNode->_right_child = createKDTree(exm_set_right);
	if(pNewNode->_right_child != NULL)
		pNewNode->_right_child->_parent = pNewNode;

	return pNewNode;
}

void KDTree::KDTree::destroyKDTree( KDTreeNode *root )
{
	if(root != NULL)
	{
		destroyKDTree(root->_left_child);
		destroyKDTree(root->_right_child);
		delete root;
	}
}

void KDTree::KDTree::destroy()
{
	destroyKDTree(_root);
}

std::pair<KDTree::_Examplar, float> KDTree::KDTree::findNearest_i( KDTreeNode *root, _Examplar target )
{
	//! 向下到达叶子节点

	KDTreeNode *pSearch = root;

	//! 堆栈用于保存搜索路径
	std::vector<KDTreeNode*> search_path;

	_Examplar nearest;

	float max_dist;

	while(pSearch != NULL)//向下到达叶子节点
	{
		search_path.push_back(pSearch);
		int s = pSearch->splitDim();
		if(target[s] <= pSearch->getDomElt()[s])
		{
			pSearch = pSearch->_left_child;
		}
		else
		{
			pSearch = pSearch->_right_child;
		}
	}

	nearest = search_path.back()->getDomElt();//返回容器最后一个元素的引用
	max_dist = Distance_exm(nearest, target);//初始最近点到搜索点的距离

	search_path.pop_back();//删除返回后的元素

	//! 回溯搜索路径
	while(!search_path.empty())
	{
		KDTreeNode *pBack = search_path.back();//再返回一个元素
		search_path.pop_back();

		if( pBack->_left_child == NULL && pBack->_right_child == NULL)//判断到达叶子节点
		{
			if( Distance_exm(nearest, target) > Distance_exm(pBack->getDomElt(), target) )
			{
				nearest = pBack->getDomElt();
				max_dist = Distance_exm(pBack->getDomElt(), target);
			}
		}
		else
		{
			int s = pBack->splitDim();
			if( abs(pBack->getDomElt()[s] - target[s]) < max_dist)
			{
				if( Distance_exm(nearest, target) > Distance_exm(pBack->getDomElt(), target) )
				{
					nearest = pBack->getDomElt();
					max_dist = Distance_exm(pBack->getDomElt(), target);
				}
				if(target[s] <= pBack->getDomElt()[s])
					pSearch = pBack->_right_child;
				else
					pSearch = pBack->_left_child;
				if(pSearch != NULL)
					search_path.push_back(pSearch);
			}
		}
	}

	std::pair<_Examplar, float> res(nearest, max_dist);

	return res;
}

std::pair<KDTree::_Examplar, float> KDTree::KDTree::findNearest( _Examplar target )
{
	std::pair<_Examplar, float> res;
	if(_root == NULL)
	{
		res.second = std::numeric_limits<float>::infinity();
		return res;
	}
	return findNearest_i(_root, target);
}

int KDTree::KDTree::findNearest( _Examplar target, float range, std::vector<std::pair<_Examplar, float>> &res_nearest )
{
	return findNearest_range(_root, target, range, res_nearest);
}


int KDTree::KDTree::findNearest_range( KDTreeNode *root, _Examplar target, float range, 
	std::vector<std::pair<_Examplar, float>> &res_nearest )
{
	if(root == NULL)
		return 0;
	float dist_sq, dx;
	int ret, added_res = 0;
	dist_sq = 0;
	dist_sq = Distance_exm(root->getDomElt(), target);

	if(dist_sq <= range) {
		std::pair<_Examplar,float> temp(root->getDomElt(), dist_sq);
		res_nearest.push_back(temp);

		//! 结果个数+1

		added_res = 1;
	}

	dx = target[root->splitDim()] - root->getDomElt()[root->splitDim()];
	//! 左子树或右子树递归的查找
	ret = findNearest_range(dx <= 0.0 ? root->_left_child : root->_right_child, target, range, res_nearest);
	if(ret >= 0 && fabs(dx) < range) {
		added_res += ret;
		ret = findNearest_range(dx <= 0.0 ? root->_right_child : root->_left_child, target, range, res_nearest);
	}

	added_res += ret;
	return added_res;
}

int KDTree::KDTree::find_K_Nearest(_Examplar target,int k,std::vector<std::pair<_Examplar, float>> &res_nearest)
{
	return findKNearest(_root, target, k, res_nearest);
}

int KDTree::KDTree::findKNearest(KDTreeNode *root,_Examplar target,
	int k,std::vector<std::pair<_Examplar, float>> &res_nearest)
{
	//! 向下到达叶子节点

	KDTreeNode *pSearch = root;

	//! 堆栈用于保存搜索路径
	std::vector<KDTreeNode*> search_path;

	_Examplar nearest;

	float max_dist;

	while(pSearch != NULL)//向下到达叶子节点
	{
		search_path.push_back(pSearch);
		int s = pSearch->splitDim();
		if(target[s] <= pSearch->getDomElt()[s])
		{
			pSearch = pSearch->_left_child;
		}
		else
		{
			pSearch = pSearch->_right_child;
		}
	}

	nearest = search_path.back()->getDomElt();//返回容器最后一个元素的引用
	max_dist = Distance_exm(nearest, target);//初始最近点到搜索点的距离
	search_path.pop_back();//删除返回后的元素

	//将搜索路径中的所有元素及到搜索点间的距离都保存下来

	std::pair<_Examplar, float> res1(nearest, max_dist);
	res_nearest.push_back(res1);//将搜索路径中最近邻的

	//! 回溯搜索路径
	while(!search_path.empty())
	{
		KDTreeNode *pBack = search_path.back();//再返回一个元素
		search_path.pop_back();

		//将搜索路径中的所有元素及到搜索点间的距离都保存下来
		float _dist=Distance_exm(pBack->getDomElt(), target);
		std::pair<_Examplar, float> _res(pBack->getDomElt(), _dist);
		res_nearest.push_back(_res);

		if( pBack->_left_child == NULL && pBack->_right_child == NULL)//判断到达叶子节点
		{
			if( Distance_exm(nearest, target) > Distance_exm(pBack->getDomElt(), target) )
			{
				nearest = pBack->getDomElt();
				max_dist = Distance_exm(pBack->getDomElt(), target);
			}
		}
		else
		{
			int s = pBack->splitDim();
			if( abs(pBack->getDomElt()[s] - target[s]) < max_dist)
			{
				if( Distance_exm(nearest, target) > Distance_exm(pBack->getDomElt(), target) )
				{
					nearest = pBack->getDomElt();
					max_dist = Distance_exm(pBack->getDomElt(), target);
				}
				if(target[s] <= pBack->getDomElt()[s])
					pSearch = pBack->_right_child;
				else
					pSearch = pBack->_left_child;
				if(pSearch != NULL)
					search_path.push_back(pSearch);
			}
		}
	}
	KDTreeCompare	cmp;
	std::sort(res_nearest.begin(), res_nearest.end(),cmp);//对数组元素按小到大排序，去前面的k位
	if(k > res_nearest.size())
		return 0;
	res_nearest.erase(res_nearest.begin()+k,res_nearest.end());//删除
	return res_nearest.size();
}
