
#pragma once

#include<stack>
#include<time.h>
#include<math.h>
#include<limits>


#include"Header.h"
#include"nanoflann.hpp"


using namespace nanoflann;

typedef KDTreeSingleIndexAdaptor < L2_Simple_Adaptor<double, PointCloud2<double>>,
	PointCloud2<double>, 3 > KDTree;


class DBSCAN
{

public:
	DBSCAN();
	DBSCAN(int minPts, double eps, std::vector<dPoint3D> points);
	~DBSCAN();

	//-------------------------------- Variables --------------------------------
public:
	
	//Points;
	dPoint3D* m_Points;
	//Number of points;
	int m_nNumPts;
	//Minimum number of points;
	int m_nMinPts;
	//Radius[Epsilon] of the clustering;
	double m_nRadius;


private:

	KDTree* m_KDTree;
	PointCloud2<double> kdPts;

	//-------------------------------- Operations --------------------------------
public:

	int initClustering(int minPts, double radius, std::vector<dPoint3D> points);
	void clusterPts();

	void performClustering();

private:

	void obtainCluster(int ptID, int clusterID);

	void assignPtType(int minPts, double radius);
	void getCluster(int ptID, int clusterID);
};




class OPTICS
{
public:
	OPTICS();
	~OPTICS();

private:

	//-------------------------------- Variables --------------------------------
public:


	dPoint3D* m_Points;

	int m_nNumPts;
	double m_nRadius;
	int m_nMinPts;

private:
	PointCloud2<double> kdPts;
	KDTree* m_KDTree;



	//-------------------------------- Operations --------------------------------

public:

	void initOptics(dPoint3D*& points, int numPts, int minPts, double radius);

	double Distance(dPoint3D ptA, dPoint3D ptB);
	dPoint3D Centroid(dPoint3D* pts, int numPts);

	std::pair<dPoint3D, double> Region();

	void obtainCoreDist();
};
