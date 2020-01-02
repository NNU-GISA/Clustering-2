
#pragma once

#include<vector>

enum PtType
{
	CORE = 0,   //Core point
	BORDER = 1, //Border point
	OUTLIER = 2, //Outlier point
};


struct dPoint3D
{
	double x;
	double y;
	double z;

	double coreDist;
	double reachDist;

	unsigned int clusterID;

	bool isPushed;
	bool isVisited;
	bool isParsed;
	PtType ptType;
	std::vector<std::pair<size_t, double>> retMatches;
	

	dPoint3D()
	{
		x = y = z = 0.0;
		clusterID = 0;
		coreDist = 0.0;
		reachDist = 0.0;
		ptType = BORDER;
		isPushed = false;
		isVisited = false;
		isParsed = false;
	}
	dPoint3D& operator=(const dPoint3D* pt)
	{
		this->x = pt->x;
		this->y = pt->y;
		this->z = pt->z;
		this->clusterID = pt->clusterID;
		this->ptType = pt->ptType;
		this->isVisited = pt->isVisited;
		this->isParsed = pt->isParsed;
		this->retMatches = pt->retMatches;

		this->coreDist = pt->coreDist;
		this->reachDist = pt->reachDist;

		return *this;
	}
};