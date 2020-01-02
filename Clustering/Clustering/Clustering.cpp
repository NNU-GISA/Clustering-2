
#include"Clustering.h"

DBSCAN::DBSCAN()
	: m_Points(NULL)
	, m_nNumPts(0)
	, m_nMinPts(0)
	, m_nRadius(0.0)
{
}


DBSCAN::DBSCAN(int minPts, double eps, std::vector<dPoint3D> points)
{
	this->m_nMinPts = minPts;
	this->m_nMinPts = points.size();
	this->m_nRadius = eps;

	if (this->m_nMinPts > 0)
	{
		this->m_Points = new dPoint3D[this->m_nMinPts];
		for (int i = 0; i < this->m_nNumPts; ++i)
		{
			this->m_Points[i] = points.at(i);
		}
	}
}


DBSCAN::~DBSCAN()
{
	if (m_Points) delete[]m_Points; m_Points = NULL;
}


int DBSCAN::initClustering(int minPts, double radius, std::vector<dPoint3D> points)
{
	this->m_nMinPts = minPts;
	this->m_nNumPts = points.size();
	this->m_nRadius = radius;

	if (this->m_nNumPts > 0)
	{
		this->m_Points = new dPoint3D[this->m_nNumPts];
		this->kdPts.pts.resize(points.size());

		for (int i = 0; i < this->m_nNumPts; ++i)
		{
			this->m_Points[i] = points.at(i);

			kdPts.pts[i].x = points.at(i).x;
			kdPts.pts[i].y = points.at(i).y;
			kdPts.pts[i].z = points.at(i).z;
		}

		this->m_KDTree = new KDTree(3, kdPts, KDTreeSingleIndexAdaptorParams(10));
		this->m_KDTree->buildIndex();
	}

	return 1;
}


void DBSCAN::clusterPts()
{
	int clusterID = 1;
	for (int i = 0; i < this->m_nNumPts; ++i)
	{
		if (this->m_Points[i].isParsed ) continue;
		if (this->m_Points[i].ptType == CORE) continue;

		this->obtainCluster(i, clusterID);

		this->m_Points[i].isParsed = true;

		clusterID++;
	}
}


void DBSCAN::obtainCluster(int ptID, int clusterID)
{
	//dPoint3D currPt = this->m_Points[ptID];

	double queryPt[3];
	queryPt[0] = this->m_Points[ptID].x;
	queryPt[1] = this->m_Points[ptID].y;
	queryPt[2] = this->m_Points[ptID].z;

	SearchParams params;
	std::vector<std::pair<size_t, double>> retMatches;
	int numMatches = this->m_KDTree->radiusSearch(&queryPt[0], m_nRadius, retMatches, params);

	if (numMatches < this->m_nMinPts)
	{
		this->m_Points[ptID].ptType = OUTLIER;
		this->m_Points[ptID].clusterID = 0;
	}
	else
	{
		this->m_Points[ptID].ptType = CORE;
		this->m_Points[ptID].clusterID = clusterID;

		std::stack<size_t> neighPts;
		
		for (int i = 1; i < numMatches; ++i)
		{
			int currPtID = retMatches[i].first;
			if (this->m_Points[currPtID].isPushed ) continue;
			if (this->m_Points[currPtID].isParsed) continue;
			//if (this->m_Points[currPtID].ptType == CORE) continue;

			neighPts.push(currPtID);
			this->m_Points[currPtID].isPushed= true;
			this->m_Points[currPtID].ptType = BORDER;
			this->m_Points[currPtID].clusterID = clusterID;
		}

		while (!neighPts.empty())
		{
			//size_t tempID = neighPts.top();
			
			double temp[3];
			temp[0] = this->m_Points[neighPts.top()].x;
			temp[1] = this->m_Points[neighPts.top()].y;
			temp[2] = this->m_Points[neighPts.top()].z;

			this->m_Points[neighPts.top()].isParsed = true;

			SearchParams tempParams;
			std::vector<std::pair<size_t, double>> tempMatches;
			int numTempMatches = this->m_KDTree->radiusSearch(&temp[0], m_nRadius, tempMatches, params);
			
			if (numTempMatches > m_nMinPts)
			{
				this->m_Points[neighPts.top()].ptType = CORE;
			}

			for (int m = 0; m < numTempMatches; ++m)
			{
				if (m_Points[tempMatches[m].first].isPushed) continue;
				//if (m_Points[tempMatches[m].first].ptType == CORE) continue;

				neighPts.push(tempMatches[m].first);
				this->m_Points[tempMatches[m].first].clusterID = clusterID;
				this->m_Points[tempMatches[m].first].isPushed = true;
				this->m_Points[tempMatches[m].first].ptType = BORDER;
			}

			neighPts.pop();

		}
	}
}



void DBSCAN::performClustering()
{
	this->assignPtType(this->m_nMinPts, this->m_nRadius);

	int clusterID = 1;
	for (int i = 0; i < this->m_nNumPts; ++i)
	{
		if (this->m_Points[i].isParsed) continue;
		if (this->m_Points[i].ptType != CORE) continue;

		this->getCluster(i, clusterID);

		this->m_Points[i].isParsed = true;
		clusterID++;
	}

	int k = 0;
}


void DBSCAN::getCluster(int ptID, int clusterID)
{
	dPoint3D currPt = this->m_Points[ptID];

	if (currPt.ptType == CORE)
	{
		this->m_Points[ptID].clusterID = clusterID;

		for (int i = 0; i < this->m_Points[ptID].retMatches.size(); ++i)
		{
			this->m_Points[this->m_Points[ptID].retMatches.at(i).first].clusterID = clusterID;
		}
	}
	else
	{
		this->m_Points[ptID].ptType = OUTLIER;
	}
}


void DBSCAN::assignPtType(int minPts, double radius)
{
	std::vector<std::pair<size_t, double>> retMatches;
	SearchParams params;
	double queryPt[3];
	for (int i = 0; i < this->m_nNumPts; ++i)
	{
		retMatches.clear();

		//dPoint3D currPt = this->m_Points[i];

		queryPt[0] = this->m_Points[i].x;
		queryPt[1] = this->m_Points[i].y;
		queryPt[2] = this->m_Points[i].z;

		int nMatches = this->m_KDTree->radiusSearch(&queryPt[0], radius, retMatches, params);

		if(nMatches >= minPts)
		{
			this->m_Points[i].ptType = CORE;
			this->m_Points[i].retMatches = retMatches;

			/*clock_t start = clock();
			for (int j = 0; j < nMatches; ++j)
			{
				size_t ptID = retMatches[j].first;
				dPoint3D tempPt = this->m_Points[ptID];
				if (tempPt.ptType == CORE) continue;
				this->m_Points[ptID].ptType = BORDER;
			}
			clock_t time = clock() - start;*/

			/*clock_t start1 = clock();
			for (int j = 0; j < nMatches; ++j)
			{
				if (m_Points[retMatches[j].first].ptType == CORE) continue;
				this->m_Points[retMatches[j].first].ptType = BORDER;
			}
			clock_t now1 = clock() - start1;*/

			int k = 0;
		}
		else
		{
			this->m_Points[i].ptType = OUTLIER;
		}
	}
}


OPTICS::OPTICS()
	: m_Points(NULL)
	, m_nRadius(0.0)
	, m_nMinPts(0)
	, m_nNumPts(0)
	, m_KDTree(NULL)
{

}

OPTICS::~OPTICS()
{

}

double OPTICS::Distance(dPoint3D ptA, dPoint3D ptB)
{
	return (sqrt(pow(ptA.x - ptB.x, 2) + pow(ptA.y - ptB.y, 2) + pow(ptA.z - ptB.z, 2)));
}


dPoint3D OPTICS::Centroid(dPoint3D* pts, int numPts)
{
	dPoint3D pt;

	for (int i = 0; i < numPts; ++i)
	{
		pt.x += pts[i].x;
		pt.y += pts[i].y;
		pt.z += pts[i].z;
	}

	pt.x /= numPts;
	pt.y /= numPts;
	pt.z /= numPts;

	return pt;
}



std::pair<dPoint3D, double> OPTICS::Region()
{
	std::pair<dPoint3D, double> currRegion;

	dPoint3D centroid = this->Centroid(this->m_Points, this->m_nNumPts);

	double maxDist = -(std::numeric_limits<double>::max)();

	for (int i = 0; i < this->m_nNumPts; ++i)
	{
		double currDist = this->Distance(centroid, this->m_Points[i]);
		if (currDist > maxDist)
		{
			maxDist + currDist;
		}
	}

	currRegion.first = centroid;
	currRegion.second = maxDist;
	return currRegion;
}



void OPTICS::initOptics(dPoint3D*& points, int numPts, int minPts, double radius)
{
	this->m_nNumPts = numPts;
	this->m_Points = new dPoint3D[this->m_nNumPts];
	
	this->kdPts.pts.resize(this->m_nNumPts);

	for (int i = 0; i < numPts; ++i)
	{
		this->m_Points[i] = points[i];
	
		this->kdPts.pts[i].x = points[i].x;
		this->kdPts.pts[i].y = points[i].y;
		this->kdPts.pts[i].z = points[i].z;
	}

	this->m_nMinPts = minPts;
	this->m_nRadius = radius;

	this->m_KDTree = new KDTree(3, kdPts, KDTreeSingleIndexAdaptorParams(10));
	this->m_KDTree->buildIndex();
}




void OPTICS::obtainCoreDist()
{
	SearchParams params;
	std::vector<std::pair<size_t, double>> retMatches;

	for (int i = 0; i < this->m_nNumPts; ++i)
	{
		if (m_Points[i].isParsed) continue;

		retMatches.clear();

		double queryPt[3];
		queryPt[0] = m_Points[i].x;
		queryPt[1] = m_Points[i].y;
		queryPt[2] = m_Points[i].z;

		size_t nMatches = this->m_KDTree->radiusSearch(&queryPt[0], m_nRadius, retMatches, params);

		//if()



	}



}