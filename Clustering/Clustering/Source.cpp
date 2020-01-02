

#include<iostream>
#include<stdio.h>
#include<vector>


#include"Clustering.h"
#include<time.h>



void main()
{

	const int minPts = 4;
	//const double radius = 0.75*0.75;
	const double radius = 0.2*0.2;


	std::vector<dPoint3D> points;
	
	/*FILE* inFile = fopen("Points.xyz", "r");
	dPoint3D pt;
	int k = 0;
	while (!feof(inFile))
	{
		fscanf(inFile, "%lf,%lf,%lf,%d\n", &pt.x, &pt.y, &pt.z, &k);
		points.push_back(pt);
	}*/

	FILE *inFile = fopen("gauss.xyz", "r");
	dPoint3D pt;
	double k = 0;
	while (!feof(inFile))
	{
		fscanf(inFile, "%lf	%lf	%lf\n", &pt.x, &pt.y, &k);
		points.push_back(pt);
	}


	DBSCAN* cluster = new DBSCAN;

	cluster->initClustering(minPts, radius, points);
//	cluster->runClustering();
	clock_t time = clock();
	cluster->clusterPts();
	//cluster->performClustering();
	clock_t endTime = clock() - time;
	std::cout << "Time consts:  " << endTime << std::endl;
	FILE* outFile = fopen("outPointsGauss.xyz", "w");
	for (int i = 0; i < cluster->m_nNumPts; ++i)
	{
		fprintf(outFile, "%lf %lf %lf %d %d\n",
			cluster->m_Points[i].x,
			cluster->m_Points[i].y,
			cluster->m_Points[i].z,
			cluster->m_Points[i].clusterID,
			cluster->m_Points[i].ptType);
	}
	fclose(outFile);

}

