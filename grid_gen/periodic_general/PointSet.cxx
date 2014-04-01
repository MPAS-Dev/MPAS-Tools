#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include "PointSet.h"
#include "DensityFunction.h"

#define MIN(A,B) (B)<(A)?(B):(A)
#define MAX(A,B) (B)>(A)?(B):(A)

void voronoi_main(PointSet *);


PointSet::PointSet()
{
	nPoints = 0;
}


PointSet::~PointSet()
{

}


int PointSet::initFromTextFile(const char * filename)
{
	ifstream fin(filename);
	double xloc, yloc;
	Point * p;
	DensityFunction d;
	ifstream new_edges("new_edges");

	assert(fin.is_open());

	fin >> xloc >> yloc;
	do {
		p = new Point(xloc, yloc, 0);
		p->setNum(nPoints);
		nPoints++;
		points.push_back(p); 
		fin >> xloc >> yloc;
	} while (!fin.eof());
}


void PointSet::print()
{
	vector<Point*>::iterator it;

	cout << "We have " << nPoints << " points" << endl;

	for (it = points.begin(); it != points.end(); it++) {
		cout << **it << endl;
	}
}


void PointSet::printToTextFile(const char * filename)
{
	ofstream fout(filename);
	vector<Point*>::iterator it;

	assert(fout.is_open());

	for (it = points.begin(); it != points.end(); it++) {
		fout << (*it)->getX() << " " << (*it)->getY() << " " << (*it)->isBoundaryPoint() << endl;
	}
}


void PointSet::addPoint(double x, double y, int boundary_point)
{
	Point * p = new Point(x, y, boundary_point);
	p->setNum(nPoints);
	nPoints++;
	
	points.push_back(p);
}


void PointSet::addPoint(Point& p)
{
	Point * pp = new Point(p);
	nPoints++;
	
	points.push_back(pp);
}


int PointSet::size()
{
	return nPoints;
}


vector<Triangle>* PointSet::getTriangulation()
{
	triangulation = new vector<Triangle>; 

	voronoi_main(this);

	return triangulation;
}


vector<Point> * PointSet::getVoronoiDiagram()
{
	vector<Triangle> * t;
	vector<Triangle>::iterator it;
	vector<Point> * voronoiCorners = new vector<Point>[nPoints];
	Point p;
	int i, n;
	int nobtuse;

	// 1) Get a triangulation
	t = PointSet::getTriangulation();

	// 2) For each triangle, compute the associated Voronoi point
	//	 Add this point to the list of Voronoi corner for each of the triangle's vertices
	nobtuse = 0;
	for (it = triangulation->begin(); it != triangulation->end(); it++) {
		if (fabs(angle(it->getVertex(0), it->getVertex(1), it->getVertex(2))) > 3.1415926536/2.0) nobtuse++;
		if (fabs(angle(it->getVertex(1), it->getVertex(2), it->getVertex(0))) > 3.1415926536/2.0) nobtuse++;
		if (fabs(angle(it->getVertex(2), it->getVertex(0), it->getVertex(1))) > 3.1415926536/2.0) nobtuse++;
		p = it->circumcenter();
		for (i=0; i<3; i++) {
			n = it->getVertex(i).getNum();
			//assert(n >= 0 && n < nPoints);
			voronoiCorners[n].push_back(p);
		}
	}

cout << nobtuse << " obtuse angles\n";

	delete t;

	// 3) For each point, order its list of Voronoi corners in ccw order
	for (i=0; i<nPoints; i++)
		orderCCW(voronoiCorners[i], *points[i]);	

	// 4) Return list lists of Voronoi corners
	return voronoiCorners;
}


vector<Point> * PointSet::getDelaunayAdjacency()
{
	vector<Triangle> * t;
	vector<Triangle>::iterator it;
	vector<Point> * adjacencyList = new vector<Point>[nPoints];
	Point p0, p1, p2;
	int i, j, found, n0, n1, n2;

	t = PointSet::getTriangulation();

	for (it = triangulation->begin(); it != triangulation->end(); it++) {
		p0 = it->getVertex(0);
		p1 = it->getVertex(1);
		p2 = it->getVertex(2);
		
		n0 = p0.getNum(); 
		n1 = p1.getNum(); 
		n2 = p2.getNum(); 

		found = 0;
		for(j=0; j<adjacencyList[n0].size(); j++) {
			if (adjacencyList[n0][j].getNum() == n1)
				found = 1;
		}
		if (!found) adjacencyList[n0].push_back(p1);

		found = 0;
		for(j=0; j<adjacencyList[n0].size(); j++) {
			if (adjacencyList[n0][j].getNum() == n2)
				found = 1;
		}
		if (!found) adjacencyList[n0].push_back(p2);


		found = 0;
		for(j=0; j<adjacencyList[n1].size(); j++) {
			if (adjacencyList[n1][j].getNum() == n0)
				found = 1;
		}
		if (!found) adjacencyList[n1].push_back(p0);

		found = 0;
		for(j=0; j<adjacencyList[n1].size(); j++) {
			if (adjacencyList[n1][j].getNum() == n2)
				found = 1;
		}
		if (!found) adjacencyList[n1].push_back(p2);
		

		found = 0;
		for(j=0; j<adjacencyList[n2].size(); j++) {
			if (adjacencyList[n2][j].getNum() == n0)
				found = 1;
		}
		if (!found) adjacencyList[n2].push_back(p0);

		found = 0;
		for(j=0; j<adjacencyList[n2].size(); j++) {
			if (adjacencyList[n2][j].getNum() == n1)
				found = 1;
		}
		if (!found) adjacencyList[n2].push_back(p1);
	}

	// MGD Could replace the above later with a blind insertion followed by call to a routine remove_duplicates,
	//	  though this would involve removing items from the list

	delete t;

	for (i=0; i<nPoints; i++)
		orderCCW(adjacencyList[i], *points[i]);	

	return adjacencyList;
}


int PointSet::nearestPoint(Point& p)
{
	double d, minD;
	double x, y;
	int idx, minIdx;
	vector<Point*>::iterator it;

	x = p.getX();
	y = p.getY();

	minD = 1.e20;
	for (it = points.begin(), idx=0; it != points.end(); it++, idx++) {
		d = pow((*it)->getX() - x, 2.0) + pow((*it)->getY() - y, 2.0);
		if (d < minD) {minD = d; minIdx = idx;}
	}

	return minIdx;
}


Point* PointSet::operator[](int i)
{
	assert(i >= 0 && i < nPoints);
	return points[i];
}


double angle(Point o, Point p1, Point p2)
{
	double P1x, P1y, mP1;
	double P2x, P2y, mP2;
	double cos_angle;

	P1x = p1.getX() - o.getX();
	P1y = p1.getY() - o.getY();

	mP1 = sqrt(P1x*P1x + P1y*P1y);

	P2x = p2.getX() - o.getX();
	P2y = p2.getY() - o.getY();

	mP2 = sqrt(P2x*P2x + P2y*P2y);

	cos_angle = (P1x*P2x + P1y*P2y) / (mP1 * mP2);

	if (((P1x * P2y) - (P1y * P2x)) >= 0.0)
		return acos(MAX(MIN(cos_angle,1.0),-1.0));
	else
		return -acos(MAX(MIN(cos_angle,1.0),-1.0));

	return 1.0;
}


void orderCCW(vector<Point>& vc, Point p)
{
	int i, j;
	int vsize;
	double * angles;
	double ftemp;
	Point ptemp;

	const double PI = 3.14159265;

	vsize = vc.size();
	angles = new double[vsize];

	angles[0] = 0.0;
	for (i=1; i<vsize; i++) {
		angles[i] = angle(p, vc[0], vc[i]);
		if (angles[i] < 0.0) angles[i] += 2.0 * PI;
	}

	for(i=1; i<vsize; i++) {
		for(j=i+1; j<vsize; j++) {
			if (angles[j] < angles[i]) {
				ftemp = angles[i];
				angles[i] = angles[j];
				angles[j] = ftemp;

				ptemp = vc[i];
				vc[i] = vc[j];
				vc[j] = ptemp;
			}
		}
	}

	delete [] angles;
}


void orderCCW_print(vector<Point>& vc, Point p)
{
	int i, j;
	int vsize;
	double * angles;
	double ftemp;
	Point ptemp;

	const double PI = 3.14159265;

	vsize = vc.size();
	angles = new double[vsize];

	angles[0] = 0.0;
	for (i=1; i<vsize; i++) {
		angles[i] = angle(p, vc[0], vc[i]);
		if (angles[i] < 0.0) angles[i] += 2.0 * PI;
//cout << p << vc[0] << vc[i] << endl;
//cout << "angle=" << angles[i] << endl;
	}

	for(i=1; i<vsize; i++) {
		for(j=i+1; j<vsize; j++) {
			if (angles[j] < angles[i]) {
				ftemp = angles[i];
				angles[i] = angles[j];
				angles[j] = ftemp;

				ptemp = vc[i];
				vc[i] = vc[j];
				vc[j] = ptemp;
			}
		}
	}

	delete [] angles;
}
