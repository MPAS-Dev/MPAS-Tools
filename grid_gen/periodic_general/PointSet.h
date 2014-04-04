#ifndef _PointSetH
#define _PointSetH
#include <vector>
#include <list>
#include "Point.h"
#include "Triangle.h"
using namespace std;

typedef struct tagFreenode
{
	struct tagFreenode * nextfree;
} Freenode ;


typedef struct tagFreelist
{
	Freenode * head;
	int nodesize;
} Freelist ;

typedef struct tagPoint
{
	double x ;
	double y ;
} VPoint ;

/* structure used both for sites and for vertices */

typedef struct tagSite
{
	VPoint coord ;
	int sitenbr ;
	int refcnt ;
} Site ;


typedef struct tagEdge
{
	double a, b, c ;
	Site * ep[2] ;
	Site * reg[2] ;
	int edgenbr ;
} Edge ;

class PointSet
{
	private:
		int nPoints;
		vector<Point*> points;
		vector<Triangle> * triangulation;

	public:
		PointSet();
		~PointSet();
		int initFromTextFile(const char *);
		void print();
		void printToTextFile(const char *);
		void addPoint(double x, double y, int boundary_point);
		void addPoint(Point& p);
		int size();
		vector<Triangle>* getTriangulation();
		vector<Point> * getVoronoiDiagram();
		vector<Point> * getDelaunayAdjacency();
		int nearestPoint(Point& p);
		Point* operator[](int i);
		friend void readsites(PointSet * p);
		friend void out_triple(PointSet * p, Site * s1, Site * s2, Site * s3);
};

double angle(Point o, Point p1, Point p2);
void orderCCW(vector<Point>& vc, Point p);
void orderCCW_normalize(vector<Point>& vc, Point p, double x_period, double y_period);
void orderCCW_normalize2(vector<Point>& vc1, vector<Point>& vc2, Point p, double x_period, double y_period);
double poly_area(vector<Point>& vc);
void orderCCW_print(vector<Point>& vc, Point p);
void periodic_normalize(vector<Point>& vc, double x_period, double y_period);
#endif
