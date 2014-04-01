#ifndef _TriangleH
#define _TriangleH
#include <iostream>
#include "Point.h"
#include "DensityFunction.h"
using namespace std;

class Triangle
{
	private:
		Point points[3];
		double det(double m[3][3]);
		void divide_segment(Point p1, Point p2, Point list[], int n);

	public:
		Triangle();
		Triangle(Point a, Point b, Point c);
		~Triangle();
		void setVertex(int i, Point p);
		Point getVertex(int i) const;
		double area();
		Point centroid();
		Point centroid(DensityFunction& d, double * mass);
		Point circumcenter();
		friend bool operator==(Triangle& lhs, Triangle& rhs);
		friend bool operator<(Triangle const& lhs, Triangle const& rhs);
};

#endif
