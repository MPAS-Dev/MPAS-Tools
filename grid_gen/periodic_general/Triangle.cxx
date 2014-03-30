#include <math.h>
#include <assert.h>
#include "Triangle.h"

Triangle::Triangle()
{
	points[0] = Point(0.0, 0.0, 0);
	points[1] = Point(0.0, 0.0, 0);
	points[2] = Point(0.0, 0.0, 0);
}

Triangle::Triangle(Point a, Point b, Point c)
{
	points[0] = a;
	points[1] = b;
	points[2] = c;
}


Triangle::~Triangle()
{

}


void Triangle::setVertex(int i, Point p)
{
	assert(i >= 0 && i <= 2);
	points[i] = p;
}


Point Triangle::getVertex(int i)
{
	assert(i >= 0 && i <= 2);
	return points[i];
}


double Triangle::area()
{
	double a, b, c, s, R;

	// Compute side lengths
	a = sqrt(pow(points[0].getX() - points[1].getX(),2.0) + pow(points[0].getY() - points[1].getY(),2.0));
	b = sqrt(pow(points[1].getX() - points[2].getX(),2.0) + pow(points[1].getY() - points[2].getY(),2.0));
	c = sqrt(pow(points[0].getX() - points[2].getX(),2.0) + pow(points[0].getY() - points[2].getY(),2.0));

	// Compute semiperimiter
	s = (a + b + c) / 2.0;

	// Compute area
	return sqrt(s*(a + b - s)*(a + c - s)*(b + c - s));
}


Point Triangle::centroid()
{
	Point p;

	p.setX((points[0].getX() + points[1].getX() + points[2].getX()) * 0.33333333);
	p.setY((points[0].getY() + points[1].getY() + points[2].getY()) * 0.33333333);

	return p;
}


void Triangle::divide_segment(Point p1, Point p2, Point list[], int n)
{
	int i;
	Point vec;
  
	list[0] = p1;
	list[n-1] = p2;

	vec.setXY(p2.getX() - p1.getX(), p2.getY() - p1.getY());

	for(i=1; i<n-1; i++) {
		list[i] = p1 + vec * ((double)i / (double)(n-1));
	}
}


Point Triangle::centroid(DensityFunction& d, double * mass)
{
	//const int GLEV = 6;  // Also change 25 if changing this value
	const int GLEV = 17;  // Also change 256 if changing this value
	//const int GLEV = 19;  // Also change 324 if changing this value
	int i, j, k;
	double density, total_weight;
	Point o, c;
	Point line[GLEV][GLEV];
	Point p1p2[GLEV];
	Point p1p3[GLEV];
	Point p[256][3];
	Triangle t(o,o,o);	// Initially, we don't care what t is

	divide_segment(points[0], points[1], p1p2, GLEV);
	divide_segment(points[0], points[2], p1p3, GLEV);

	line[0][0] = points[0];
	line[1][0] = p1p2[1];
	line[1][1] = p1p3[1];

	for (i=2; i<GLEV; i++)
		divide_segment(p1p2[i], p1p3[i], line[i], i+1);

	k = 0;
	for(i=0; i<GLEV-1; i++) {
		for(j=0; j<=i; j++) {
			p[k][0] = line[i][j];
			p[k][1] = line[i+1][j];
			p[k][2] = line[i+1][j+1];
			k++;
		}
	}

	for(i=GLEV-1; i>=2; i--) {
		for(j=1; j<i; j++) {
			p[k][0] = line[i][j];
			p[k][1] = line[i-1][j];
			p[k][2] = line[i-1][j-1];
			k++;
		}
	}

	o.setXY(0.0, 0.0);
	total_weight = 0.0;
	for(i=0; i<256; i++) {
		t = Triangle(p[i][0], p[i][1], p[i][2]);
		c = t.centroid();
		density = d.evaluate(c);
		total_weight += density * t.area();
		o = o + (c * density * t.area());
	}
	o = o * (1.0 / total_weight);

	*mass = total_weight;

	return o;
}


Point Triangle::circumcenter()
{
	double x[3], y[3], xy2[3];
	double m[3][3];
	double a, bx, by;
	Point c;

	x[0] = points[0].getX();
	x[1] = points[1].getX();
	x[2] = points[2].getX();

	y[0] = points[0].getY();
	y[1] = points[1].getY();
	y[2] = points[2].getY();

	xy2[0] = x[0]*x[0] + y[0]*y[0];
	xy2[1] = x[1]*x[1] + y[1]*y[1];
	xy2[2] = x[2]*x[2] + y[2]*y[2];

	m[0][0] = x[0];
	m[0][1] = y[0];
	m[0][2] = 1.0;
	m[1][0] = x[1];
	m[1][1] = y[1];
	m[1][2] = 1.0;
	m[2][0] = x[2];
	m[2][1] = y[2];
	m[2][2] = 1.0;

	a = det(m);

	m[0][0] = xy2[0];
	m[0][1] = y[0];
	m[0][2] = 1.0;
	m[1][0] = xy2[1];
	m[1][1] = y[1];
	m[1][2] = 1.0;
	m[2][0] = xy2[2];
	m[2][1] = y[2];
	m[2][2] = 1.0;

	bx = -det(m);

	m[0][0] = xy2[0];
	m[0][1] = x[0];
	m[0][2] = 1.0;
	m[1][0] = xy2[1];
	m[1][1] = x[1];
	m[1][2] = 1.0;
	m[2][0] = xy2[2];
	m[2][1] = x[2];
	m[2][2] = 1.0;

	by = det(m);

	c.setXY(-bx/(2.0*a), -by/(2.0*a));
	c.setBoundaryPoint(0);

	return c; 
}


double Triangle::det(double m[3][3])
{
	return m[0][0] * (m[1][1]*m[2][2] - m[1][2]*m[2][1]) - m[0][1] * (m[1][0]*m[2][2] - m[1][2]*m[2][0]) + m[0][2] * (m[1][0]*m[2][1] - m[1][1]*m[2][0]);
}
