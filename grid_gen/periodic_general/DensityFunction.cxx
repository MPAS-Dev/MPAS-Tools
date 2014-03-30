#include <cstdlib>
#include <math.h>
#include "DensityFunction.h"


DensityFunction::DensityFunction()
{
	minX = minY = 0.0;
	maxX = 40.0;
	maxY = 40.0 * 0.866025403784439;
}


DensityFunction::~DensityFunction()
{

}


double DensityFunction::f(double x, double y)
{
	const double R = 0.125;	// Radius of bell
	const double R1 = 0.0875;	// Radius of bell
	const double R2 = 0.1750;	// Radius of bell
	const double PI = 3.141592653;

	double r;
	double x0 = 0.5*(minX + maxX);
	double y0 = 0.5*(minY + maxY);
	double xd, yd;

	xd = (x-x0) / (maxX-minY);
//	yd = (y-y0) / (maxY-minY);
	yd = (y-y0) / (maxX-minX);
	xd = fabs(xd);
	yd = fabs(yd);
	r = sqrt(xd*xd + yd*yd);

	if (r < R1)
		return 81.0;
	else if (r >= R1 && r < R2)
		return pow(1.0 + 1.0*(1.0 + cos(PI*(r-R1)/(R2-R1))),4.0);
	else
		return 1.0;
}


double DensityFunction::evaluate(Point& p)
{
	return f(p.getX(), p.getY());
}


Point * DensityFunction::randomPoint()
{
	Point * p;

	p = new Point();
	DensityFunction::randomPoint(*p);
	return p;
}


void DensityFunction::randomPoint(Point& p)
{
	const double rhoMax = 3.0;	// maximum value of density function

	double x, y, U;
	double rrm;

	rrm = 1.0 / (double)RAND_MAX;

	x = (double)rand() * rrm;
	y = (double)rand() * rrm;
	// for general pdf, scale x and y to lie in domain of pdf
	U = (double)rand() * rrm;
	while(rhoMax * U >= DensityFunction::f(x, y)) {
		x = (double)rand() * rrm;
		y = (double)rand() * rrm;
		U = (double)rand() * rrm;
	}

	p.setX(minX + (maxX - minX)*x);
	p.setY(minY + (maxY - minY)*y);
}
