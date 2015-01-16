#include <cstdlib>
#include <math.h>
#include "DensityFunction.h"


DensityFunction::DensityFunction(double X_PERIOD, double Y_PERIOD)
{
	minX = minY = 0.0;
	maxX = X_PERIOD;
	maxY = Y_PERIOD;
}


DensityFunction::~DensityFunction()
{

}


double DensityFunction::f(double x, double y)
{
	const double R1 = 0.03;	// inner Radius of bell
	const double R2 = 0.30;	// outer Radius of bell
	const double PI = 2.0 * acos(0.0);
	const double D1 = 0.25; // cell spacing inside bell
	const double D2 = 1.0; // cell spacing outside bell - should be 1.0

	double r;
	double x0 = 0.5*(minX + maxX);
	double y0 = 0.5*(minY + maxY);
	double xd, yd;
	double weight;

	xd = (x-x0) / (maxX-minY);
//	yd = (y-y0) / (maxY-minY);
	yd = (y-y0) / (maxX-minX);
	xd = fabs(xd);
	yd = fabs(yd);
	r = sqrt(xd*xd + yd*yd);

	if (r < R1)
		return pow(D1, 4.0);
	else if (r >= R1 && r < R2) {
		weight = cos(PI/2.0*(r-R1)/(R2-R1));
		return pow(weight*D1 + (1.0-weight)*D2, 4.0);
//		return pow(1.0 + cos(PI/2.0*(r-R1)/(R2-R1)),4.0);
		}
	else
		return D2;
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
