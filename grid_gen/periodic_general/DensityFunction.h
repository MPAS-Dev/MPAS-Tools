#ifndef _DensityFunctionH
#define _DensityFunctionH
#include "Point.h"

using namespace std;

class DensityFunction
{
	private:
		double minX, maxX, minY, maxY;
		double f(double x, double y);
	public:
		DensityFunction();
		~DensityFunction();
		double evaluate(Point& p);
		Point * randomPoint();
		void randomPoint(Point& p);
};
#endif
