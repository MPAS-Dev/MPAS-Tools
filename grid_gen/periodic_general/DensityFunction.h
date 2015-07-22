#ifndef _DensityFunctionH
#define _DensityFunctionH
#include "Point.h"

using namespace std;

class DensityFunction
{
	private:
		double minX, maxX, minY, maxY;
		double f(double x, double y);
		double AnalyticDensityFunction(double x, double y);
		double DataDensityFunction(double x, double y);
		double *xPosDG, *yPosDG, *densityDG;  // The x (1d), y (1d), and density (2d) values of the data density function (regular grid)
		int dxDG, dyDG; // grid spacing on the regular data density grid
		int nxDG, nyDG; // number of cells on regular data density grid
		void read_density_netcdf(double **xPosDG, double **yPosDG, double **densityDG, int dxDG, int dyDG);
		double UniformValue(double x, double y);
		int use_data_density;
	public:
		DensityFunction(double X_PERIOD, double Y_PERIOD, int USE_DATA_DENSITY);
		~DensityFunction();
		double evaluate(Point& p);
		Point * randomPoint();
		void randomPoint(Point& p);
};


#endif
