#include <cstdlib>
#include <math.h>
#include "DensityFunction.h"
#include "netcdf.h"
#include <valarray>

DensityFunction::DensityFunction(double X_PERIOD, double Y_PERIOD, int USE_DATA_DENSITY)
{

	minX = minY = 0.0;
	maxX = X_PERIOD;
	maxY = Y_PERIOD;
	use_data_density = USE_DATA_DENSITY;

	if (use_data_density == 1){
		read_density_netcdf(&xPosDG, &yPosDG, &densityDG, dxDG, dyDG);

		dxDG = xPosDG[1] - xPosDG[0];
		cout << "  dx=" << dxDG <<endl;

		dyDG = yPosDG[1] - yPosDG[0];
		cout << "  dy=" << dyDG <<endl;
	}
}


DensityFunction::~DensityFunction()
{

}


double DensityFunction::f(double x, double y)
{
// Densities should vary from 0 to 1.

	if (use_data_density == 1){
		return DataDensityFunction(x, y);
	} else {
		return AnalyticDensityFunction(x, y);
	}
}

double DensityFunction::AnalyticDensityFunction(double x, double y)
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



double DensityFunction::DataDensityFunction(double x, double y)
{
//	return UniformValue(x, y);
	return BilinearInterp(x, y);
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
	const double rhoMax = 1.0;	// maximum value of density function

	double x, y, U;
	double rrm;

	rrm = 1.0 / (double)RAND_MAX;

	x = (minX + (maxX - minX) * (double)rand() * rrm) ;
	y = (minY + (maxY - minY) * (double)rand() * rrm) ;
	// for general pdf, scale x and y to lie in domain of pdf
	U = (double)rand() * rrm;
        // The power of 0.5 is needed to get the proper density of points.
        // The logic is that grid spacing in 2d is approximately equal to
        // density to the 1/4 power.  However, cell area is approximately
        // grid spacing squared, leading to a 1/2 power.
	while(rhoMax * U >= pow(DensityFunction::f(x, y), 0.5)) {
		x = (minX + (maxX - minX) * (double)rand() * rrm) ;
		y = (minY + (maxY - minY) * (double)rand() * rrm) ;
		U = (double)rand() * rrm;
	}

	p.setX(x);
	p.setY(y);
}



/* ***** Read density function data ***** */
void DensityFunction::read_density_netcdf(double **xPosDG, double **yPosDG, double **densityDG, int dxDG, int dyDG)
{
	int ncerr;
	int ncid;
	int x_dimID, y_dimID;
	int x_dim, y_dim;
	int x_id, y_id, density_id;
	size_t temp;
        double max_density;

	cout << "Reading in data density function from density.nc" <<endl;

	// Open file
	ncerr = nc_open("density.nc", NC_NOWRITE, &ncid);
        if (ncerr>0) {
           cout << "Error reading density.nc.  Aborting." << endl;
           exit(1);
        }

	// Get needed dimensions
	ncerr = nc_inq_dimid(ncid, "x", &x_dimID);
	ncerr = nc_inq_dimlen(ncid, x_dimID, &temp);
	x_dim = (int)temp;
	ncerr = nc_inq_dimid(ncid, "y", &y_dimID);
	ncerr = nc_inq_dimlen(ncid, y_dimID, &temp);
	y_dim = (int)temp;
	cout << "  Got dimensions from file." <<endl;

	// Allocate variables to read into
	*xPosDG = (double *)malloc(sizeof(double) * (x_dim));
	*yPosDG = (double *)malloc(sizeof(double) * (y_dim));
	*densityDG = (double *)malloc(sizeof(double) * (x_dim) * (y_dim));
	cout << "  Allocated internal variables." <<endl;

	// Get needed variables
	ncerr = nc_inq_varid (ncid, "x", &x_id);
	ncerr = nc_get_var_double(ncid, x_id, *xPosDG);

	ncerr = nc_inq_varid (ncid, "y", &y_id);
	ncerr = nc_get_var_double(ncid, y_id, *yPosDG);

	ncerr = nc_inq_varid (ncid, "density", &density_id);
	ncerr = nc_get_var_double(ncid, density_id, *densityDG);
	cout << "  Got variables from file." <<endl;

        // Normalize data density to have a max of 1.0
        max_density = 0.0;
        for( int i = 0; i < x_dim*y_dim; i++ ) max_density = std::max(max_density, *(*densityDG+i));
        cout << "  Max density value is: " << max_density << ". Normalizing density field to have max of 1.0." << endl;
        for( int i = 0; i < x_dim*y_dim; i++ )  *(*densityDG+i) = *(*densityDG+i) / max_density;
        //max_density = 0.0;
        //for( int i = 0; i < x_dim*y_dim; i++ ) max_density = std::max(max_density, *(*densityDG+i));
        //cout << "  new Max density value is: " << max_density << endl;

	// Close file
	ncerr = nc_close(ncid);

	nxDG = x_dim;
	cout << "  nx=" << nxDG <<endl;

	nyDG = y_dim;
	cout << "  ny=" << nyDG <<endl;

	cout << "Completed reading data density." <<endl;

}


// Data density helpers

double DensityFunction::UniformValue(double x, double y)
{
// Gives the value of the density function cell that x,y falls inside of.

	int xpos, ypos; // the cells that the point falls in

	xpos = (int) floor( (x - xPosDG[0]) / dxDG);  // floor should not be needed since c++ will truncate...
	if (xpos < 0) {
		xpos = 0;
	} else if (xpos >= nxDG - 1) {
		xpos = nxDG - 2;
	}
	ypos = (int) floor( (y - yPosDG[0]) / dyDG);
	if (ypos < 0) {
		ypos = 0;
	} else if (ypos >= nyDG - 1) {
		ypos = nyDG - 2;
	}

	return densityDG[ypos * nxDG +  xpos];
}



double DensityFunction::BilinearInterp(double x, double y)
{
// Gives the value of the density function at x,y using Bilinear Interpolation

	int xpos, ypos; // the cells that the point falls in
	double value;

	xpos = (int) floor( (x - xPosDG[0]) / dxDG);  // floor should not be needed since c++ will truncate...
	if (xpos < 0) {
		xpos = 0;
	} else if (xpos >= nxDG - 1) {
		xpos = nxDG - 2;
	}
	ypos = (int) floor( (y - yPosDG[0]) / dyDG);
	if (ypos < 0) {
		ypos = 0;
	} else if (ypos >= nyDG - 1) {
		ypos = nyDG - 2;
	}

	value = (
		densityDG[ypos * nxDG +  xpos] * (xPosDG[xpos+1] - x) * (yPosDG[ypos+1] - y) +
		densityDG[(ypos+1) * nxDG +  xpos] * (xPosDG[xpos+1] - x) * (y - yPosDG[ypos]) +
		densityDG[ypos * nxDG +  xpos+1] * (x - xPosDG[xpos]) * (yPosDG[ypos+1] - y) +
		densityDG[(ypos+1) * nxDG +  xpos+1] * (x - xPosDG[xpos]) * (y - yPosDG[ypos])
		)	/ (dxDG * dyDG);

	return value;
}

