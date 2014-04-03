#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <math.h>
#include "PointSet.h"
#include "Triangle.h"
#include "DensityFunction.h"
#include "netcdf.h"
using namespace std;

#define EPS 1.0e-7

#define ALLOC_INT2D(ARR,I,J) (ARR) = new int*[(I)]; for(int i=0; i<(I); i++) (ARR)[i] = new int[(J)];
#define DEALLOC_INT2D(ARR,I,J) for(int i=0; i<(I); i++) delete [] (ARR)[i]; delete [] (ARR);

#define ALLOC_REAL2D(ARR,I,J) (ARR) = new double*[(I)]; for(int i=0; i<(I); i++) (ARR)[i] = new double[(J)];
#define DEALLOC_REAL2D(ARR,I,J) for(int i=0; i<(I); i++) delete [] (ARR)[i]; delete [] (ARR);

void read_netcdf(int *nCells, int *nVertices, int *vertexDegree, 
		double **xCell, double **yCell, double **zCell, 
		double **xVertex, double **yVertex, double **zVertex, 
		double **meshDensity, int **cellsOnVertex,
		double *x_period, double *y_period);

Point segment_intersect(Point& p0, Point &p1, Point &q0, Point&q1);

int main(int argc, char ** argv)
{
//	int i, ii, jj, n, iter, idx, npts, np;
//	DensityFunction f;
//	PointSet pset;
//	PointSet out_pset;
//	Point * cells;
//	Point * temp_p;
//	Point * temp_pp;
//	Point p3;
//	Triangle t;
//	Point p, p2;
//	vector<Point> * clist;
//	vector<Triangle> * triangulation;
//	vector<Triangle>::iterator it;
//	set<Triangle> delaunay_tri;
//	set<Triangle>::iterator dti;
//	list<Triangle> norm_dt;
//	list<Triangle>::iterator norm_dti;
//	vector< vector<Point> > vertices_on_cell;
//	vector< vector<Point> > cells_on_cell;
//	vector< set<Point> > coc;
//	set<Point>::iterator cell_iter;
//	vector< vector<Point> > cv_on_cell;
//	Triangle * tri;
//	double xcell, ycell;
//	double x, y;
//	double total_mass, mass; 
//	FILE * restart;
	int nCells, nVertices, vertexDegree;
	double *xCell, *yCell, *zCell, *xVertex, *yVertex, *zVertex, *meshDensity;
	int *cellsOnVertex;
	double x_period, y_period;


	/*
	 * Read basic grid info from NetCDF file
	 */
	read_netcdf(&nCells, &nVertices, &vertexDegree, &xCell, &yCell, &zCell, &xVertex, &yVertex, &zVertex, &meshDensity, &cellsOnVertex, &x_period, &y_period);

	cout << "Read from input file:" << endl;
	cout << "   nCells    = " << nCells << endl;
	cout << "   nVertices = " << nVertices << endl;
	cout << "   vertexDegree = " << vertexDegree << endl;
	cout << "   x_period = " << x_period << endl;
	cout << "   y_period = " << y_period << endl;
	cout << endl;

#if 0
	vertices_on_cell.resize(nCells);
	for (i=0; i<nVertices; i++) {
		temp_p = new Point(xVertex[i], yVertex[i], 0);
		temp_p->setNum(i);
		for (int j=0; j<3; j++) {
			vertices_on_cell[cellsOnVertex[3*i+j]-1].push_back(*temp_p);
		}
	}

	for (i=0; i<nCells; i++)
	{
		p.setX(xCell[i]);
		p.setY(yCell[i]);
		orderCCW_normalize(vertices_on_cell[i], p, x_period, y_period);
/*
		for (int j=0; j<vertices_on_cell[i].size(); j++)
			cout << vertices_on_cell[i][j] << endl;
		cout << vertices_on_cell[i][0] << endl;
		cout << endl;
*/
	}

	coc.resize(nCells);
	for (i=0; i<nVertices; i++) {

		/* First cell on vertex */
		temp_p = new Point(xCell[cellsOnVertex[3*i]-1], yCell[cellsOnVertex[3*i]-1], 0);
		temp_p->setNum(cellsOnVertex[3*i]-1);

		cell_iter = coc[cellsOnVertex[3*i+1]-1].find(*temp_p);
		if (cell_iter == coc[cellsOnVertex[3*i+1]-1].end()) {
			coc[cellsOnVertex[3*i+1]-1].insert(*temp_p);	
		}

		cell_iter = coc[cellsOnVertex[3*i+2]-1].find(*temp_p);
		if (cell_iter == coc[cellsOnVertex[3*i+2]-1].end()) {
			coc[cellsOnVertex[3*i+2]-1].insert(*temp_p);	
		}

		/* Second cell on vertex */
		temp_p = new Point(xCell[cellsOnVertex[3*i+1]-1], yCell[cellsOnVertex[3*i+1]-1], 0);
		temp_p->setNum(cellsOnVertex[3*i+1]-1);

		cell_iter = coc[cellsOnVertex[3*i]-1].find(*temp_p);
		if (cell_iter == coc[cellsOnVertex[3*i]-1].end()) {
			coc[cellsOnVertex[3*i]-1].insert(*temp_p);	
		}

		cell_iter = coc[cellsOnVertex[3*i+2]-1].find(*temp_p);
		if (cell_iter == coc[cellsOnVertex[3*i+2]-1].end()) {
			coc[cellsOnVertex[3*i+2]-1].insert(*temp_p);	
		}

		/* Third cell on vertex */
		temp_p = new Point(xCell[cellsOnVertex[3*i+2]-1], yCell[cellsOnVertex[3*i+2]-1], 0);
		temp_p->setNum(cellsOnVertex[3*i+2]-1);

		cell_iter = coc[cellsOnVertex[3*i]-1].find(*temp_p);
		if (cell_iter == coc[cellsOnVertex[3*i]-1].end()) {
			coc[cellsOnVertex[3*i]-1].insert(*temp_p);	
		}

		cell_iter = coc[cellsOnVertex[3*i+1]-1].find(*temp_p);
		if (cell_iter == coc[cellsOnVertex[3*i+1]-1].end()) {
			coc[cellsOnVertex[3*i+1]-1].insert(*temp_p);	
		}
	}

	cells_on_cell.resize(nCells);
	for (i=0; i<nCells; i++) {
		for (cell_iter = coc[i].begin(); cell_iter != coc[i].end(); cell_iter++) {
				cells_on_cell[i].push_back(*cell_iter);
		}
	}

	cv_on_cell.resize(nCells);
	for (i=0; i<nCells; i++) {
		for (int j=0; j<cells_on_cell[i].size(); j++) {
			cv_on_cell[i].push_back(cells_on_cell[i][j]);
		}
		for (int j=0; j<vertices_on_cell[i].size(); j++) {
			cv_on_cell[i].push_back(vertices_on_cell[i][j]);
		}
		orderCCW_normalize(cv_on_cell[i], *pset[i], x_period, y_period);
	}

	for (i=0; i<nCells; i++)
	{
		for (int j=0; j<cv_on_cell[i].size(); j++)
			cout << cv_on_cell[i][j] << endl;
		cout << cv_on_cell[i][0] << endl;
		cout << endl;
	}
#endif

        free(xCell);
        free(yCell);
        free(zCell);
        free(xVertex);
        free(yVertex);
        free(zVertex);
        free(meshDensity);
        free(cellsOnVertex);

	return 0;
}


Point segment_intersect(Point& p0, Point &p1, Point &q0, Point&q1)
{
	Point retval;

	Point u = (p1 - p0);
	Point v = (q1 - q0);
	Point w = (p0 - q0);

	double s;
	
	s = (v.getY()*w.getX() - v.getX()*w.getY())/(v.getX()*u.getY() - v.getY()*u.getX());

	retval = p0 + u*s;

	return retval;
}


void read_netcdf(int *nCells, int *nVertices, int *vertexDegree, 
		double **xCell, double **yCell, double **zCell,
		double **xVertex, double **yVertex, double **zVertex,
	 	double **meshDensity, int **cellsOnVertex,
		double *x_period, double *y_period
	  )
{
	int ncerr;
	int ncid;
	int dimIDnCells, dimIDnVertices, dimIDvertexDegree;
	int varIDxCell, varIDyCell, varIDzCell;
	int varIDxVertex, varIDyVertex, varIDzVertex;
	int varIDcellsOnVertex, varIDmeshDensity;
	size_t temp;


	ncerr = nc_open("grid.nc", NC_SHARE, &ncid);

	ncerr = nc_inq_dimid(ncid, "nCells", &dimIDnCells);
	ncerr = nc_inq_dimid(ncid, "nVertices", &dimIDnVertices);
	ncerr = nc_inq_dimid(ncid, "vertexDegree", &dimIDvertexDegree);

	ncerr = nc_inq_dimlen(ncid, dimIDnCells, &temp);
	*nCells = (int)temp;

	ncerr = nc_inq_dimlen(ncid, dimIDnVertices, &temp);
	*nVertices = (int)temp;

	ncerr = nc_inq_dimlen(ncid, dimIDvertexDegree, &temp);
	*vertexDegree = (int)temp;

	ncerr = nc_inq_varid(ncid, "xCell", &varIDxCell);
	ncerr = nc_inq_varid(ncid, "yCell", &varIDyCell);
	ncerr = nc_inq_varid(ncid, "zCell", &varIDzCell);
	ncerr = nc_inq_varid(ncid, "meshDensity", &varIDmeshDensity);
	ncerr = nc_inq_varid(ncid, "xVertex", &varIDxVertex);
	ncerr = nc_inq_varid(ncid, "yVertex", &varIDyVertex);
	ncerr = nc_inq_varid(ncid, "zVertex", &varIDzVertex);
	ncerr = nc_inq_varid(ncid, "cellsOnVertex", &varIDcellsOnVertex);

	ncerr = nc_get_att_double(ncid, NC_GLOBAL, "x_period", x_period);
	ncerr = nc_get_att_double(ncid, NC_GLOBAL, "y_period", y_period);

	*xCell = (double *)malloc(sizeof(double) * (*nCells));
	*yCell = (double *)malloc(sizeof(double) * (*nCells));
	*zCell = (double *)malloc(sizeof(double) * (*nCells));
	*meshDensity = (double *)malloc(sizeof(double) * (*nCells));
	*xVertex = (double *)malloc(sizeof(double) * (*nVertices));
	*yVertex = (double *)malloc(sizeof(double) * (*nVertices));
	*zVertex = (double *)malloc(sizeof(double) * (*nVertices));
	*cellsOnVertex = (int *)malloc(sizeof(int) * (*nVertices) * (*vertexDegree));

	ncerr = nc_get_var_double(ncid, varIDxCell, *xCell);
	ncerr = nc_get_var_double(ncid, varIDyCell, *yCell);
	ncerr = nc_get_var_double(ncid, varIDzCell, *zCell);
	ncerr = nc_get_var_double(ncid, varIDmeshDensity, *meshDensity);
	ncerr = nc_get_var_double(ncid, varIDxVertex, *xVertex);
	ncerr = nc_get_var_double(ncid, varIDyVertex, *yVertex);
	ncerr = nc_get_var_double(ncid, varIDzVertex, *zVertex);
	ncerr = nc_get_var_int(ncid, varIDcellsOnVertex, *cellsOnVertex);

	ncerr = nc_close(ncid);
}
