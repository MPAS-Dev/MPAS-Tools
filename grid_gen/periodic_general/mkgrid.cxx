#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
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

void write_netcdf(int nCells, int nVertices, int vertexDegree, 
		double *xCell, double *yCell, double *zCell, 
		double *xVertex, double *yVertex, double *zVertex, 
		double *meshDensity, int *cellsOnVertex,
		double x_period, double y_period);

Point segment_intersect(Point& p0, Point &p1, Point &q0, Point&q1);

int main(int argc, char ** argv)
{
	int i, ii, jj, n, iter, idx, npts, np;
	DensityFunction f;
	PointSet pset;
	PointSet out_pset;
	vector<Point> * vcs;
	Point * cells;
	Point * temp_p;
	Point * temp_pp;
	Point p3;
	Triangle t;
	Point p, p2;
	vector<Point> * clist;
	vector<Triangle> * triangulation;
	vector<Triangle>::iterator it;
	set<Triangle> delaunay_tri;
	set<Triangle>::iterator dti;
	list<Triangle> norm_dt;
	list<Triangle>::iterator norm_dti;
	vector< vector<Point> > vertices_on_cell;
	vector< vector<Point> > cells_on_cell;
	vector< set<Point> > coc;
	set<Point>::iterator cell_iter;
	vector< vector<Point> > cv_on_cell;
	Triangle * tri;
	vector<Point> * vlist;
	vector<Point> * elist;
	double xcell, ycell;
	double x, y;
	double total_mass, mass; 
	FILE * restart;
	int nCells, nVertices, vertexDegree;
	double *xCell, *yCell, *zCell, *xVertex, *yVertex, *zVertex, *meshDensity;
	int *cellsOnVertex;
	double x_period, y_period;


	/*
	 * Read basic grid info from NetCDF file
	 */
	write_netcdf(nCells, nVertices, vertexDegree, xCell, yCell, zCell, xVertex, yVertex, zVertex, meshDensity, cellsOnVertex, x_period, y_period);

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


void write_netcdf(int nCells, int nVertices, int vertexDegree, 
		double * xCell, double * yCell, double * zCell,
		double * xVertex, double * yVertex, double * zVertex,
	 	double * meshDensity, int * cellsOnVertex,
		double x_period, double y_period
	  )
{
	int i, j, k;
	int ncerr;
	int ncid;
	int dimIDnCells, dimIDnVertices, dimIDvertexDegree;
	int varIDxCell, varIDyCell, varIDzCell;
	int varIDxVertex, varIDyVertex, varIDzVertex;
	int varIDcellsOnVertex, varIDmeshDensity;
	
	int dimids1[1];
	int dimids2[2];
	int dimids3[3];
	size_t start1[1], count1[1];
	size_t start2[2], count2[2];
	size_t start3[3], count3[3];

	double sphere_radius = 0.0;


	ncerr = nc_create("grid.nc", NC_SHARE, &ncid);

	ncerr = nc_def_dim(ncid, "nCells", (size_t)nCells, &dimIDnCells);
	ncerr = nc_def_dim(ncid, "nVertices", (size_t)nVertices, &dimIDnVertices);
	ncerr = nc_def_dim(ncid, "vertexDegree", (size_t)vertexDegree, &dimIDvertexDegree);

	dimids1[0] = dimIDnCells;
	ncerr = nc_def_var(ncid, "xCell", NC_DOUBLE, 1, dimids1, &varIDxCell);
	ncerr = nc_def_var(ncid, "yCell", NC_DOUBLE, 1, dimids1, &varIDyCell);
	ncerr = nc_def_var(ncid, "zCell", NC_DOUBLE, 1, dimids1, &varIDzCell);
	ncerr = nc_def_var(ncid, "meshDensity", NC_DOUBLE, 1, dimids1, &varIDmeshDensity);
	dimids1[0] = dimIDnVertices;
	ncerr = nc_def_var(ncid, "xVertex", NC_DOUBLE, 1, dimids1, &varIDxVertex);
	ncerr = nc_def_var(ncid, "yVertex", NC_DOUBLE, 1, dimids1, &varIDyVertex);
	ncerr = nc_def_var(ncid, "zVertex", NC_DOUBLE, 1, dimids1, &varIDzVertex);
	dimids2[0] = dimIDnVertices;
	dimids2[1] = dimIDvertexDegree;
	ncerr = nc_def_var(ncid, "cellsOnVertex", NC_INT, 2, dimids2, &varIDcellsOnVertex);

	ncerr = nc_put_att_text(ncid, NC_GLOBAL, "on_a_sphere", 16, "NO              ");
	ncerr = nc_put_att_double(ncid, NC_GLOBAL, "sphere_radius", NC_DOUBLE, 1, &sphere_radius);
	ncerr = nc_put_att_double(ncid, NC_GLOBAL, "x_period", NC_DOUBLE, 1, &x_period);
	ncerr = nc_put_att_double(ncid, NC_GLOBAL, "y_period", NC_DOUBLE, 1, &y_period);

	ncerr = nc_enddef(ncid);

	start1[0] = 0;
	start2[0] = 0;
	start2[1] = 0;
	count1[0] = nCells;
	ncerr = nc_put_vara_double(ncid, varIDxCell, start1, count1, xCell);
	ncerr = nc_put_vara_double(ncid, varIDyCell, start1, count1, yCell);
	ncerr = nc_put_vara_double(ncid, varIDzCell, start1, count1, zCell);
	ncerr = nc_put_vara_double(ncid, varIDmeshDensity, start1, count1, meshDensity);
	count1[0] = nVertices;
	ncerr = nc_put_vara_double(ncid, varIDxVertex, start1, count1, xVertex);
	ncerr = nc_put_vara_double(ncid, varIDyVertex, start1, count1, yVertex);
	ncerr = nc_put_vara_double(ncid, varIDzVertex, start1, count1, zVertex);
	count2[0] = nVertices;
	count2[1] = vertexDegree;
	ncerr = nc_put_vara_int(ncid, varIDcellsOnVertex, start2, count2, cellsOnVertex);

	ncerr = nc_close(ncid);
}
