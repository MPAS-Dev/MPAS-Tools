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
	int i, j, k, ii, jj;
//	DensityFunction f;
//	PointSet out_pset;
//	Point * cells;
	Point * temp_p;
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
	vector< set<Point> > cellsOnCell_temp;
	vector< vector<Point> > cellsOnVertex_v;
	vector< vector<Point> > edgesOnVertex_v;
	vector< vector<Point> > verticesOnCell_v;
	vector< vector<Point> > cellsOnCell_v;
	vector< vector<Point> > cellsOnEdge_v;
	vector< vector<Point> > verticesOnEdge_v;
	vector< vector<Point> > edgesOnCell_v;
	vector<Point> cells_v;
	vector<Point> vertices_v;
	vector<Point> edges_v;
	vector<Point> edge_segments;
	set<Point>::iterator cell_iter;        /* TESTING CODE */
//	vector< vector<Point> > cv_on_cell;    /* TESTING CODE */
//	Triangle * tri;
//	double xcell, ycell;
//	double x, y;
//	double total_mass, mass; 
//	FILE * restart;
	int nCells, nVertices, nEdges, vertexDegree;
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

	/*
	 * vector of cells
	 */
	cells_v.resize(nCells);
	for (i=0; i<nCells; i++) {
		temp_p = new Point(xCell[i], yCell[i], 0);
		temp_p->setNum(i);
		cells_v[i] = *temp_p;
	}


	/*
	 * vector of vertices
	 */
	vertices_v.resize(nVertices);
	for (i=0; i<nVertices; i++) {
		temp_p = new Point(xVertex[i], yVertex[i], 0);
		temp_p->setNum(i);
		vertices_v[i] = *temp_p;
	}


	/*
	 * cellsOnVertex
	 */
	cellsOnVertex_v.resize(nVertices);
	for (i=0; i<nVertices; i++) {
		cellsOnVertex_v[i].resize(3);
		for (j=0; j<3; j++) {
			cellsOnVertex_v[i][j] = cells_v[cellsOnVertex[3*i+j]];
		}
	}	


	/*
	 * verticesOnCell
	 */
	verticesOnCell_v.resize(nCells);
	for (i=0; i<nVertices; i++) {
		for (int j=0; j<3; j++) {
			verticesOnCell_v[cellsOnVertex_v[i][j].getNum()].push_back(vertices_v[i]);
		}
	}


	/*
	 * cellsOnCell
	 */
	cellsOnCell_temp.resize(nCells);
	for (i=0; i<nVertices; i++) {
		/* First cell on vertex */
		cell_iter      = cellsOnCell_temp[cellsOnVertex_v[i][0].getNum()].find(cellsOnVertex_v[i][1]);
		if (cell_iter == cellsOnCell_temp[cellsOnVertex_v[i][0].getNum()].end()) {
		                 cellsOnCell_temp[cellsOnVertex_v[i][0].getNum()].insert(cellsOnVertex_v[i][1]);
		}

		cell_iter      = cellsOnCell_temp[cellsOnVertex_v[i][0].getNum()].find(cellsOnVertex_v[i][2]);
		if (cell_iter == cellsOnCell_temp[cellsOnVertex_v[i][0].getNum()].end()) {
		                 cellsOnCell_temp[cellsOnVertex_v[i][0].getNum()].insert(cellsOnVertex_v[i][2]);
		}

		/* Second cell on vertex */
		cell_iter      = cellsOnCell_temp[cellsOnVertex_v[i][1].getNum()].find(cellsOnVertex_v[i][0]);
		if (cell_iter == cellsOnCell_temp[cellsOnVertex_v[i][1].getNum()].end()) {
		                 cellsOnCell_temp[cellsOnVertex_v[i][1].getNum()].insert(cellsOnVertex_v[i][0]);
		}

		cell_iter      = cellsOnCell_temp[cellsOnVertex_v[i][1].getNum()].find(cellsOnVertex_v[i][2]);
		if (cell_iter == cellsOnCell_temp[cellsOnVertex_v[i][1].getNum()].end()) {
		                 cellsOnCell_temp[cellsOnVertex_v[i][1].getNum()].insert(cellsOnVertex_v[i][2]);
		}

		/* Third cell on vertex */
		cell_iter      = cellsOnCell_temp[cellsOnVertex_v[i][2].getNum()].find(cellsOnVertex_v[i][0]);
		if (cell_iter == cellsOnCell_temp[cellsOnVertex_v[i][2].getNum()].end()) {
		                 cellsOnCell_temp[cellsOnVertex_v[i][2].getNum()].insert(cellsOnVertex_v[i][0]);
		}
		else {
		}

		cell_iter      = cellsOnCell_temp[cellsOnVertex_v[i][2].getNum()].find(cellsOnVertex_v[i][1]);
		if (cell_iter == cellsOnCell_temp[cellsOnVertex_v[i][2].getNum()].end()) {
		                 cellsOnCell_temp[cellsOnVertex_v[i][2].getNum()].insert(cellsOnVertex_v[i][1]);
		}

	}

	cellsOnCell_v.resize(nCells);
	for (i=0; i<nCells; i++) {
		for (cell_iter = cellsOnCell_temp[i].begin(); cell_iter != cellsOnCell_temp[i].end(); cell_iter++) {
				cellsOnCell_v[i].push_back(*cell_iter);
		}
	}

	/* Place cellsOnCell, verticesOnCell in CCW order */
	for (i=0; i<nCells; i++) {
		orderCCW_normalize(cellsOnCell_v[i], cells_v[i], x_period, y_period);
		orderCCW_normalize(verticesOnCell_v[i], cells_v[i], x_period, y_period);
	}


	/*
	 * cellsOnEdge, verticesOnEdge
	 */
	nEdges = nVertices + nCells;
	edges_v.resize(nEdges);
	cellsOnEdge_v.resize(nEdges);
	verticesOnEdge_v.resize(nEdges);
	edgesOnCell_v.resize(nCells);
	edgesOnVertex_v.resize(nVertices);
	edge_segments.resize(4);   /* Used to hold the two vertices and two cells that determine the edge location */
	ii = 0;
	for (i=0; i<nCells; i++) {
		edge_segments[0] = cells_v[i];
		for (j=0; j<cellsOnCell_v[i].size(); j++) {
			if (cells_v[i].getNum() < cellsOnCell_v[i][j].getNum()) {

				/* Need to scan through all vertices and find the two that are shared with cellsOnCell_v[i][j] */
				edge_segments[1] = cellsOnCell_v[i][j];
				for (k=0; k<verticesOnCell_v[i].size(); k++) {
					for (int kk=0; kk<verticesOnCell_v[cellsOnCell_v[i][j].getNum()].size(); kk++) {
						if (verticesOnCell_v[i][k].getNum() == verticesOnCell_v[cellsOnCell_v[i][j].getNum()][kk].getNum()) {
							edge_segments[2] = verticesOnCell_v[i][k];
							goto foo;
						}
					}
				}
				cerr << "No matching vertex 1..." << endl;
foo:
				cout << " Found matching vertex 1 as " << edge_segments[2].getNum() << endl;
				for (k=k+1; k<verticesOnCell_v[i].size(); k++) {
					for (int kk=0; kk<verticesOnCell_v[cellsOnCell_v[i][j].getNum()].size(); kk++) {
						if (verticesOnCell_v[i][k].getNum() == verticesOnCell_v[cellsOnCell_v[i][j].getNum()][kk].getNum()) {
							edge_segments[3] = verticesOnCell_v[i][k];
							goto foo2;
						}
					}
				}
				cerr << "No matching vertex 2..." << endl;
foo2:
				cout << " Found matching vertex 2 as " << edge_segments[3].getNum() << endl;

				periodic_normalize(edge_segments, x_period, y_period);
				cellsOnEdge_v[ii].resize(2);
				cellsOnEdge_v[ii][0] = edge_segments[0];
				cellsOnEdge_v[ii][1] = edge_segments[1];
				verticesOnEdge_v[ii].resize(2);
				verticesOnEdge_v[ii][0] = edge_segments[2];
				verticesOnEdge_v[ii][1] = edge_segments[3];
				edges_v[ii] = segment_intersect(edge_segments[0], edge_segments[1], edge_segments[2], edge_segments[3]);
				edges_v[ii].setNum(ii);
				edgesOnCell_v[i].push_back(edges_v[ii]);
				edgesOnVertex_v[edge_segments[2].getNum()].push_back(edges_v[ii]);
				edgesOnVertex_v[edge_segments[3].getNum()].push_back(edges_v[ii]);
				ii++;

	/* Should also add the edge to cellOnCell[i][j]'s list of edges */
			}
		}
	}


	/*
	 * dcEdge, dvEdge
	 */
	for (i=0; i<nEdges; i++) {

	}

cout << "nEdges = " << nEdges << " " << ii << endl;
for (i=0; i<nEdges; i++)
{
	cout << edges_v[i] << endl;
}
cout << endl;

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
