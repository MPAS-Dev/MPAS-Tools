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

void write_netcdf(int nCells, int nEdges, int nVertices, int maxEdges, int vertexDegree,
		int * indexToCellID, int * indexToEdgeID, int * indexToVertexID,
		double * xCell, double * yCell, double * zCell, double * latCell, double * lonCell,
		double * xEdge, double * yEdge, double * zEdge, double * latEdge, double * lonEdge,
		double * xVertex, double * yVertex, double * zVertex, double * latVertex, double * lonVertex,
		int * nEdgesOnCell, int * nEdgesOnEdge,
		int ** cellsOnCell, int ** edgesOnCell, int ** verticesOnCell,
		int ** cellsOnEdge, int ** verticesOnEdge, int ** edgesOnEdge,
		int ** edgesOnVertex, int ** cellsOnVertex, double ** kiteAreasOnVertex,
		double * fEdge, double * fVertex, double * dvEdge, double * dcEdge, double * areaCell, double * areaTriangle, double * angleEdge,
		double ** weightsOnEdge);

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
	Triangle t;
	Point p;
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
	vector<double> areaCell_v, areaTriangle_v;
	vector<double> dcEdge_v, dvEdge_v;
	vector<double> angleEdge_v;
	vector<int> nEdgesOnCell_v;
	vector<Point> cells_v;
	vector<Point> vertices_v;
	vector<Point> edges_v;
	vector<Point> edge_segments;
	set<Point>::iterator cell_iter;        /* TESTING CODE */
//	vector< vector<Point> > cv_on_cell;    /* TESTING CODE */
//	Triangle * tri;
//	double xcell, ycell;
	double x, y;
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
	 * areaCell
	 */
	areaCell_v.resize(nCells);
	for (i=0; i<nCells; i++) {
		areaCell_v[i] = poly_area(verticesOnCell_v[i]);
	}


	/*
	 * areaTriangle
	 */
	areaTriangle_v.resize(nVertices);
	for (i=0; i<nVertices; i++) {
		orderCCW_normalize(cellsOnVertex_v[i], vertices_v[i], x_period, y_period);
		t.setVertex(0,cellsOnVertex_v[i][0]);
		t.setVertex(1,cellsOnVertex_v[i][1]);
		t.setVertex(2,cellsOnVertex_v[i][2]);
		areaTriangle_v[i] = t.area();
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
	dcEdge_v.resize(nEdges);
	dvEdge_v.resize(nEdges);
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
//				cout << " Found matching vertex 1 as " << edge_segments[2].getNum() << endl;
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
//				cout << " Found matching vertex 2 as " << edge_segments[3].getNum() << endl;

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
				edgesOnCell_v[cellsOnCell_v[i][j].getNum()].push_back(edges_v[ii]);

				edgesOnVertex_v[edge_segments[2].getNum()].push_back(edges_v[ii]);
				edgesOnVertex_v[edge_segments[3].getNum()].push_back(edges_v[ii]);

				p = verticesOnEdge_v[ii][1] - verticesOnEdge_v[ii][0];
				x = p.getX();
				y = p.getY();
				dvEdge_v[ii] = sqrt(x*x + y*y);

				p = cellsOnEdge_v[ii][1] - cellsOnEdge_v[ii][0];
				x = p.getX();
				y = p.getY();
				dcEdge_v[ii] = sqrt(x*x + y*y);

				ii++;
			}
		}
	}

	
	/*
	 * nEdgesOnCell
	 */
	nEdgesOnCell_v.resize(nCells);
	for (i=0; i<nCells; i++) {
		orderCCW_normalize(edgesOnCell_v[i], cells_v[i], x_period, y_period);
		nEdgesOnCell_v[i] = edgesOnCell_v[i].size();
	}


	/*
	 * angleEdge
	 */
	angleEdge_v.resize(nEdges);
	for (i=0; i<nEdges; i++) {
		p = verticesOnEdge_v[i][0];
		p.setY(p.getY() + dvEdge_v[i]);
		angleEdge_v[i] = angle(verticesOnEdge_v[i][0], p, verticesOnEdge_v[i][1]);
	}

	
	/*
	 * Before writing output NetCDF file, ensure that all connectivity arrays are ordered properly
	 */


	/*
	 * Write NetCDF grid file
	 */

	/* Decide on a sane value for maxEdges... */


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


void write_netcdf(int nCells, int nEdges, int nVertices, int maxEdges, int vertexDegree,
                  int * indexToCellID, int * indexToEdgeID, int * indexToVertexID,
                  double * xCell, double * yCell, double * zCell, double * latCell, double * lonCell,
                  double * xEdge, double * yEdge, double * zEdge, double * latEdge, double * lonEdge,
                  double * xVertex, double * yVertex, double * zVertex, double * latVertex, double * lonVertex,
                  int * nEdgesOnCell, int * nEdgesOnEdge,
                  int ** cellsOnCell, int ** edgesOnCell, int ** verticesOnCell,
                  int ** cellsOnEdge, int ** verticesOnEdge, int ** edgesOnEdge,
                  int ** edgesOnVertex, int ** cellsOnVertex, double ** kiteAreasOnVertex,
                  double * fEdge, double * fVertex, double * dvEdge, double * dcEdge, double * areaCell, double * areaTriangle, double * angleEdge,
                  double ** weightsOnEdge
                 )
{
   int i, j, k;
   int ncerr;
   int ncid;
   int dimIDnCells, dimIDnEdges, dimIDnVertices, dimIDmaxEdges, dimIDmaxEdges2, dimIDvertexDegree, dimIDTWO;
   int varIDindexToCellID;
   int varIDxCell, varIDyCell, varIDzCell;
   int varIDlatCell, varIDlonCell;
   int varIDindexToEdgeID;
   int varIDxEdge, varIDyEdge, varIDzEdge;
   int varIDlatEdge, varIDlonEdge;
   int varIDindexToVertexID;
   int varIDxVertex, varIDyVertex, varIDzVertex;
   int varIDlatVertex, varIDlonVertex;
   int varIDnEdgesOnCell, varIDnEdgesOnEdge;
   int varIDcellsOnCell, varIDedgesOnCell, varIDverticesOnCell;
   int varIDcellsOnEdge, varIDverticesOnEdge, varIDedgesOnEdge, varIDweightsOnEdge;
   int varIDedgesOnVertex, varIDcellsOnVertex, varIDkiteAreasOnVertex;
   int varIDfEdge, varIDfVertex, varIDdvEdge, varIDdcEdge, varIDareaCell, varIDareaTriangle, varIDangleEdge;

   int dimids1[1];
   int dimids2[2];
   int dimids3[3];
   size_t start1[1], count1[1];
   size_t start2[2], count2[2];
   size_t start3[3], count3[3];

   int cellsOnCell1d[nCells*maxEdges];
   int edgesOnCell1d[nCells*maxEdges];
   int verticesOnCell1d[nCells*maxEdges];
   int cellsOnEdge1d[nEdges*2];
   int verticesOnEdge1d[nEdges*2];
   int edgesOnEdge1d[nEdges*2*maxEdges];
   double weightsOnEdge1d[nEdges*2*maxEdges];
   int edgesOnVertex1d[nVertices*vertexDegree];
   int cellsOnVertex1d[nVertices*vertexDegree];
   double kiteAreasOnVertex1d[nVertices*vertexDegree];

   double sphere_radius = 0.0;

   k = 0;
   for(i=0; i<nCells; i++) {
      for(j=0; j<maxEdges; j++) {
         cellsOnCell1d[k] = cellsOnCell[i][j];
         edgesOnCell1d[k] = edgesOnCell[i][j];
         verticesOnCell1d[k] = verticesOnCell[i][j];
         k++;
      }
   }

   k = 0;
   for(i=0; i<nEdges; i++) {
      for(j=0; j<2; j++) {
         cellsOnEdge1d[k] = cellsOnEdge[i][j];
         verticesOnEdge1d[k] = verticesOnEdge[i][j];
         k++;
      }
   }

   k = 0;
   for(i=0; i<nEdges; i++) {
      for(j=0; j<2*maxEdges; j++) {
         edgesOnEdge1d[k] = edgesOnEdge[i][j];
         weightsOnEdge1d[k] = weightsOnEdge[i][j];
         k++;
      }
   }

   k = 0;
   for(i=0; i<nVertices; i++) {
      for(j=0; j<vertexDegree; j++) {
         edgesOnVertex1d[k] = edgesOnVertex[i][j];
         cellsOnVertex1d[k] = cellsOnVertex[i][j];
         kiteAreasOnVertex1d[k] = kiteAreasOnVertex[i][j];
         k++;
      }
   }

   ncerr = nc_create("grid.nc", NC_SHARE, &ncid);

   ncerr = nc_def_dim(ncid, "nCells", (size_t)nCells, &dimIDnCells);
   ncerr = nc_def_dim(ncid, "nEdges", (size_t)nEdges, &dimIDnEdges);
   ncerr = nc_def_dim(ncid, "nVertices", (size_t)nVertices, &dimIDnVertices);
   ncerr = nc_def_dim(ncid, "maxEdges", (size_t)maxEdges, &dimIDmaxEdges);
   ncerr = nc_def_dim(ncid, "maxEdges2", (size_t)(2*maxEdges), &dimIDmaxEdges2);
   ncerr = nc_def_dim(ncid, "vertexDegree", (size_t)vertexDegree, &dimIDvertexDegree);
   ncerr = nc_def_dim(ncid, "TWO", (size_t)2, &dimIDTWO);

   dimids1[0] = dimIDnCells;
   ncerr = nc_def_var(ncid, "indexToCellID", NC_INT, 1, dimids1, &varIDindexToCellID);
   ncerr = nc_def_var(ncid, "xCell", NC_DOUBLE, 1, dimids1, &varIDxCell);
   ncerr = nc_def_var(ncid, "yCell", NC_DOUBLE, 1, dimids1, &varIDyCell);
   ncerr = nc_def_var(ncid, "zCell", NC_DOUBLE, 1, dimids1, &varIDzCell);
   ncerr = nc_def_var(ncid, "latCell", NC_DOUBLE, 1, dimids1, &varIDlatCell);
   ncerr = nc_def_var(ncid, "lonCell", NC_DOUBLE, 1, dimids1, &varIDlonCell);
   ncerr = nc_def_var(ncid, "nEdgesOnCell", NC_INT, 1, dimids1, &varIDnEdgesOnCell);
   ncerr = nc_def_var(ncid, "areaCell", NC_DOUBLE, 1, dimids1, &varIDareaCell);
   dimids1[0] = dimIDnEdges;
   ncerr = nc_def_var(ncid, "indexToEdgeID", NC_INT, 1, dimids1, &varIDindexToEdgeID);
   ncerr = nc_def_var(ncid, "xEdge", NC_DOUBLE, 1, dimids1, &varIDxEdge);
   ncerr = nc_def_var(ncid, "yEdge", NC_DOUBLE, 1, dimids1, &varIDyEdge);
   ncerr = nc_def_var(ncid, "zEdge", NC_DOUBLE, 1, dimids1, &varIDzEdge);
   ncerr = nc_def_var(ncid, "latEdge", NC_DOUBLE, 1, dimids1, &varIDlatEdge);
   ncerr = nc_def_var(ncid, "lonEdge", NC_DOUBLE, 1, dimids1, &varIDlonEdge);
   ncerr = nc_def_var(ncid, "nEdgesOnEdge", NC_INT, 1, dimids1, &varIDnEdgesOnEdge);
   ncerr = nc_def_var(ncid, "fEdge", NC_DOUBLE, 1, dimids1, &varIDfEdge);
   ncerr = nc_def_var(ncid, "dvEdge", NC_DOUBLE, 1, dimids1, &varIDdvEdge);
   ncerr = nc_def_var(ncid, "dcEdge", NC_DOUBLE, 1, dimids1, &varIDdcEdge);
   ncerr = nc_def_var(ncid, "angleEdge", NC_DOUBLE, 1, dimids1, &varIDangleEdge);
   dimids1[0] = dimIDnVertices;
   ncerr = nc_def_var(ncid, "indexToVertexID", NC_INT, 1, dimids1, &varIDindexToVertexID);
   ncerr = nc_def_var(ncid, "xVertex", NC_DOUBLE, 1, dimids1, &varIDxVertex);
   ncerr = nc_def_var(ncid, "yVertex", NC_DOUBLE, 1, dimids1, &varIDyVertex);
   ncerr = nc_def_var(ncid, "zVertex", NC_DOUBLE, 1, dimids1, &varIDzVertex);
   ncerr = nc_def_var(ncid, "latVertex", NC_DOUBLE, 1, dimids1, &varIDlatVertex);
   ncerr = nc_def_var(ncid, "lonVertex", NC_DOUBLE, 1, dimids1, &varIDlonVertex);
   ncerr = nc_def_var(ncid, "fVertex", NC_DOUBLE, 1, dimids1, &varIDfVertex);
   ncerr = nc_def_var(ncid, "areaTriangle", NC_DOUBLE, 1, dimids1, &varIDareaTriangle);
   dimids2[0] = dimIDnCells;
   dimids2[1] = dimIDmaxEdges;
   ncerr = nc_def_var(ncid, "cellsOnCell", NC_INT, 2, dimids2, &varIDcellsOnCell);
   ncerr = nc_def_var(ncid, "edgesOnCell", NC_INT, 2, dimids2, &varIDedgesOnCell);
   ncerr = nc_def_var(ncid, "verticesOnCell", NC_INT, 2, dimids2, &varIDverticesOnCell);
   dimids2[0] = dimIDnEdges;
   dimids2[1] = dimIDTWO;
   ncerr = nc_def_var(ncid, "cellsOnEdge", NC_INT, 2, dimids2, &varIDcellsOnEdge);
   ncerr = nc_def_var(ncid, "verticesOnEdge", NC_INT, 2, dimids2, &varIDverticesOnEdge);
   dimids2[0] = dimIDnEdges;
   dimids2[1] = dimIDmaxEdges2;
   ncerr = nc_def_var(ncid, "edgesOnEdge", NC_INT, 2, dimids2, &varIDedgesOnEdge);
   ncerr = nc_def_var(ncid, "weightsOnEdge", NC_DOUBLE, 2, dimids2, &varIDweightsOnEdge);
   dimids2[0] = dimIDnVertices;
   dimids2[1] = dimIDvertexDegree;
   ncerr = nc_def_var(ncid, "edgesOnVertex", NC_INT, 2, dimids2, &varIDedgesOnVertex);
   ncerr = nc_def_var(ncid, "cellsOnVertex", NC_INT, 2, dimids2, &varIDcellsOnVertex);
   ncerr = nc_def_var(ncid, "kiteAreasOnVertex", NC_DOUBLE, 2, dimids2, &varIDkiteAreasOnVertex);

   ncerr = nc_put_att_text(ncid, NC_GLOBAL, "on_a_sphere", 16, "NO              ");
   ncerr = nc_put_att_double(ncid, NC_GLOBAL, "sphere_radius", NC_DOUBLE, 1, &sphere_radius);

   ncerr = nc_enddef(ncid);

   start1[0] = 0;
   start2[0] = 0;
   start2[1] = 0;
   count1[0] = nCells;
   ncerr = nc_put_vara_int(ncid, varIDindexToCellID, start1, count1, indexToCellID);
   ncerr = nc_put_vara_double(ncid, varIDxCell, start1, count1, xCell);
   ncerr = nc_put_vara_double(ncid, varIDyCell, start1, count1, yCell);
   ncerr = nc_put_vara_double(ncid, varIDzCell, start1, count1, zCell);
   ncerr = nc_put_vara_double(ncid, varIDlatCell, start1, count1, latCell);
   ncerr = nc_put_vara_double(ncid, varIDlonCell, start1, count1, lonCell);
   ncerr = nc_put_vara_int(ncid, varIDnEdgesOnCell, start1, count1, nEdgesOnCell);
   ncerr = nc_put_vara_double(ncid, varIDareaCell, start1, count1, areaCell);
   count1[0] = nEdges;
   ncerr = nc_put_vara_int(ncid, varIDindexToEdgeID, start1, count1, indexToEdgeID);
   ncerr = nc_put_vara_double(ncid, varIDxEdge, start1, count1, xEdge);
   ncerr = nc_put_vara_double(ncid, varIDyEdge, start1, count1, yEdge);
   ncerr = nc_put_vara_double(ncid, varIDzEdge, start1, count1, zEdge);
   ncerr = nc_put_vara_double(ncid, varIDlatEdge, start1, count1, latEdge);
   ncerr = nc_put_vara_double(ncid, varIDlonEdge, start1, count1, lonEdge);
   ncerr = nc_put_vara_int(ncid, varIDnEdgesOnEdge, start1, count1, nEdgesOnEdge);
   ncerr = nc_put_vara_double(ncid, varIDfEdge, start1, count1, fEdge);
   ncerr = nc_put_vara_double(ncid, varIDdvEdge, start1, count1, dvEdge);
   ncerr = nc_put_vara_double(ncid, varIDdcEdge, start1, count1, dcEdge);
   ncerr = nc_put_vara_double(ncid, varIDangleEdge, start1, count1, angleEdge);
   count1[0] = nVertices;
   ncerr = nc_put_vara_int(ncid, varIDindexToVertexID, start1, count1, indexToVertexID);
   ncerr = nc_put_vara_double(ncid, varIDxVertex, start1, count1, xVertex);
   ncerr = nc_put_vara_double(ncid, varIDyVertex, start1, count1, yVertex);
   ncerr = nc_put_vara_double(ncid, varIDzVertex, start1, count1, zVertex);
   ncerr = nc_put_vara_double(ncid, varIDlatVertex, start1, count1, latVertex);
   ncerr = nc_put_vara_double(ncid, varIDlonVertex, start1, count1, lonVertex);
   ncerr = nc_put_vara_double(ncid, varIDfVertex, start1, count1, fVertex);
   ncerr = nc_put_vara_double(ncid, varIDareaTriangle, start1, count1, areaTriangle);
   count2[0] = nCells;
   count2[1] = maxEdges;
   ncerr = nc_put_vara_int(ncid, varIDcellsOnCell, start2, count2, cellsOnCell1d);
   ncerr = nc_put_vara_int(ncid, varIDedgesOnCell, start2, count2, edgesOnCell1d);
   ncerr = nc_put_vara_int(ncid, varIDverticesOnCell, start2, count2, verticesOnCell1d);
   count2[0] = nEdges;
   count2[1] = 2;
   ncerr = nc_put_vara_int(ncid, varIDcellsOnEdge, start2, count2, cellsOnEdge1d);
   ncerr = nc_put_vara_int(ncid, varIDverticesOnEdge, start2, count2, verticesOnEdge1d);
   count2[0] = nEdges;
   count2[1] = 2*maxEdges;
   ncerr = nc_put_vara_int(ncid, varIDedgesOnEdge, start2, count2, edgesOnEdge1d);
   ncerr = nc_put_vara_double(ncid, varIDweightsOnEdge, start2, count2, weightsOnEdge1d);
   count2[0] = nVertices;
   count2[1] = vertexDegree;
   ncerr = nc_put_vara_int(ncid, varIDedgesOnVertex, start2, count2, edgesOnVertex1d);
   ncerr = nc_put_vara_int(ncid, varIDcellsOnVertex, start2, count2, cellsOnVertex1d);
   ncerr = nc_put_vara_double(ncid, varIDkiteAreasOnVertex, start2, count2, kiteAreasOnVertex1d);

   ncerr = nc_close(ncid);
}
