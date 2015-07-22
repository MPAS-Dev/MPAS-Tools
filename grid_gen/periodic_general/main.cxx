#include <iostream>
#include <fstream>
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

#define ALLOC_INT2D(ARR,I,J) (ARR) = new int*[(I)]; for(int i=0; i<(I); i++) (ARR)[i] = new int[(J)];
#define DEALLOC_INT2D(ARR,I,J) for(int i=0; i<(I); i++) delete [] (ARR)[i]; delete [] (ARR);

#define ALLOC_REAL2D(ARR,I,J) (ARR) = new double*[(I)]; for(int i=0; i<(I); i++) (ARR)[i] = new double[(J)];
#define DEALLOC_REAL2D(ARR,I,J) for(int i=0; i<(I); i++) delete [] (ARR)[i]; delete [] (ARR);

int obtuse_triangle(Triangle &t);

void write_netcdf(int nCells, int nVertices, int vertexDegree, 
		double *xCell, double *yCell, double *zCell, 
		double *xVertex, double *yVertex, double *zVertex, 
		double *meshDensity, int *cellsOnVertex,
		double x_period, double y_period);

void readParamsFile();

Point segment_intersect(Point& p0, Point &p1, Point &q0, Point&q1);


// run-time parameters
	double EPS = 1.0e-7;
	double X_PERIOD = 1.0;
	double Y_PERIOD = 1.0;
	double X_BUFFER_FRAC = 0.05;
	double Y_BUFFER_FRAC = 0.05;
	double X_BUFFER_W = 1.0;  // Buffer width in length units - these get calculated automatically from X/Y_BUFFER_FRAC
	double Y_BUFFER_W = 1.0;
	int NUMPOINTS = 200;
	int MAXITR = 100;
	int USE_MC = 1; // 1=true, 0=read from file
	int USE_DATA_DENSITY = 0; // 1=true, 0=analytic density function in DensityFunction.cxx



int main(int argc, char ** argv)
{

// read user-specified settings
	readParamsFile();

	int i, ii, jj, n, iter, idx, npts, np;
	DensityFunction f(X_PERIOD, Y_PERIOD, USE_DATA_DENSITY);
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






	if (USE_MC == 1) {
		cout << "MC" <<endl;
		pset.makeMCPoints(NUMPOINTS, X_PERIOD, Y_PERIOD, USE_DATA_DENSITY);
	} else {
		cout << "file" <<endl;
		pset.initFromTextFile(X_PERIOD, Y_PERIOD, USE_DATA_DENSITY, "centroids.txt");
	}


	/*
	 * Set flags in point set for immovable "boundary" points
	 */
	npts = pset.size();
	for (i=0; i<npts; i++) {
		if (pset[i]->getX() < (double)( X_BUFFER_W ) || pset[i]->getX() > (double)( X_PERIOD - X_BUFFER_W )) 
			pset[i]->setBoundaryPoint(1);
		if (pset[i]->getY() < (double)( Y_BUFFER_W ) || pset[i]->getY() > (double)( Y_PERIOD - Y_BUFFER_W )) 
			pset[i]->setBoundaryPoint(1);
	}


	/*
	 * Lloyd iteration
	 */
	for (iter=0; iter<MAXITR; iter++) {
		cout << "Iteration " << iter << endl;
		vcs = pset.getVoronoiDiagram();
		npts = pset.size();
		for (i=0; i<npts; i++) {
			if (!pset[i]->isBoundaryPoint()) {
				total_mass = 0.0;
				p.setXY(0.0, 0.0);
				for (int j=0; j<vcs[i].size(); j++) {
					t = Triangle(*pset[i], vcs[i][j], vcs[i][(j+1)%vcs[i].size()]);
					p2 = t.centroid(f, &mass);
					p = p + p2 * mass;
					total_mass += mass;
				}
 
				p = p * (1.0 / total_mass);
				pset[i]->setXY(p.getX(), p.getY());
 
				/* If point has drifted into boundary region, push it back... */
				pset[i]->setX(pset[i]->getX() < (double)( X_BUFFER_W ) ? (double)( X_BUFFER_W ) : pset[i]->getX());
				pset[i]->setX(pset[i]->getX() > (double)( X_PERIOD - X_BUFFER_W ) ? (double)( X_PERIOD - X_BUFFER_W ) : pset[i]->getX());
				pset[i]->setY(pset[i]->getY() < (double)( Y_BUFFER_W ) ? (double)( Y_BUFFER_W ) : pset[i]->getY());
				pset[i]->setY(pset[i]->getY() > (double)( Y_PERIOD - Y_BUFFER_W ) ? (double)( Y_PERIOD - Y_BUFFER_W ) : pset[i]->getY());
			}
		}
		delete [] vcs;
	}

	
	restart = fopen("restart.txt","w");
	for(i=0; i<pset.size(); i++) {
		fprintf(restart, "%lf %lf\n", pset[i]->getX(), pset[i]->getY());
	}
	fclose(restart);


	/*
	 * To get a triangulation of the points, we'll need to make copies of the boundary points
	 */
	npts = pset.size();
	for (i=0; i<npts; i++) {
		temp_p = new Point(pset[i]->getX(), pset[i]->getY(), 0);
		temp_p->setNum(pset[i]->getNum());
		if (pset[i]->isBoundaryPoint())
			temp_p->setBoundaryPoint(1);
		out_pset.addPoint(*temp_p);

		/* If this is a boundary point, add it again in a periodic way */
		if (temp_p->isBoundaryPoint()) {

			if (temp_p->getX() < (double)( X_BUFFER_W )) {

				/* RIGHT SIDE */
				temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
				temp_pp->setNum(-1 * (temp_p->getNum() + 1));           /* Bdy points have negative indices */
				temp_pp->setX(temp_pp->getX() + (double)( X_PERIOD ));
				out_pset.addPoint(*temp_pp);

				if (temp_p->getY() < (double)( Y_BUFFER_W )) {

					/* UPPER-RIGHT CORNER */
					temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
					temp_pp->setNum(-1 * (temp_p->getNum() + 1));           /* Bdy points have negative indices */
					temp_pp->setX(temp_pp->getX() + (double)( X_PERIOD ));
					temp_pp->setY(temp_pp->getY() + (double)( Y_PERIOD ));
					out_pset.addPoint(*temp_pp);
				}
				else if (temp_p->getY() > (double)( Y_PERIOD - Y_BUFFER_W )) {

					/* LOWER-RIGHT CORNER */
					temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
					temp_pp->setNum(-1 * (temp_p->getNum() + 1));           /* Bdy points have negative indices */
					temp_pp->setX(temp_pp->getX() + (double)( X_PERIOD ));
					temp_pp->setY(temp_pp->getY() - (double)( Y_PERIOD ));
					out_pset.addPoint(*temp_pp);
				}
			}
			else if (temp_p->getX() > (double)( X_PERIOD - X_BUFFER_W )) {

				/* LEFT SIDE */
				temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
				temp_pp->setNum(-1 * (temp_p->getNum() + 1));           /* Bdy points have negative indices */
				temp_pp->setX(temp_pp->getX() - (double)( X_PERIOD ));
				out_pset.addPoint(*temp_pp);

				if (temp_p->getY() < (double)( Y_BUFFER_W )) {

					/* UPPER-LEFT CORNER */
					temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
					temp_pp->setNum(-1 * (temp_p->getNum() + 1));           /* Bdy points have negative indices */
					temp_pp->setX(temp_pp->getX() - (double)( X_PERIOD ));
					temp_pp->setY(temp_pp->getY() + (double)( Y_PERIOD ));
					out_pset.addPoint(*temp_pp);
				}
				else if (temp_p->getY() > (double)( Y_PERIOD - Y_BUFFER_W )) {

					/* LOWER-LEFT CORNER */
					temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
					temp_pp->setNum(-1 * (temp_p->getNum() + 1));           /* Bdy points have negative indices */
					temp_pp->setX(temp_pp->getX() - (double)( X_PERIOD ));
					temp_pp->setY(temp_pp->getY() - (double)( Y_PERIOD ));
					out_pset.addPoint(*temp_pp);
				}
			}

			if (temp_p->getY() < (double)( Y_BUFFER_W )) {

				/* TOP SIDE */
				temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
				temp_pp->setNum(-1 * (temp_p->getNum() + 1));           /* Bdy points have negative indices */
				temp_pp->setY(temp_pp->getY() + (double)( Y_PERIOD ));
				out_pset.addPoint(*temp_pp);
			}
			else if (temp_p->getY() > (double)( Y_PERIOD - Y_BUFFER_W )) {

				/* BOTTOM SIDE */
				temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
				temp_pp->setNum(-1 * (temp_p->getNum() + 1));           /* Bdy points have negative indices */
				temp_pp->setY(temp_pp->getY() - (double)( Y_PERIOD ));
				out_pset.addPoint(*temp_pp);
			}
			
		}
		
	}


	/*
	 * Having obtained a triangulation of "real" generating points as well as "ghost" points,
	 *    we need to scan through the triangles and keep a unique set that triangulates a truly
	 *    doubly-periodic grid
	 */
	triangulation = out_pset.getTriangulation();
        for (it = triangulation->begin(); it != triangulation->end(); it++) {

		/* 
		 * Ghost/halo points have a negative index; if all of the vertices of a triangle
		 *    are negative, the triangle is redundant
		 */
		ii = 0;
		for (int j=0; j<3; j++)
			if ( it->getVertex(j).getNum() >= 0 )
				ii++;
		
		/* 
		 * If at least one corner of the triangle is non-negative, we consider keeping it,
		 *    but only if it isn't redundant with another triangle already added to the set
		 */
		if ( ii > 0 ) {
			tri = new Triangle();

			for (int j=0; j<3; j++) {
				temp_p = new Point(it->getVertex(j).getX(), it->getVertex(j).getY(), 0);
				temp_p->setNum(it->getVertex(j).getNum());
	
				/* Set point number back to positive value */
				if (temp_p->getNum() < 0)
					temp_p->setNum(-1 * (temp_p->getNum() + 1));
				tri->setVertex(j, *temp_p);
			}

			dti = delaunay_tri.find(*tri);
			if (dti == delaunay_tri.end()) 
				delaunay_tri.insert(*tri);
			else
				delete tri;
		}
	}


	/*
	 * Scan through triangles and ensure that corner locations are in the range (0,X_PERIOD],(0,Y_PERIOD]
	 */
        for (dti = delaunay_tri.begin(); dti != delaunay_tri.end(); dti++) {
		t = *dti;
		t.normalizeVertices((double)( EPS ), (double)( X_PERIOD + EPS ), (double)( EPS ), (double)( Y_PERIOD + EPS ));
		norm_dt.push_back(t);
	}


	delete triangulation;


	/*
	 * Generate {x,y,z}{Cell,Vertex}, meshDensity, and cellsOnVertex fields into simple arrays
	 */
	nCells = pset.size();
	nVertices = norm_dt.size();
	vertexDegree = 3;
	cout << "nCells = " << nCells << endl;
	cout << "nVertices = " << nVertices << endl;

	xCell = (double *)malloc(sizeof(double) * (size_t)nCells);
	yCell = (double *)malloc(sizeof(double) * (size_t)nCells);
	zCell = (double *)malloc(sizeof(double) * (size_t)nCells);

	xVertex = (double *)malloc(sizeof(double) * (size_t)nVertices);
	yVertex = (double *)malloc(sizeof(double) * (size_t)nVertices);
	zVertex = (double *)malloc(sizeof(double) * (size_t)nVertices);

	meshDensity = (double *)malloc(sizeof(double) * (size_t)nCells);

	cellsOnVertex = (int *)malloc(sizeof(int) * (size_t)nVertices * (size_t)vertexDegree);
	
	npts = pset.size();
	for (i=0; i<npts; i++) {
		xCell[i] = pset[i]->getX();
		yCell[i] = pset[i]->getY();
		zCell[i] = 0.0;
		meshDensity[i] = f.evaluate(*pset[i]);
	}

	i = 0;
	ii = 0;
        for (norm_dti = norm_dt.begin(); norm_dti != norm_dt.end(); norm_dti++) {
		p =  norm_dti->circumcenter();
		xVertex[i] = p.getX();
		yVertex[i] = p.getY();
		zVertex[i] = 0.0;
		for (int j=0; j<3; j++)
			cellsOnVertex[ii++] = norm_dti->getVertex(j).getNum() + 1;   /* indices are 1-based in MPAS */
//			cellsOnVertex[ii++] = norm_dti->getVertex(j).getNum();       /* Do not use the 0-based indices when making meshes for MPAS */
		i++;
	}


	/*
	 * Write fields to NetCDF file
	 */
	write_netcdf(nCells, nVertices, vertexDegree, xCell, yCell, zCell, xVertex, yVertex, zVertex, meshDensity, cellsOnVertex, (double)( X_PERIOD ), (double)( Y_PERIOD ));


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


int obtuse_triangle(Triangle &t) 
{
	int i;
	Point p[3];
	double PI = 2.0 * acos(0.0);

	p[0] = t.getVertex(0);
	p[1] = t.getVertex(1);
	p[2] = t.getVertex(2);

	for(i=0; i<3; i++) {
		if (fabs(angle(p[i], p[(i+1)%3], p[(i+2)%3])) > PI/2.0) {
cout << p[i] << " " << p[(i+1)%3] << " " << p[(i+2)%3] << endl;
			return i+1;
		}
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
	ncerr = nc_put_att_text(ncid, NC_GLOBAL, "is_periodic", 16, "YES             ");
	ncerr = nc_put_att_double(ncid, NC_GLOBAL, "sphere_radius", NC_DOUBLE, 1, &sphere_radius);
	ncerr = nc_put_att_double(ncid, NC_GLOBAL, "x_offset", NC_DOUBLE, 1, &x_period);
	ncerr = nc_put_att_double(ncid, NC_GLOBAL, "y_offset", NC_DOUBLE, 1, &y_period);

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



/* ***** Setup Routines ***** */
void readParamsFile(){
	//Read in parameters from Params.
	//If Params doesn't exist, write out Params with a default set of parameters
	string junk;
	ifstream params("Params.txt");
	int temp_restart_mode;
	int temp_fileio_mode;

	if(!params){
		cout << "Error opening Params.txt file." << endl;
		cout << "Writing a default Params.txt file." << endl;
		cout << "Exiting, please set up Params.txt, and rerun." << endl;
		ofstream pout("Params.txt");
		pout << "Convergence tolerance to use:" << endl;
		pout << EPS << endl;
		pout << "Maximum number of iterations to perform:" << endl;
		pout << MAXITR << endl;
		pout << "How to get initial pointset. 0=from file; 1=Monte Carlo points from density function" << endl;
		pout << USE_MC << endl;
		pout << "If using Monte Carlo points, how many do you want?" << endl;
		pout << NUMPOINTS << endl;
		pout << "Domain width (x)" << endl;
		pout << X_PERIOD << endl;
		pout << "Domain height (y)" << endl;
		pout << Y_PERIOD << endl;
		pout << "Fraction of domain to set as a buffer in which initial point locations remain fixed, x-direction" << endl;
		pout << X_BUFFER_FRAC << endl;
		pout << "Fraction of domain to set as a buffer in which initial point locations remain fixed, y-direction" << endl;
		pout << Y_BUFFER_FRAC << endl;
		pout << "Use data density in file named density.nc with variables x, y, density. 1=true, 0=analytic density function in DensityFunction.cxx" << endl;
		pout << USE_DATA_DENSITY << endl;

		pout.close();

		exit(1);
	}
	

	getline(params,junk);
	params >> EPS;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> MAXITR;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> USE_MC;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> NUMPOINTS;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> X_PERIOD;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> Y_PERIOD;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> X_BUFFER_FRAC;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> Y_BUFFER_FRAC;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> USE_DATA_DENSITY;
	params.ignore(10000,'\n');

	params.close();

	cout << "=== Specified settings are: ===" << endl;
	cout << "Convergence tolerance to use:" << endl;
	cout << EPS << endl;
	cout << "Maximum number of iterations to perform:" << endl;
	cout << MAXITR << endl;
	cout << "How to get initial pointset. 0=from file; 1=Monte Carlo points from density function" << endl;
	cout << USE_MC << endl;
	cout << "If using Monte Carlo points, how many do you want?" << endl;
	cout << NUMPOINTS << endl;
	cout << "Domain width (x)" << endl;
	cout << X_PERIOD << endl;
	cout << "Domain height (y)" << endl;
	cout << Y_PERIOD << endl;
	cout << "Fraction of domain to set as a buffer in which initial point locations remain fixed, x-direction" << endl;
	cout << X_BUFFER_FRAC << endl;
	cout << "Fraction of domain to set as a buffer in which initial point locations remain fixed, y-direction" << endl;
	cout << Y_BUFFER_FRAC << endl;
	cout << "Use data density in file named density.nc with variables x, y, density. 1=true, 0=analytic density function in DensityFunction.cxx" << endl;
	cout << USE_DATA_DENSITY << endl;

	X_BUFFER_W = X_PERIOD * X_BUFFER_FRAC;
	Y_BUFFER_W = Y_PERIOD * Y_BUFFER_FRAC;
}
