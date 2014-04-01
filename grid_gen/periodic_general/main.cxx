#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include "PointSet.h"
#include "Triangle.h"
#include "DensityFunction.h"
#include "netcdf.h"
using namespace std;

/* Sets the period of the grid in x and y */
#define X_PERIOD 40.0
#define Y_PERIOD 40.0*0.866025403784439

/* Sets the width of the zone of cells that are immovable along the x and y boundaries */
#define X_BUFFER_W 5.0
#define Y_BUFFER_W 5.0

#define ALLOC_INT2D(ARR,I,J) (ARR) = new int*[(I)]; for(int i=0; i<(I); i++) (ARR)[i] = new int[(J)];
#define DEALLOC_INT2D(ARR,I,J) for(int i=0; i<(I); i++) delete [] (ARR)[i]; delete [] (ARR);

#define ALLOC_REAL2D(ARR,I,J) (ARR) = new double*[(I)]; for(int i=0; i<(I); i++) (ARR)[i] = new double[(J)];
#define DEALLOC_REAL2D(ARR,I,J) for(int i=0; i<(I); i++) delete [] (ARR)[i]; delete [] (ARR);

void compute_grid_meta(int nPoints, Point * cells, vector<Point> * alist, vector<Point> * clist, vector<Point> * elist);
int obtuse_triangle(Triangle &t);
Point segment_intersect(Point& p0, Point &p1, Point &q0, Point&q1);
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
	  );


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
	Triangle * tri;
	vector<Point> * vlist;
	vector<Point> * elist;
	double xcell, ycell;
	double total_mass, mass; 
	ifstream cellsOnCell("cellsOnCell.txt");
	ifstream verticesOnCell("verticesOnCell.txt");
	ifstream edgesOnCell("edgesOnCell.txt");
	FILE * restart;

	const int MAXITR = 100;

	pset.initFromTextFile("centroids.txt");


	/*
	 * Set flags in point set for immovable "boundary" points
	 */
	npts = pset.size();
	for (i=0; i<npts; i++) {
		if (pset[i]->getX() < (float)( X_BUFFER_W ) || pset[i]->getX() > (float)( X_PERIOD - X_BUFFER_W )) 
			pset[i]->setBoundaryPoint(1);
		if (pset[i]->getY() < (float)( Y_BUFFER_W ) || pset[i]->getY() > (float)( Y_PERIOD - Y_BUFFER_W )) 
			pset[i]->setBoundaryPoint(1);
	}


	/*
	 * Lloyd iteration
	 */
	for (iter=0; iter<MAXITR; iter++) {
		cout << "Iteration " << iter << endl;
		fprintf(stderr, "Iteration %i\n", iter);
		vcs = pset.getVoronoiDiagram();
		fprintf(stderr, "	got Voronoi diagram\n");
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
				pset[i]->setX(pset[i]->getX() < (float)( X_BUFFER_W ) ? (float)( X_BUFFER_W ) : pset[i]->getX());
				pset[i]->setX(pset[i]->getX() > (float)( X_PERIOD - X_BUFFER_W ) ? (float)( X_PERIOD - X_BUFFER_W ) : pset[i]->getX());
				pset[i]->setY(pset[i]->getY() < (float)( Y_BUFFER_W ) ? (float)( Y_BUFFER_W ) : pset[i]->getY());
				pset[i]->setY(pset[i]->getY() > (float)( Y_PERIOD - Y_BUFFER_W ) ? (float)( Y_PERIOD - Y_BUFFER_W ) : pset[i]->getY());
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

			if (temp_p->getX() < (float)( X_BUFFER_W )) {

				/* RIGHT SIDE */
				temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
				temp_pp->setNum(-1 * (temp_p->getNum() + 1));
				temp_pp->setX(temp_pp->getX() + (float)( X_PERIOD ));
				out_pset.addPoint(*temp_pp);

				if (temp_p->getY() < (float)( Y_BUFFER_W )) {

					/* UPPER-RIGHT CORNER */
					temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
					temp_pp->setNum(-1 * (temp_p->getNum() + 1));
					temp_pp->setX(temp_pp->getX() + (float)( X_PERIOD ));
					temp_pp->setY(temp_pp->getY() + (float)( Y_PERIOD ));
					out_pset.addPoint(*temp_pp);
				}
				else if (temp_p->getY() > (float)( Y_PERIOD - Y_BUFFER_W )) {

					/* LOWER-RIGHT CORNER */
					temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
					temp_pp->setNum(-1 * (temp_p->getNum() + 1));
					temp_pp->setX(temp_pp->getX() + (float)( X_PERIOD ));
					temp_pp->setY(temp_pp->getY() - (float)( Y_PERIOD ));
					out_pset.addPoint(*temp_pp);
				}
			}
			else if (temp_p->getX() > (float)( X_PERIOD - X_BUFFER_W )) {

				/* LEFT SIDE */
				temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
				temp_pp->setNum(-1 * (temp_p->getNum() + 1));
				temp_pp->setX(temp_pp->getX() - (float)( X_PERIOD ));
				out_pset.addPoint(*temp_pp);

				if (temp_p->getY() < (float)( Y_BUFFER_W )) {

					/* UPPER-LEFT CORNER */
					temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
					temp_pp->setNum(-1 * (temp_p->getNum() + 1));
					temp_pp->setX(temp_pp->getX() - (float)( X_PERIOD ));
					temp_pp->setY(temp_pp->getY() + (float)( Y_PERIOD ));
					out_pset.addPoint(*temp_pp);
				}
				else if (temp_p->getY() > (float)( Y_PERIOD - Y_BUFFER_W )) {

					/* LOWER-LEFT CORNER */
					temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
					temp_pp->setNum(-1 * (temp_p->getNum() + 1));
					temp_pp->setX(temp_pp->getX() - (float)( X_PERIOD ));
					temp_pp->setY(temp_pp->getY() - (float)( Y_PERIOD ));
					out_pset.addPoint(*temp_pp);
				}
			}

			if (temp_p->getY() < (float)( Y_BUFFER_W )) {

				/* TOP SIDE */
				temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
				temp_pp->setNum(-1 * (temp_p->getNum() + 1));
				temp_pp->setY(temp_pp->getY() + (float)( Y_PERIOD ));
				out_pset.addPoint(*temp_pp);
			}
			else if (temp_p->getY() > (float)( Y_PERIOD - Y_BUFFER_W )) {

				/* BOTTOM SIDE */
				temp_pp = new Point(temp_p->getX(), temp_p->getY(), 0);
				temp_pp->setNum(-1 * (temp_p->getNum() + 1));
				temp_pp->setY(temp_pp->getY() - (float)( Y_PERIOD ));
				out_pset.addPoint(*temp_pp);
			}
			
		}
		
	}

	out_pset.printToTextFile("debug.dat");

	triangulation = out_pset.getTriangulation();
        for (it = triangulation->begin(); it != triangulation->end(); it++) {
		ii = 0;
		for (int j=0; j<3; j++) {
			if (it->getVertex(j).getNum() > 0)
				ii++;
		}

		if (ii > 0) {
			tri = new Triangle();

			for (int j=0; j<3; j++) {
				tri->setVertex(j, it->getVertex(j));
				if (tri->getVertex(j).getNum() < 0)
					tri->getVertex(j).setNum(-1 * (tri->getVertex(j).getNum() - 1));
			}
			dti = delaunay_tri.find(*tri);
			if (dti == delaunay_tri.end()) 
				delaunay_tri.insert(*tri);
		}
	}

        for (dti = delaunay_tri.begin(); dti != delaunay_tri.end(); dti++) {
		cout << "TRI " << dti->getVertex(0).getNum() << " " << dti->getVertex(1).getNum() << " " << dti->getVertex(2).getNum() << endl;
	}

	delete triangulation;

#if 0
	pset.printToTextFile("cvt.dat");


	cells = new Point[pset.size()];
	clist = pset.getDelaunayAdjacency();
	vlist = pset.getVoronoiDiagram();
	elist = new vector<Point>[pset.size()];

	assert(cellsOnCell.is_open());
	assert(verticesOnCell.is_open());
	assert(edgesOnCell.is_open());

	for(i=0; i<pset.size(); i++) {
		cells[i] = *pset[i];
		cells[i].setNum(i);
		if (pset[i]->isBoundaryPoint() == 2) {
			clist[i].clear();
			for(int j=0; j<6; j++) {
				cellsOnCell >> idx;
				idx--;
				assert(idx >= 0 && idx < pset.size()); 
				clist[i].push_back(*pset[idx]);
			}

			vlist[i].clear();
			for(int j=0; j<6; j++) {
				verticesOnCell >> xcell >> ycell;
				p3.setXY(xcell, ycell);
				vlist[i].push_back(p3);
			}

		}
	}

	for(i=0; i<pset.size(); i++) {
		if (pset[i]->isBoundaryPoint() == 2) {
			elist[i].clear();
			for(int j=0; j<6; j++) {
				edgesOnCell >> xcell >> ycell;
				p3.setXY(xcell, ycell);
				elist[i].push_back(p3);
			}
		}
		else {
			// NB: we don't need to worry about periodic cells, since this block only executes for interior cells
			for(int j=0; j<clist[i].size(); j++) {
				p3 = (cells[i] + clist[i][j]) * 0.5;
				elist[i].push_back(p3);
			}
		}
	}

	cellsOnCell.close();
	verticesOnCell.close();
	edgesOnCell.close();

	compute_grid_meta(pset.size(), cells, clist, vlist, elist);
	
	delete [] cells;
	delete [] clist;
	delete [] vlist;
	delete [] elist;
#endif

	return 0;
}


void compute_grid_meta(int nPoints, Point * points, vector<Point> * clist, vector<Point> * vlist, vector<Point> * elist)
{
	int found;
	int shared_edge, shared_vtx1, shared_vtx2;
	double d;
	double area, area_sum, kite_area, s;
	int i, j, k, ii, jj, kk, last;
	const int maxEdges = 9;
	int nCells, nEdges, nVertices;
	double *latCell, *lonCell, *latEdge, *lonEdge, *latVertex, *lonVertex;
	double *xCell, *yCell, *zCell, *xEdge, *yEdge, *zEdge, *xVertex, *yVertex, *zVertex;
	double *dcEdge, *dvEdge;
	int *indexToCellID, *indexToEdgeID, *indexToVertexID;
	int *nEdgesOnCell, *nEdgesOnEdge;
	int **cellsOnCell, **verticesOnCell, **edgesOnCell;
	int **cellsOnVertex, **edgesOnVertex;
	int **cellsOnEdge, **verticesOnEdge, **edgesOnEdge;
	double *areaCell, *areaTriangle;
	double **kiteAreasOnVertex, **weightsOnEdge;
	double *angleEdge;
	double *fEdge, *fVertex;
	Point *cells, *edges, *vertices;
	vector<Point> evec;
	Point o, p1, p2, p3, cell1, cell2, vtx1, vtx2;
	Triangle t(p1, p2, p3);  // We don't care about initial triangle vertices
	ofstream graph_info("graph.info");

#if OLDCODE
	nCells = nPoints;
cout << "nCells=" << nCells << endl;
  
	nVertices = 0;
	for(i=0; i<nPoints; i++)
		nVertices += clist[i].size();
	assert(nVertices % 3 == 0);
	nVertices = nVertices / 3;
cout << "nVertices=" << nVertices << endl;

	nEdges = nVertices + nCells;	// Euler characteristic 0
cout << "nEdges=" << nEdges << endl;

	indexToCellID = new int[nCells];
	indexToEdgeID = new int[nEdges];
	indexToVertexID = new int[nVertices];

	for(i=0; i<nCells; i++)
		indexToCellID[i] = i+1;

	for(i=0; i<nEdges; i++)
		indexToEdgeID[i] = i+1;

	for(i=0; i<nVertices; i++)
		indexToVertexID[i] = i+1;

	latCell = new double[nCells]; 
	lonCell = new double[nCells]; 

	for(i=0; i<nCells; i++) {
		latCell[i] = 0.0;
		lonCell[i] = 0.0;
	}

	latEdge = new double[nEdges];
	lonEdge = new double[nEdges];

	for(i=0; i<nEdges; i++) {
		latEdge[i] = 0.0;
		lonEdge[i] = 0.0;
	}

	latVertex = new double[nVertices];
	lonVertex = new double[nVertices];

	for(i=0; i<nVertices; i++) {
		latVertex[i] = 0.0;
		lonVertex[i] = 0.0;
	}

	cells = new Point[nCells];
	edges = new Point[nEdges];
	vertices = new Point[nVertices];

	for(i=0; i<nCells; i++) {
		cells[i] = points[i];
	}


	last = -1;
	for(i=0; i<nCells; i++) {
		for(j=0; j<elist[i].size(); j++) {
			found = 0; 
			for(k=0; k<=last; k++) {
				if (elist[i][j].distance(edges[k]) < 0.001) {found = 1; break;}
			}
			if (!found) {edges[++last] = elist[i][j]; edges[last].setNum(last);}
		}
	}
cout << "last=" << last << endl;

	last = -1;
	for(i=0; i<nCells; i++) {
		for(j=0; j<vlist[i].size(); j++) {
			found = 0; 
			for(k=0; k<=last; k++) {
				if (vlist[i][j].distance(vertices[k]) < 0.001) {found = 1; break;}
			}
			if (!found) {vertices[++last] = vlist[i][j]; vertices[last].setNum(last);}
		}
	}
cout << "last=" << last << endl;


	xCell = new double[nCells];
	yCell = new double[nCells];
	zCell = new double[nCells];

	xEdge = new double[nEdges];
	yEdge = new double[nEdges];
	zEdge = new double[nEdges];

	xVertex = new double[nVertices];
	yVertex = new double[nVertices];
	zVertex = new double[nVertices];

	for(i=0; i<nCells; i++) {
		xCell[i] = cells[i].getX();
		yCell[i] = cells[i].getY();
		zCell[i] = 0.0;
	}

	for(i=0; i<nEdges; i++) {
		xEdge[i] = edges[i].getX();
		yEdge[i] = edges[i].getY();
		zEdge[i] = 0.0;
	}

	for(i=0; i<nVertices; i++) {
		xVertex[i] = vertices[i].getX();
		yVertex[i] = vertices[i].getY();
		zVertex[i] = 0.0;
	}


	nEdgesOnCell = new int[nCells];
	ALLOC_INT2D(cellsOnCell, nCells, maxEdges)
	ALLOC_INT2D(edgesOnCell, nCells, maxEdges)
	ALLOC_INT2D(verticesOnCell, nCells, maxEdges)

	for(i=0; i<nCells; i++) {
		nEdgesOnCell[i] = clist[i].size();
		for(j=0; j<clist[i].size(); j++) {
			for(k=0; k<nCells; k++) {
				if (cells[k].distance(clist[i][j]) < 0.001) {cellsOnCell[i][j] = k+1; break;}
			}
			assert(k<nCells);
		}
	}

	for(i=0; i<nCells; i++) {
		for(j=0; j<elist[i].size(); j++) {
			for(k=0; k<nEdges; k++) {
				if (edges[k].distance(elist[i][j]) < 0.001) {edgesOnCell[i][j] = k+1; break;}
			}
			assert(k<nEdges);
		}
	}

	for(i=0; i<nCells; i++) {
		for(j=0; j<vlist[i].size(); j++) {
			for(k=0; k<nVertices; k++) {
				if (vertices[k].distance(vlist[i][j]) < 0.001) {verticesOnCell[i][j] = k+1; break;}
			}
			assert(k<nVertices);
		}
	}

	ALLOC_INT2D(cellsOnEdge, nEdges, 2)
	ALLOC_INT2D(cellsOnVertex, nVertices, 3)

	for(i=0; i<nEdges; i++) {
		cellsOnEdge[i][0] = 0;
		cellsOnEdge[i][1] = 0;
	}

	for(i=0; i<nVertices; i++) {
		cellsOnVertex[i][0] = 0;
		cellsOnVertex[i][1] = 0;
		cellsOnVertex[i][2] = 0;
	}

	for(i=0; i<nCells; i++) {
		for(j=0; j<nEdgesOnCell[i]; j++) {
			if (cellsOnEdge[edgesOnCell[i][j]-1][0] == 0) {
				cellsOnEdge[edgesOnCell[i][j]-1][0] = i+1;
			}
			else if (cellsOnEdge[edgesOnCell[i][j]-1][1] == 0)
				cellsOnEdge[edgesOnCell[i][j]-1][1] = i+1;
			else {
				assert(0);
			}
		}
	}

	for(i=0; i<nCells; i++) {
		for(j=0; j<nEdgesOnCell[i]; j++) {
			if (cellsOnVertex[verticesOnCell[i][j]-1][0] == 0)
				cellsOnVertex[verticesOnCell[i][j]-1][0] = i+1;
			else if (cellsOnVertex[verticesOnCell[i][j]-1][1] == 0)
				cellsOnVertex[verticesOnCell[i][j]-1][1] = i+1;
			else if (cellsOnVertex[verticesOnCell[i][j]-1][2] == 0)
				cellsOnVertex[verticesOnCell[i][j]-1][2] = i+1;
			else
				assert(0);
		}
	}

for(i=0; i<nVertices; i++) {
	assert(cellsOnVertex[i][2] != 0);
}

	ALLOC_INT2D(verticesOnEdge, nEdges, 2)
	ALLOC_INT2D(edgesOnVertex, nVertices, 3)

	for(i=0; i<nEdges; i++) {
		verticesOnEdge[i][0] = 0;
		verticesOnEdge[i][1] = 0;
	}

	for(i=0; i<nVertices; i++) {
		edgesOnVertex[i][0] = 0;
		edgesOnVertex[i][1] = 0;
		edgesOnVertex[i][2] = 0;
	}

	for(i=0; i<nCells; i++) {
		for(j=0; j<nEdgesOnCell[i]; j++) {
			if ((i+1) == cellsOnEdge[edgesOnCell[i][j]-1][0]) ii = cellsOnEdge[edgesOnCell[i][j]-1][1] - 1;
			else ii = cellsOnEdge[edgesOnCell[i][j]-1][0] - 1;
// Search through edge lists of i and ii to find common edge and common vertices

			shared_edge = 0;
			for(k=0; k<nEdgesOnCell[i]; k++) {
				for(kk=0; kk<nEdgesOnCell[ii]; kk++) {
					if (edgesOnCell[i][k] == edgesOnCell[ii][kk]) {
						shared_edge = edgesOnCell[i][k]; break;
					}
				}
				if (shared_edge) break;
			}
assert(shared_edge);

			shared_vtx1 = 0;
			for(k=0; k<nEdgesOnCell[i]; k++) {
				for(kk=0; kk<nEdgesOnCell[ii]; kk++) {
					if (verticesOnCell[i][k] == verticesOnCell[ii][kk]) {
						shared_vtx1 = verticesOnCell[i][k]; break;
					}
				}
				if (shared_vtx1) break;
			}
assert(shared_vtx1);

			shared_vtx2 = 0;
			for(k=0; k<nEdgesOnCell[i]; k++) {
				for(kk=0; kk<nEdgesOnCell[ii]; kk++) {
					if ((verticesOnCell[i][k] == verticesOnCell[ii][kk]) && (verticesOnCell[i][k] != shared_vtx1)) {
						shared_vtx2 = verticesOnCell[i][k]; break;
					}
				}
				if (shared_vtx2) break;
			}
assert(shared_vtx2);

			verticesOnEdge[shared_edge-1][0] = shared_vtx1;
			verticesOnEdge[shared_edge-1][1] = shared_vtx2;
	
			// Check that order of verticesOnEdge is correct (maintains the right-hand rule)
			cell1 = cells[cellsOnEdge[shared_edge-1][0]-1];
			cell2 = cells[cellsOnEdge[shared_edge-1][1]-1];
			vtx1 = vertices[verticesOnEdge[shared_edge-1][0]-1];
			vtx2 = vertices[verticesOnEdge[shared_edge-1][1]-1];
			if (cell1.distance(cell2) > 20000.0) {
				if (cell1.getX() - cell2.getX() > 20000.0)
					cell1.setX(cell1.getX() - 1000.0*NCOLS);
				else if (cell2.getX() - cell1.getX() > 20000.0)
					cell2.setX(cell2.getX() - 1000.0*NCOLS);
				if (cell1.getY() - cell2.getY() > 20000.0)
					cell1.setY(cell1.getY() - 1000.0*NCOLS);
				else if (cell2.getY() - cell1.getY() > 20000.0)
					cell2.setY(cell2.getY() - 1000.0*NCOLS);
			}
			if (cell1.distance(vtx1) > 20000.0) {
				if (vtx1.getX() - cell1.getX() > 20000.0)
					vtx1.setX(vtx1.getX() - 1000.0*NCOLS);
				else if (cell1.getX() - vtx1.getX() > 20000.0)
					vtx1.setX(vtx1.getX() + 1000.0*NCOLS);
				if (vtx1.getY() - cell1.getY() > 20000.0)
					vtx1.setY(vtx1.getY() - 1000.0*NCOLS);
				else if (cell1.getY() - vtx1.getY() > 20000.0)
					vtx1.setY(vtx1.getY() + 1000.0*NCOLS);
			}
			if (cell1.distance(vtx2) > 20000.0) {
				if (vtx2.getX() - cell1.getX() > 20000.0)
					vtx2.setX(vtx2.getX() - 1000.0*NCOLS);
				else if (cell1.getX() - vtx2.getX() > 20000.0)
					vtx2.setX(vtx2.getX() + 1000.0*NCOLS);
				if (vtx2.getY() - cell1.getY() > 20000.0)
					vtx2.setY(vtx2.getY() - 1000.0*NCOLS);
				else if (cell1.getY() - vtx2.getY() > 20000.0)
					vtx2.setY(vtx2.getY() + 1000.0*NCOLS);
			}
			o.setXY(0.0, 0.0);
			p1 = cell2 - cell1;
			p2 = vtx2 - vtx1;

			if (angle(o, p1, p2) < 0.0) {
				verticesOnEdge[shared_edge-1][0] = shared_vtx2;
				verticesOnEdge[shared_edge-1][1] = shared_vtx1;
// Could even redo check here just to be sure...
			}

			if (edgesOnVertex[shared_vtx1-1][0] != shared_edge &&
				 edgesOnVertex[shared_vtx1-1][1] != shared_edge && 
				 edgesOnVertex[shared_vtx1-1][2] != shared_edge) {
				if (edgesOnVertex[shared_vtx1-1][0] == 0) edgesOnVertex[shared_vtx1-1][0] = shared_edge;
				else if (edgesOnVertex[shared_vtx1-1][1] == 0) edgesOnVertex[shared_vtx1-1][1] = shared_edge;
				else if (edgesOnVertex[shared_vtx1-1][2] == 0) edgesOnVertex[shared_vtx1-1][2] = shared_edge;
				else assert(0);
			}

			if (edgesOnVertex[shared_vtx2-1][0] != shared_edge &&
				 edgesOnVertex[shared_vtx2-1][1] != shared_edge && 
				 edgesOnVertex[shared_vtx2-1][2] != shared_edge) {
				if (edgesOnVertex[shared_vtx2-1][0] == 0) edgesOnVertex[shared_vtx2-1][0] = shared_edge;
				else if (edgesOnVertex[shared_vtx2-1][1] == 0) edgesOnVertex[shared_vtx2-1][1] = shared_edge;
				else if (edgesOnVertex[shared_vtx2-1][2] == 0) edgesOnVertex[shared_vtx2-1][2] = shared_edge;
				else assert(0);
			}

		}
	}

	ALLOC_REAL2D(kiteAreasOnVertex, nVertices, 3)

	for(i=0; i<nVertices; i++) {
		kiteAreasOnVertex[i][0] = 0.0;
		kiteAreasOnVertex[i][1] = 0.0;
		kiteAreasOnVertex[i][2] = 0.0;
	}
#endif

#if 0
	// This code is not correct, in particular, when a ccw sort of edges and vertices leaves
	//	 two edges or vertices in the list consecutively (the j%2 test breaks down)
	for(i=0; i<nCells; i++) {
		evec.clear();
		o = cells[i];
		for(j=0; j<nEdgesOnCell[i]; j++) {
			p1 = edges[edgesOnCell[i][j]-1];
			if (p1.distance(o) > 20000.0) {
				if (p1.getX() - o.getX() > 20000.0)
					p1.setX(p1.getX() - 1000.0*NCOLS);
				else if (o.getX() - p1.getX() > 20000.0)
					p1.setX(p1.getX() + 1000.0*NCOLS);
				if (p1.getY() - o.getY() > 20000.0)
					p1.setY(p1.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
				else if (o.getY() - p1.getY() > 20000.0)
					p1.setY(p1.getY() + 1000.0*NROWS*sqrt(3.0)/2.0);
			}
			evec.push_back(p1);

			p1 = vertices[verticesOnCell[i][j]-1];
			if (p1.distance(o) > 20000.0) {
				if (p1.getX() - o.getX() > 20000.0)
					p1.setX(p1.getX() - 1000.0*NCOLS);
				else if (o.getX() - p1.getX() > 20000.0)
					p1.setX(p1.getX() + 1000.0*NCOLS);
				if (p1.getY() - o.getY() > 20000.0)
					p1.setY(p1.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
				else if (o.getY() - p1.getY() > 20000.0)
					p1.setY(p1.getY() + 1000.0*NROWS*sqrt(3.0)/2.0);
			}
			evec.push_back(p1);
		}
		orderCCW(evec, o);

		t.setVertex(0, o); 
		for(j=0; j<evec.size(); j++) {
			if (j%2 == 0) area = 0.0;
			t.setVertex(1, evec[j]); 
			t.setVertex(2, evec[(j+1)%evec.size()]); 
			area += t.area();
			if (j%2 == 1) {
				kk = evec[j].getNum();
				if (cellsOnVertex[kk][0] == (i+1)) kiteAreasOnVertex[kk][0] = area;
				else if (cellsOnVertex[kk][1] == (i+1)) kiteAreasOnVertex[kk][1] = area;
				else if (cellsOnVertex[kk][2] == (i+1)) kiteAreasOnVertex[kk][2] = area;
				else assert(0);
//				if (kiteAreasOnVertex[kk][0] == 0) kiteAreasOnVertex[kk][0] = area;
//				else if (kiteAreasOnVertex[kk][1] == 0) kiteAreasOnVertex[kk][1] = area;
//				else if (kiteAreasOnVertex[kk][2] == 0) kiteAreasOnVertex[kk][2] = area;
//				else assert(0);
			} 
		}
	}
#endif


#ifdef OLDCODE
	ALLOC_INT2D(edgesOnEdge, nEdges, 2*maxEdges)
	nEdgesOnEdge = new int[nEdges];
  
	for(i=0; i<nEdges; i++) {
		for(j=0; j<2*maxEdges; j++)
			edgesOnEdge[i][j] = 0;
	}

	evec.clear();
	for(i=0; i<nEdges; i++) {
		k = 0;
		o = cells[cellsOnEdge[i][0]-1];
//		evec.push_back(edges[i]);
		p1 = edges[i];
		if (p1.distance(o) > 20000.0) {
			if (p1.getX() - o.getX() > 20000.0)
				p1.setX(p1.getX() - 1000.0*NCOLS);
			else if (o.getX() - p1.getX() > 20000.0)
				p1.setX(p1.getX() + 1000.0*NCOLS);
			if (p1.getY() - o.getY() > 20000.0)
				p1.setY(p1.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
			else if (o.getY() - p1.getY() > 20000.0)
				p1.setY(p1.getY() + 1000.0*NROWS*sqrt(3.0)/2.0);
		}
		evec.push_back(p1);

		kk = cellsOnEdge[i][0]-1;
		for(j=0; j<nEdgesOnCell[kk]; j++) {
			if (edgesOnCell[kk][j] != (i+1)) {
				p1 = edges[edgesOnCell[kk][j]-1];
				if (p1.distance(o) > 20000.0) {
					if (p1.getX() - o.getX() > 20000.0)
						p1.setX(p1.getX() - 1000.0*NCOLS);
					else if (o.getX() - p1.getX() > 20000.0)
						p1.setX(p1.getX() + 1000.0*NCOLS);
					if (p1.getY() - o.getY() > 20000.0)
						p1.setY(p1.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
					else if (o.getY() - p1.getY() > 20000.0)
						p1.setY(p1.getY() + 1000.0*NROWS*sqrt(3.0)/2.0);
				}
				evec.push_back(p1);
			}
		}
		orderCCW_print(evec, o);
		for(ii=1; ii<evec.size(); ii++)
			edgesOnEdge[i][k++] = evec[ii].getNum() + 1;
		evec.clear();

		o = cells[cellsOnEdge[i][1]-1];
//		evec.push_back(edges[i]);
		p1 = edges[i];
		if (p1.distance(o) > 20000.0) {
			if (p1.getX() - o.getX() > 20000.0)
				p1.setX(p1.getX() - 1000.0*NCOLS);
			else if (o.getX() - p1.getX() > 20000.0)
				p1.setX(p1.getX() + 1000.0*NCOLS);
			if (p1.getY() - o.getY() > 20000.0)
				p1.setY(p1.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
			else if (o.getY() - p1.getY() > 20000.0)
				p1.setY(p1.getY() + 1000.0*NROWS*sqrt(3.0)/2.0);
		}
		evec.push_back(p1);

		kk = cellsOnEdge[i][1]-1;
		for(j=0; j<nEdgesOnCell[kk]; j++) {
			if (edgesOnCell[kk][j] != (i+1)) {
				p1 = edges[edgesOnCell[kk][j]-1];
				if (p1.distance(o) > 20000.0) {
					if (p1.getX() - o.getX() > 20000.0)
						p1.setX(p1.getX() - 1000.0*NCOLS);
					else if (o.getX() - p1.getX() > 20000.0)
						p1.setX(p1.getX() + 1000.0*NCOLS);
					if (p1.getY() - o.getY() > 20000.0)
						p1.setY(p1.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
					else if (o.getY() - p1.getY() > 20000.0)
						p1.setY(p1.getY() + 1000.0*NROWS*sqrt(3.0)/2.0);
				}
				evec.push_back(p1);
			}
		}
		orderCCW_print(evec, o);
		for(ii=1; ii<evec.size(); ii++)
			edgesOnEdge[i][k++] = evec[ii].getNum() + 1;
		evec.clear();
		
		nEdgesOnEdge[i] = k;
//		cout << "Edge " << i+1 << " has " << k << " neighbors" << endl;
//		cout << "	on cells " << cellsOnEdge[i][0] << " " << cellsOnEdge[i][1] << endl;

//cout << "	Neighboring edges are: ";
//		for(j=0; j<nEdgesOnEdge[i]; j++)
//			cout << edgesOnEdge[i][j] << " "; 
//cout << endl;
	}

	dcEdge = new double[nEdges];
	dvEdge = new double[nEdges];

	for(i=0; i<nEdges; i++) {
		p1 = cells[cellsOnEdge[i][0]-1];
		p2 = cells[cellsOnEdge[i][1]-1];
		dcEdge[i] = p1.distance(p2);
//MGD remove the magic number 20000 later...
		if (dcEdge[i] > 20000.0) {
			if (p1.getX() - p2.getX() > 20000.0)
				p1.setX(p1.getX() - 1000.0*NCOLS);
			else if (p2.getX() - p1.getX() > 20000.0)
				p2.setX(p2.getX() - 1000.0*NCOLS);
			if (p1.getY() - p2.getY() > 20000.0)
				p1.setY(p1.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
			else if (p2.getY() - p1.getY() > 20000.0)
				p2.setY(p2.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
		}
		dcEdge[i] = p1.distance(p2);
//assert(fabs(dcEdge[i] - 1000.0) < 0.001);

		p1 = vertices[verticesOnEdge[i][0]-1];
		p2 = vertices[verticesOnEdge[i][1]-1];
		dvEdge[i] = p1.distance(p2);
//MGD remove the magic number 20000 later...
		if (dvEdge[i] > 20000.0) {
			if (p1.getX() - p2.getX() > 20000.0)
				p1.setX(p1.getX() - 1000.0*NCOLS);
			else if (p2.getX() - p1.getX() > 20000.0)
				p2.setX(p2.getX() - 1000.0*NCOLS);
			if (p1.getY() - p2.getY() > 20000.0)
				p1.setY(p1.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
			else if (p2.getY() - p1.getY() > 20000.0)
				p2.setY(p2.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
		}
		dvEdge[i] = p1.distance(p2);
//assert(fabs(dvEdge[i] - 577.35) < 0.001);
	}

  
	areaCell = new double[nCells];
	areaTriangle = new double[nVertices];

	for(i=0; i<nVertices; i++) {
		p1 = cells[cellsOnVertex[i][0]-1];
		p2 = cells[cellsOnVertex[i][1]-1];
		p3 = cells[cellsOnVertex[i][2]-1];
	
		if (p1.distance(p2) > 20000.0) {
			if (p1.getX() - p2.getX() > 20000.0)
				p1.setX(p1.getX() - 1000.0*NCOLS);
			else if (p2.getX() - p1.getX() > 20000.0)
				p2.setX(p2.getX() - 1000.0*NCOLS);
			if (p1.getY() - p2.getY() > 20000.0)
				p1.setY(p1.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
			else if (p2.getY() - p1.getY() > 20000.0)
				p2.setY(p2.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
		}
		if (p1.distance(p3) > 20000.0) {
			if (p3.getX() - p1.getX() > 20000.0)
				p3.setX(p3.getX() - 1000.0*NCOLS);
			else if (p1.getX() - p3.getX() > 20000.0)
				p3.setX(p3.getX() + 1000.0*NCOLS);
			if (p3.getY() - p1.getY() > 20000.0)
				p3.setY(p3.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
			else if (p1.getY() - p3.getY() > 20000.0)
				p3.setY(p3.getY() + 1000.0*NROWS*sqrt(3.0)/2.0);
		}

		t.setVertex(0, p1);
		t.setVertex(1, p2);
		t.setVertex(2, p3);

		areaTriangle[i] = t.area();

//assert(fabs(areaTriangle[i] - 433012.705) < 0.001);
	}

	for(i=0; i<nCells; i++) {
		areaCell[i] = 0.0;
		for(j=0; j<nEdgesOnCell[i]; j++) {
			p1 = cells[i];
			p2 = vertices[verticesOnCell[i][j]-1];
			p3 = vertices[verticesOnCell[i][(j+1)%nEdgesOnCell[i]]-1];

			if (p1.distance(p2) > 20000.0) {
				if (p1.getX() - p2.getX() > 20000.0)
					p1.setX(p1.getX() - 1000.0*NCOLS);
				else if (p2.getX() - p1.getX() > 20000.0)
					p2.setX(p2.getX() - 1000.0*NCOLS);
				if (p1.getY() - p2.getY() > 20000.0)
					p1.setY(p1.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
				else if (p2.getY() - p1.getY() > 20000.0)
					p2.setY(p2.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
			}
			if (p1.distance(p3) > 20000.0) {
				if (p3.getX() - p1.getX() > 20000.0)
					p3.setX(p3.getX() - 1000.0*NCOLS);
				else if (p1.getX() - p3.getX() > 20000.0)
					p3.setX(p3.getX() + 1000.0*NCOLS);
				if (p3.getY() - p1.getY() > 20000.0)
					p3.setY(p3.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
				else if (p1.getY() - p3.getY() > 20000.0)
					p3.setY(p3.getY() + 1000.0*NROWS*sqrt(3.0)/2.0);
			}

			t.setVertex(0, p1);
			t.setVertex(1, p2);
			t.setVertex(2, p3);

			areaCell[i] += t.area();
			
		}

//assert(fabs(areaCell[i] - 866025.4) < 0.001);
	}

//MGD maybe set kiteAreasOnVertex here...
	for(i=0; i<nVertices; i++) {
		p1 = cells[cellsOnVertex[i][0]-1];
		p2 = cells[cellsOnVertex[i][1]-1];
		p3 = cells[cellsOnVertex[i][2]-1];

		if (p1.distance(p2) > 20000.0) {
			if (p1.getX() - p2.getX() > 20000.0)
				p1.setX(p1.getX() - 1000.0*NCOLS);
			else if (p2.getX() - p1.getX() > 20000.0)
				p2.setX(p2.getX() - 1000.0*NCOLS);
			if (p1.getY() - p2.getY() > 20000.0)
				p1.setY(p1.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
			else if (p2.getY() - p1.getY() > 20000.0)
				p2.setY(p2.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
		}
		if (p1.distance(p3) > 20000.0) {
			if (p3.getX() - p1.getX() > 20000.0)
				p3.setX(p3.getX() - 1000.0*NCOLS);
			else if (p1.getX() - p3.getX() > 20000.0)
				p3.setX(p3.getX() + 1000.0*NCOLS);
			if (p3.getY() - p1.getY() > 20000.0)
				p3.setY(p3.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
			else if (p1.getY() - p3.getY() > 20000.0)
				p3.setY(p3.getY() + 1000.0*NROWS*sqrt(3.0)/2.0);
		}

		t.setVertex(0, p1);
		t.setVertex(1, p2);
		t.setVertex(2, p3);
		evec.clear();
		if ((ii = obtuse_triangle(t)) > 0) {
			ii--;
			printf("Vertex %i of triangle %i is obtuse\n", ii, i);

			evec.push_back(t.getVertex(ii));
			evec.push_back(t.getVertex((ii+1)%3));
			evec.push_back(t.getVertex((ii+2)%3));

			evec.push_back(edges[edgesOnVertex[i][0]-1]);
			evec.push_back(edges[edgesOnVertex[i][1]-1]);
			evec.push_back(edges[edgesOnVertex[i][2]-1]);

			if (evec[3].distance(evec[4]) > 20000.0) {
				if (evec[3].getX() - evec[4].getX() > 20000.0)
					evec[3].setX(evec[3].getX() - 1000.0*NCOLS);
				else if (evec[4].getX() - evec[3].getX() > 20000.0)
					evec[4].setX(evec[4].getX() - 1000.0*NCOLS);
				if (evec[3].getY() - evec[4].getY() > 20000.0)
					evec[3].setY(evec[3].getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
				else if (evec[4].getY() - evec[3].getY() > 20000.0)
					evec[4].setY(evec[4].getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
			}
			if (evec[3].distance(evec[5]) > 20000.0) {
				if (evec[5].getX() - evec[3].getX() > 20000.0)
					evec[5].setX(evec[5].getX() - 1000.0*NCOLS);
				else if (evec[3].getX() - evec[5].getX() > 20000.0)
					evec[5].setX(evec[5].getX() + 1000.0*NCOLS);
				if (evec[5].getY() - evec[3].getY() > 20000.0)
					evec[5].setY(evec[5].getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
				else if (evec[3].getY() - evec[5].getY() > 20000.0)
					evec[5].setY(evec[5].getY() + 1000.0*NROWS*sqrt(3.0)/2.0);
			}

			o = t.centroid();
			orderCCW(evec,o);

			p1 = segment_intersect(evec[5], vertices[i], evec[4], evec[2]);
			p2 = segment_intersect(evec[1], vertices[i], evec[4], evec[2]);
cout << "Intersection points are " << p1 << " and " << p2 << endl;

			t.setVertex(0, evec[4]);
			t.setVertex(1, evec[5]);
			t.setVertex(2, p1);
			area = fabs(t.area());
			area_sum = area;
			if (cellsOnVertex[i][0]-1 == evec[4].getNum()) kiteAreasOnVertex[i][0] = area;
			else if (cellsOnVertex[i][1]-1 == evec[4].getNum()) kiteAreasOnVertex[i][1] = area;
			else if (cellsOnVertex[i][2]-1 == evec[4].getNum()) kiteAreasOnVertex[i][2] = area;
			else assert(1 == 2);

			t.setVertex(0, evec[2]);
			t.setVertex(1, evec[1]);
			t.setVertex(2, p2);
			area = fabs(t.area());
			area_sum += area;
			if (cellsOnVertex[i][0]-1 == evec[2].getNum()) kiteAreasOnVertex[i][0] = area;
			else if (cellsOnVertex[i][1]-1 == evec[2].getNum()) kiteAreasOnVertex[i][1] = area;
			else if (cellsOnVertex[i][2]-1 == evec[2].getNum()) kiteAreasOnVertex[i][2] = area;
			else assert(1 == 2);

			area = areaTriangle[i] - area_sum;
			if (cellsOnVertex[i][0]-1 == evec[0].getNum()) kiteAreasOnVertex[i][0] = area;
			else if (cellsOnVertex[i][1]-1 == evec[0].getNum()) kiteAreasOnVertex[i][1] = area;
			else if (cellsOnVertex[i][2]-1 == evec[0].getNum()) kiteAreasOnVertex[i][2] = area;
			else assert(1 == 2);

		}
		else {
			evec.push_back(edges[edgesOnVertex[i][0]-1]);
			evec.push_back(edges[edgesOnVertex[i][1]-1]);
			evec.push_back(edges[edgesOnVertex[i][2]-1]);

			if (evec[0].distance(evec[1]) > 20000.0) {
				if (evec[0].getX() - evec[1].getX() > 20000.0)
					evec[0].setX(evec[0].getX() - 1000.0*NCOLS);
				else if (evec[1].getX() - evec[0].getX() > 20000.0)
					evec[1].setX(evec[1].getX() - 1000.0*NCOLS);
				if (evec[0].getY() - evec[1].getY() > 20000.0)
					evec[0].setY(evec[0].getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
				else if (evec[1].getY() - evec[0].getY() > 20000.0)
					evec[1].setY(evec[1].getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
			}
			if (evec[0].distance(evec[2]) > 20000.0) {
				if (evec[2].getX() - evec[0].getX() > 20000.0)
					evec[2].setX(evec[2].getX() - 1000.0*NCOLS);
				else if (evec[0].getX() - evec[2].getX() > 20000.0)
					evec[2].setX(evec[2].getX() + 1000.0*NCOLS);
				if (evec[2].getY() - evec[0].getY() > 20000.0)
					evec[2].setY(evec[2].getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
				else if (evec[0].getY() - evec[2].getY() > 20000.0)
					evec[2].setY(evec[2].getY() + 1000.0*NROWS*sqrt(3.0)/2.0);
			}

			evec.push_back(p1);
			evec.push_back(p2);
			evec.push_back(p3);
			o = t.centroid();
			orderCCW(evec,o);
			t.setVertex(0, o);
			for(j=0; j<6; j+=2) {
				t.setVertex(1,evec[j]);
				t.setVertex(2,evec[j+1]);
				area = fabs(t.area());
				t.setVertex(1,evec[j+1]);
				t.setVertex(2,evec[(j+2)%6]);
				area += fabs(t.area());
				if (cellsOnVertex[i][0]-1 == evec[j+1].getNum()) kiteAreasOnVertex[i][0] = area;
				else if (cellsOnVertex[i][1]-1 == evec[j+1].getNum()) kiteAreasOnVertex[i][1] = area;
				else if (cellsOnVertex[i][2]-1 == evec[j+1].getNum()) kiteAreasOnVertex[i][2] = area;
				else assert(1 == 2);
			}
		}
	}


	fEdge = new double[nEdges];
	fVertex = new double[nVertices];

	for(i=0; i<nEdges; i++)
		fEdge[i] = 0.0;

	for(i=0; i<nVertices; i++)
		fVertex[i] = 0.0;

	angleEdge = new double[nEdges];

	for(i=0; i<nEdges; i++) {
		p1 = vertices[verticesOnEdge[i][0]-1];
		p2 = vertices[verticesOnEdge[i][1]-1];
		if (p1.distance(p2) > 20000.0) {
			if (p1.getX() - p2.getX() > 20000.0)
				p1.setX(p1.getX() - 1000.0*NCOLS);
			else if (p2.getX() - p1.getX() > 20000.0)
				p2.setX(p2.getX() - 1000.0*NCOLS);
			if (p1.getY() - p2.getY() > 20000.0)
				p1.setY(p1.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
			else if (p2.getY() - p1.getY() > 20000.0)
				p2.setY(p2.getY() - 1000.0*NROWS*sqrt(3.0)/2.0);
		}

		o.setXY(p1.getX(), p1.getY() + 1000.0);
		angleEdge[i] = angle(p1, o, p2);
	}

	ALLOC_REAL2D(weightsOnEdge, nEdges, 2*maxEdges);

	for(i=0; i<nEdges; i++) {
		for(j=0; j<2*maxEdges; j++) 
			weightsOnEdge[i][j] = 0.0;
	}
	for(i=0; i<nEdges; i++) {
		ii = cellsOnEdge[i][0]-1;
		shared_vtx1 = verticesOnEdge[i][1];
		area_sum = 0.0;

		jj = 0;
		for(j=0; j<nEdgesOnCell[ii]-1; j++) {
			if (cellsOnVertex[shared_vtx1-1][0] == (ii+1))		 kite_area = kiteAreasOnVertex[shared_vtx1-1][0];
			else if (cellsOnVertex[shared_vtx1-1][1] == (ii+1))	kite_area = kiteAreasOnVertex[shared_vtx1-1][1];
			else if (cellsOnVertex[shared_vtx1-1][2] == (ii+1))	kite_area = kiteAreasOnVertex[shared_vtx1-1][2];
			else assert(0);
			area_sum += kite_area;

			if (cellsOnEdge[edgesOnEdge[i][jj]-1][0] == (ii+1)) s = 1.0;
			else s = -1.0;
	 
			weightsOnEdge[i][jj] = s * (0.5 - area_sum / areaCell[ii]) * dvEdge[edgesOnEdge[i][jj]-1] / dcEdge[i];

			if (verticesOnEdge[edgesOnEdge[i][jj]-1][0] != shared_vtx1) shared_vtx1 = verticesOnEdge[edgesOnEdge[i][jj]-1][0];
			else shared_vtx1 = verticesOnEdge[edgesOnEdge[i][jj]-1][1];
			jj++;
		}

		ii = cellsOnEdge[i][1]-1;
		area_sum = 0.0;

		for(j=0; j<nEdgesOnCell[ii]-1; j++) {
			if (cellsOnVertex[shared_vtx1-1][0] == (ii+1))		 kite_area = kiteAreasOnVertex[shared_vtx1-1][0];
			else if (cellsOnVertex[shared_vtx1-1][1] == (ii+1))	kite_area = kiteAreasOnVertex[shared_vtx1-1][1];
			else if (cellsOnVertex[shared_vtx1-1][2] == (ii+1))	kite_area = kiteAreasOnVertex[shared_vtx1-1][2];
			else assert(0);
			area_sum += kite_area;

			if (cellsOnEdge[edgesOnEdge[i][jj]-1][0] == (ii+1)) s = -1.0;
			else s = 1.0;
	 
			weightsOnEdge[i][jj] = s * (0.5 - area_sum / areaCell[ii]) * dvEdge[edgesOnEdge[i][jj]-1] / dcEdge[i];

			if (verticesOnEdge[edgesOnEdge[i][jj]-1][0] != shared_vtx1) shared_vtx1 = verticesOnEdge[edgesOnEdge[i][jj]-1][0];
			else shared_vtx1 = verticesOnEdge[edgesOnEdge[i][jj]-1][1];
			jj++;
		}
assert(jj == nEdgesOnEdge[i]);
	}

for(i=0; i<nEdges; i++) {
	if (sqrt(pow(xEdge[i]-50000.0,2.0) + pow(yEdge[i]-43301.25,2.0)) < 23687.0)
		cout << "NEW EDGE" << xEdge[i] << " " << yEdge[i] << endl;
}


	write_netcdf(nCells, nEdges, nVertices, maxEdges, 3,
					 indexToCellID, indexToEdgeID, indexToVertexID,
					 xCell, yCell, zCell, latCell, lonCell,
					 xEdge, yEdge, zEdge, latEdge, lonEdge,
					 xVertex, yVertex, zVertex, latVertex, lonVertex,
					 nEdgesOnCell, nEdgesOnEdge,
					 cellsOnCell, edgesOnCell, verticesOnCell,
					 cellsOnEdge, verticesOnEdge, edgesOnEdge,
					 edgesOnVertex, cellsOnVertex, kiteAreasOnVertex,
					 fEdge, fVertex, dvEdge, dcEdge, areaCell, areaTriangle, angleEdge,
					 weightsOnEdge
					);

	assert(graph_info.is_open());
	graph_info << nCells << " " << nEdges << endl;
	for(i=0; i<nCells; i++) {
		for(j=0; j<nEdgesOnCell[i]; j++)
			graph_info << cellsOnCell[i][j] << " ";
		graph_info << endl;
	}
	graph_info.close();

	delete [] indexToCellID;
	delete [] indexToEdgeID;
	delete [] indexToVertexID;
	delete [] latCell;
	delete [] lonCell;
	delete [] latEdge;
	delete [] lonEdge;
	delete [] latVertex;
	delete [] lonVertex;
	delete [] cells;
	delete [] edges;
	delete [] vertices;
	delete [] xCell;
	delete [] yCell;
	delete [] zCell;
	delete [] xEdge;
	delete [] yEdge;
	delete [] zEdge;
	delete [] xVertex;
	delete [] yVertex;
	delete [] zVertex;
	delete [] nEdgesOnCell;
	delete [] nEdgesOnEdge;
	delete [] dcEdge;
	delete [] dvEdge;
	delete [] areaCell;
	delete [] areaTriangle;
	delete [] fEdge;
	delete [] fVertex;
	delete [] angleEdge;
	DEALLOC_INT2D(cellsOnCell, nCells, maxEdges)
	DEALLOC_INT2D(edgesOnCell, nCells, maxEdges)
	DEALLOC_INT2D(verticesOnCell, nCells, maxEdges)
	DEALLOC_INT2D(cellsOnEdge, nEdges, 2)
	DEALLOC_INT2D(verticesOnEdge, nEdges, 2)
	DEALLOC_INT2D(cellsOnVertex, nVertices, 3)
	DEALLOC_INT2D(edgesOnVertex, nVertices, 3)
	DEALLOC_INT2D(edgesOnEdge, nEdges, 2*maxEdges)
	DEALLOC_REAL2D(kiteAreasOnVertex, nVertices, 3)
	DEALLOC_REAL2D(weightsOnEdge, nEdges, 2*maxEdges);
#endif
}


int obtuse_triangle(Triangle &t) {
	int i;
	Point p[3];

	p[0] = t.getVertex(0);
	p[1] = t.getVertex(1);
	p[2] = t.getVertex(2);

	for(i=0; i<3; i++) {
		if (fabs(angle(p[i], p[(i+1)%3], p[(i+2)%3])) > 3.1415926536/2.0) {
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

	ncerr = nc_put_att_text(ncid, NC_GLOBAL, "on_a_sphere", 16, "NO				  ");
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
