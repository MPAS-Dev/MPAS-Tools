//**************************************************
// points-mpas.cpp
//
//  Purpose:
//   
//   points-mpas.cpp is supposed to take in a triangulation defined on a sphere and a point set defined on a sphere, and create
//   a mpas grid file out of the two.
//
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license. 
// 
//  Modified:
//
//    03 December 2010
//
//  Author:
//
//    Doug Jacobsen
//
//**************************************************


#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <iostream>
#include <fstream>
#include <tr1/unordered_set>
#include <vector>
#include <utility>
#include <math.h>
#include <assert.h>

//#include <netcdfcpp.h>
#include "netcdfcpp.h"
#include "triangulation.h"

#define SEED	3729

using namespace std;
using namespace tr1;

struct int_hasher {/*{{{*/
	  size_t operator()(const int v) const { return v; }
};/*}}}*/

int pt_type;
int tri_base;
int vert_levs;
int num_tracers;
double radius;
double eps;
int vordraw;

/*{{{*/ // Grid information, points, triangles, ccenters, edges
//unordered_set<pnt,pnt::hasher> edges;
unordered_set<pnt,pnt::edge_hasher> edges;
vector<pnt> edge_vec;
vector<pnt> points;
vector<pnt> ccenters;
vector<tri> triangles;
/*}}}*/
/*{{{*/ // Real connectivity arrays
vector<vector<int> > cellsOnCell, cellsOnEdge, cellsOnVertex;
vector<vector<int> > edgesOnCell, edgesOnEdge, edgesOnVertex;
vector<vector<int> > verticesOnCell, verticesOnEdge;
vector<vector<double> > weightsOnEdge;
/*}}}*/
/*{{{*/ // Unique connectivity holders
vector<unordered_set<int, int_hasher> > cellsOnCell_u;
vector<unordered_set<int, int_hasher> > cellsOnEdge_u;
vector<unordered_set<int, int_hasher> > cellsOnVertex_u;
vector<unordered_set<int, int_hasher> > edgesOnCell_u;
vector<unordered_set<int, int_hasher> > edgesOnEdge_u;
vector<unordered_set<int, int_hasher> > edgesOnVertex_u;
vector<unordered_set<int, int_hasher> > verticesOnCell_u;
vector<unordered_set<int, int_hasher> > verticesOnEdge_u;
/*}}}*/
/*{{{*/ // Grid Parameters
vector<vector<double> > kiteAreasOnVertex;
vector<double> areaCell;
vector<double> areaTriangle;
vector<double> angleEdge;
vector<double> dcEdge;
vector<double> dvEdge;
vector<double> fCell;
vector<double> fEdge;
vector<double> fVertex;
/*}}}*/
/*{{{*/ // Iterators for STL containers
vector<vector<int> >::iterator vec_int_itr;
vector<int>::iterator int_itr;
vector<vector<double> >::iterator vec_dbl_itr;
vector<double>::iterator dbl_itr;
unordered_set<int, int_hasher>::iterator u_int_itr;
vector<unordered_set<int, int_hasher> >::iterator us_itr;
unordered_set<pnt,pnt::hasher>::iterator edge_itr;
vector<pnt>::iterator point_itr;
vector<tri>::iterator tri_itr;
/*}}}*/

/*{{{*/ // Function Declarations
void readParameters();
void readPoints();
void readTriangulation();
void triangulatePoints();
void buildConnectivityArrays();
void orderConnectivityArrays();
void makeWeightsOnEdge();
int outputGridDimensions(const string outputFilename);
int outputGridAttributes(const string outputFilename);
int outputGridCoordinates(const string outputFilename);
int outputCellConnectivity(const string outputFilename);
int outputEdgeConnectivity(const string outputFilename);
int outputVertexConnectivity(const string outputFilename);
int outputCellParameters(const string outputFilename);
int outputVertexParameters(const string outputFilename);
int outputEdgeParameters(const string outputFilename);
int outputInitialConditions(const string outputFilename);
int outputMeshDensity(const string outputFilename);
int outputVordrawArrays(const string outputFilename);
int writeGraphFile();
double coriolisParameter(const pnt &p);
/*}}}*/

int main(){
	int error_code;
	string name = "grid.nc";

	cout << "Reading in paramters" << endl;
	readParameters();

	cout << " --- Points are defined on a sphere --- Radius = " << radius << endl;

	cout << "Reading in points" << endl;
	readPoints();
	cout << "Reading in triangles" << endl;
	readTriangulation();

	cout << "Building connectivity arrays" << endl;
	buildConnectivityArrays();
	cout << "Ordering connectivity arrays" << endl;
	orderConnectivityArrays();
	cout << "Making weights on edge" << endl;
	makeWeightsOnEdge();

	cout << "Writing grid dimensions" << endl;
	if(error_code = outputGridDimensions(name)){
		cout << "Error - " << error_code << endl;
		exit(error_code);
	}
	cout << "Writing grid attributes" << endl;
	if(error_code = outputGridAttributes(name)){
		cout << "Error - " << error_code << endl;
		exit(error_code);
	}
	cout << "Writing grid coordinates" << endl;
	if(error_code = outputGridCoordinates(name)){
		cout << "Error - " << error_code << endl;
		exit(error_code);
	}
	cout << "Writing cell connectivity" << endl;
	if(error_code = outputCellConnectivity(name)){
		cout << "Error - " << error_code << endl;
		exit(error_code);
	}
	cout << "Writing edge connectivity" << endl;
	if(error_code = outputEdgeConnectivity(name)){
		cout << "Error - " << error_code << endl;
		exit(error_code);
	}
	cout << "Writing vertex connectivity" << endl;
	if(error_code = outputVertexConnectivity(name)){
		cout << "Error - " << error_code << endl;
		exit(error_code);
	}
	cout << "Writing cell parameters" << endl;
	if(error_code = outputCellParameters(name)){
		cout << "Error - " << error_code << endl;
		exit(error_code);
	}
	cout << "Writing edge parameters" << endl;
	if(error_code = outputEdgeParameters(name)){
		cout << "Error - " << error_code << endl;
		exit(error_code);
	}
	cout << "Writing vertex parameters" << endl;
	if(error_code = outputVertexParameters(name)){
		cout << "Error - " << error_code << endl;
		exit(error_code);
	}
	cout << "Making and writing initial conditions" << endl;
	if(error_code = outputInitialConditions(name)){
		cout << "Error - " << error_code << endl;
		exit(error_code);
	}
	cout << "Reading and writing meshDensity" << endl;
	if(error_code = outputMeshDensity(name)){
		cout << "Error - " << error_code << endl;
		exit(error_code);
	}
	if(vordraw){
		cout << "Writing arrays for Vordraw" << endl;
		if(error_code = outputVordrawArrays(name)){
			cout << "Error - " << error_code << endl;
			exit(error_code);
		}
	}

	cout << "Writing graph.info file" << endl;
	writeGraphFile();
	cout << points.size() << " cells." << endl;
	cout << edge_vec.size() << " edges." << endl;
	cout << ccenters.size() << " vertices." << endl;
	
	return 0;
}

void readParameters(){/*{{{*/
	ifstream params("Params");
	if(!params){
		cout << "Params file not found. Writing default, and using default values." << endl;
		radius = 1.0;
		vert_levs = 1;
		num_tracers = 1;
		eps = 0.0;
		vordraw = 0;
		params.close();

		ofstream pout("Params");
		pout << "Is the input Cartesian or Latitude-Longitude (0 - Cartesian, 1 - Lat-lon)" << endl << "0" << endl;
		pout << "Are the triangles base zero or base one? (0 - base 0, 1 - base 1)" << endl << "0" << endl;
		pout << "What is the radius of the sphere these points are defined on?" << endl << "1.0" << endl;
		pout << "How many vertical levels do you want in the output grid?" << endl << "1" << endl;
		pout << "How many tracers do you want in the output grid?" << endl << "1" << endl;
		pout << "What was the convergence criteria used to make this grid?" << endl << "0.0" << endl;
		pout << "Should this grid be vordraw compatible?" << endl << "0" << endl;
		pout.close();
	} else {
		params.ignore(10000,'\n');
		params >> pt_type;
		params.ignore(10000,'\n');
		params.ignore(10000,'\n');
		params >> tri_base;
		params.ignore(10000,'\n');
		params.ignore(10000,'\n');
		params >> radius;
		params.ignore(10000,'\n');
		params.ignore(10000,'\n');
		params >> vert_levs;
		params.ignore(10000,'\n');
		params.ignore(10000,'\n');
		params >> num_tracers;
		params.ignore(10000,'\n');
		params.ignore(10000,'\n');
		params >> eps;	
		params.ignore(10000,'\n');
		params.ignore(10000,'\n');
		params >> vordraw;
		params.close();
	}

}/*}}}*/
void readPoints(){/*{{{*/
	/******************************************************************
	 *
	 * This function reads in the point set from a file, and inserts it
	 * into a vector.
	 *
	 ******************************************************************/
	pnt p;
	double lat, lon;
	int i;
	ifstream pt_start("SaveVertices");

#ifdef _DEBUG
	cout << "Setting up points." << endl;
#endif

	i = 0;
	while(!pt_start.eof()){
		if(pt_type){
			pt_start >> lat >> lon;
			p = pntFromLatLon(lat,lon);
		} else {
			pt_start >> p;
		}
		p.idx = i;
		pt_start.ignore(10000,'\n');

		if(pt_start.good()){
			points.push_back(p);
		}
		i++;
	}
	pt_start.close();
}/*}}}*/
void readTriangulation(){/*{{{*/
	/*****************************************************************
	 *
	 * This function reads in the triangulation from a file.
	 * It computes all of the circumcenters, and edges and adds them into
	 * corresponding vectors and hash tables.
	 *
	 * A hash table is used to store the edges, to ensure only one copy
	 * of each edge is kept.
	 *
	 *****************************************************************/

	pnt a, b, c, ccent, edge;
	pair<unordered_set<pnt,pnt::hasher>::iterator,bool> out_pair;
	double jv1, jv2, jv3;
	int vi1, vi2, vi3;
	int ei1, ei2, ei3;
	int i, j, junk;
	int min_vi;

	ifstream tris("SaveTriangles");

	tri t;

	i = 0;
	j = 0;
	while(!tris.eof()){
		tris >> vi1 >> vi2 >> vi3;
		ei1 = -1;
		ei2 = -1;
		ei3 = -1;

		if(tri_base == 1){
			vi1--;
			vi2--;
			vi3--;
		}

		tris.ignore(1000,'\n');

		if(tris.good()){
			a = points.at(vi1);
			b = points.at(vi2);
			c = points.at(vi3);

			vi1 = a.idx;
			vi2 = b.idx;
			vi3 = c.idx;

			if(!isCcw(a,b,c)){
				junk = vi2;
				vi2 = vi3;
				vi3 = junk;

				b = points.at(vi2);
				c = points.at(vi3);
			}

			circumcenter(a,b,c,ccent);

			ccent.normalize();

			ccent.idx = i;

			ccenters.push_back(ccent);

			edge = (a+b)/2.0;
			edge.idx = j;
			edge.isBdry = 0;
			edge.normalize();
			if(b.idx > a.idx){
				edge.vert_idx1 = a.idx;
				edge.vert_idx2 = b.idx;
			} else {
				edge.vert_idx2 = a.idx;
				edge.vert_idx1 = b.idx;
			}

			out_pair = edges.insert(edge);
			if(out_pair.second){
				edge_vec.push_back(edge);
				j++;
			}
			ei1 = (*out_pair.first).idx;

			edge = (b+c)/2.0;
			edge.idx = j;
			edge.isBdry = 0;
			edge.normalize();
			if(c.idx > b.idx){
				edge.vert_idx1 = b.idx;
				edge.vert_idx2 = c.idx;
			} else {
				edge.vert_idx2 = b.idx;
				edge.vert_idx1 = c.idx;
			}

			out_pair = edges.insert(edge);
			if(out_pair.second){
				edge_vec.push_back(edge);
				j++;
			}
			ei2 = (*out_pair.first).idx;

			edge = (c+a)/2.0;
			edge.idx = j;
			edge.isBdry = 0;
			edge.normalize();
			if(a.idx > c.idx){
				edge.vert_idx1 = c.idx;
				edge.vert_idx2 = a.idx;
			} else {
				edge.vert_idx2 = c.idx;
				edge.vert_idx1 = a.idx;
			}

			out_pair = edges.insert(edge);
			if(out_pair.second){
				edge_vec.push_back(edge);
				j++;
			}
			ei3 = (*out_pair.first).idx; 

			t = tri(vi1,vi2,vi3,i);

			t.ei1 = ei1;
			t.ei2 = ei2;
			t.ei3 = ei3;
			triangles.push_back(t);
			i++;
		}
	}

	edges.clear();
	tris.close();
}/*}}}*/
void buildConnectivityArrays(){/*{{{*/
	/*************************************************************************
	 *
	 * This function takes the triangulation and point set previously read in
	 * from files, and computes the unique connectivity arrays for use in the
	 * ordering function.
	 *
	 * This is done by adding every item into a hash table, which takes care
	 * of duplicates for us. Later, we will order this data, and transfer it
	 * into a vector
	 *
	 *************************************************************************/
	pnt a, b, c;
	pnt edge1, edge2, edge3;
	double area_temp;
	double angle;
	int vi1, vi2, vi3, ei1, ei2, ei3;

	cellsOnCell.resize(points.size());
	edgesOnCell.resize(points.size());
	verticesOnCell.resize(points.size());
	cellsOnCell_u.resize(points.size());
	edgesOnCell_u.resize(points.size());
	verticesOnCell_u.resize(points.size());
	areaCell.resize(points.size());
	fCell.resize(points.size());

	cellsOnEdge.resize(edge_vec.size());
	edgesOnEdge.resize(edge_vec.size());
	verticesOnEdge.resize(edge_vec.size());
	cellsOnEdge_u.resize(edge_vec.size());
	edgesOnEdge_u.resize(edge_vec.size());
	verticesOnEdge_u.resize(edge_vec.size());
	angleEdge.resize(edge_vec.size());
	dcEdge.resize(edge_vec.size());
	dvEdge.resize(edge_vec.size());
	fEdge.resize(edge_vec.size());

	cellsOnVertex_u.resize(ccenters.size());
	edgesOnVertex_u.resize(ccenters.size());
	cellsOnVertex.resize(ccenters.size());
	edgesOnVertex.resize(ccenters.size());
	areaTriangle.resize(ccenters.size());
	kiteAreasOnVertex.resize(ccenters.size());
	fVertex.resize(ccenters.size());

	for(int i = 0; i < points.size(); i ++){
		areaCell[i] = 0.0;
	}

	for(int i = 0; i < ccenters.size(); i ++){
		areaTriangle[i] = 0.0;
		kiteAreasOnVertex[i].resize(3);
	}

	for(tri_itr = triangles.begin(); tri_itr != triangles.end(); ++tri_itr){
		vi1 = (*tri_itr).vi1;
		vi2 = (*tri_itr).vi2;
		vi3 = (*tri_itr).vi3;

		a = points.at(vi1);
		b = points.at(vi2);
		c = points.at(vi3);

		ei1 = (*tri_itr).ei1;
		ei2 = (*tri_itr).ei2;
		ei3 = (*tri_itr).ei3;

		edge1 = edge_vec.at(ei1);
		edge2 = edge_vec.at(ei2);
		edge3 = edge_vec.at(ei3);

		cellsOnCell_u[vi1].insert(vi2);
		cellsOnCell_u[vi1].insert(vi3);
		cellsOnCell_u[vi2].insert(vi3);
		cellsOnCell_u[vi2].insert(vi1);
		cellsOnCell_u[vi3].insert(vi1);
		cellsOnCell_u[vi3].insert(vi2);

		cellsOnEdge_u[ei1].insert(vi1);
		cellsOnEdge_u[ei1].insert(vi2);
		cellsOnEdge_u[ei2].insert(vi2);
		cellsOnEdge_u[ei2].insert(vi3);
		cellsOnEdge_u[ei3].insert(vi3);
		cellsOnEdge_u[ei3].insert(vi1);

		cellsOnVertex_u[(*tri_itr).idx].insert(vi1);
		cellsOnVertex_u[(*tri_itr).idx].insert(vi2);
		cellsOnVertex_u[(*tri_itr).idx].insert(vi3);

		edgesOnCell_u[vi1].insert(ei1);
		edgesOnCell_u[vi1].insert(ei3);
		edgesOnCell_u[vi2].insert(ei1);
		edgesOnCell_u[vi2].insert(ei2);
		edgesOnCell_u[vi3].insert(ei2);
		edgesOnCell_u[vi3].insert(ei3);

		edgesOnEdge_u[ei1].insert(ei2);
		edgesOnEdge_u[ei1].insert(ei3);
		edgesOnEdge_u[ei2].insert(ei1);
		edgesOnEdge_u[ei2].insert(ei3);
		edgesOnEdge_u[ei3].insert(ei1);
		edgesOnEdge_u[ei3].insert(ei2);

		edgesOnVertex_u[(*tri_itr).idx].insert(ei1);
		edgesOnVertex_u[(*tri_itr).idx].insert(ei2);
		edgesOnVertex_u[(*tri_itr).idx].insert(ei3);

		verticesOnCell_u[vi1].insert((*tri_itr).idx);
		verticesOnCell_u[vi2].insert((*tri_itr).idx);
		verticesOnCell_u[vi3].insert((*tri_itr).idx);

		verticesOnEdge_u[ei1].insert((*tri_itr).idx);
		verticesOnEdge_u[ei2].insert((*tri_itr).idx);
		verticesOnEdge_u[ei3].insert((*tri_itr).idx);

		//areaCell
		area_temp = triArea(points.at(vi1),ccenters.at((*tri_itr).idx),edge3);
		area_temp += triArea(points.at(vi1),edge1,ccenters.at((*tri_itr).idx));
		areaCell[vi1] += area_temp;
		area_temp = triArea(points.at(vi2),ccenters.at((*tri_itr).idx),edge1);
		area_temp += triArea(points.at(vi2),edge2,ccenters.at((*tri_itr).idx));
		areaCell[vi2] += area_temp;
		area_temp = triArea(points.at(vi3),ccenters.at((*tri_itr).idx),edge2);
		area_temp += triArea(points.at(vi3),edge3,ccenters.at((*tri_itr).idx));
		areaCell[vi3] += area_temp;

		//Coriolis Parameters for cells and vertices
		fCell[vi1] = coriolisParameter(points.at(vi1));
		fCell[vi2] = coriolisParameter(points.at(vi2));
		fCell[vi3] = coriolisParameter(points.at(vi3));
		fVertex[(*tri_itr).idx] = coriolisParameter(ccenters.at((*tri_itr).idx));
	}
}/*}}}*/
void orderConnectivityArrays(){/*{{{*/
	/******************************************************************
	 *
	 * This function takes all of the hash tables that uniquely define the connectivity
	 * arrays, and orders them to be CCW and adjacent.
	 *
	 ******************************************************************/
	int i, j, k;

	// Order cellsOnEdge and verticesOnEdge, and make angleEdge
	j = 0;
	for(i = 0; i < edge_vec.size(); i++){/*{{{*/
		int cell1, cell2, vert1, vert2, tmp;
		bool flipped = false;
		pnt u, v, cross;
		pnt np;
		pnt edge;
		double sign;

		angleEdge[i] = 0;
		// Ensure that u (cellsOnEdge2 - cellsOnEdge1) crossed with
		// v (verticesOnEdge2 - verticesOnEdge1) is a right handed
		// cross product
		//
		// The vectors u and v represent the perpendicular and parallel velocity
		// directions along an edge
		assert(cellsOnEdge_u[i].size() == 2);
		u_int_itr = cellsOnEdge_u[i].begin();	
		cell1 = (*u_int_itr);
		u_int_itr++;
		cell2 = (*u_int_itr);

		assert(verticesOnEdge_u[i].size() == 2);
		u_int_itr = verticesOnEdge_u[i].begin();
		vert1 = (*u_int_itr);
		u_int_itr++;
		vert2 = (*u_int_itr);

		if(cell2 < cell1){
			tmp = cell1;
			cell1 = cell2;
			cell2 = tmp;
		}

		flipped = flip_vertices(points.at(cell1), points.at(cell2), ccenters.at(vert1), ccenters.at(vert2));

		edge = edge_vec.at(i);
		edge = gcIntersect(points.at(cell1),points.at(cell2),ccenters.at(vert1),ccenters.at(vert2));
		edge.idx = i;
		edge.isBdry = 0;
		edge.normalize();
		edge_vec.at(i) = edge;

		dcEdge[i] = points.at(cell1).sphereDistance(points.at(cell2));
		dvEdge[i] = ccenters.at(vert1).sphereDistance(ccenters.at(vert2));

		fEdge[i] = coriolisParameter(edge);

		cellsOnEdge[i].push_back(cell1);
		cellsOnEdge[i].push_back(cell2);


		// angleEdge is either:
		// 1. The angle the positive tangential direction (v)
		//    makes with the local northward direction.
		//    or
		// 2. The angles the positive normal direction (u)
		// 	  makes with the local eastward direction.
		np = pntFromLatLon(edge.getLat()+0.05, edge.getLon());
		np.normalize();

		angleEdge[i] = (ccenters.at(vert2).getLat() - ccenters.at(vert1).getLat())/dvEdge[i];
		angleEdge[i] = std::max( std::min( angleEdge[i], 1.0), -1.0);
		angleEdge[i] = acos(angleEdge[i]);

		sign = planeAngle(edge, np, ccenters.at(vert2), edge);
		if(sign != 0.0){
			sign = sign/fabs(sign);
		} else {
			sign = 1.0;
		}
		angleEdge[i] = angleEdge[i] * sign;
		if (angleEdge[i] > M_PI) angleEdge[i] = angleEdge[i] - 2.0*M_PI;
		if (angleEdge[i] < -M_PI) angleEdge[i] = angleEdge[i] + 2.0*M_PI;


		if(flipped) {
			angleEdge[i] = angleEdge[i] + M_PI;
			if (angleEdge[i] > M_PI) angleEdge[i] = angleEdge[i] - 2.0*M_PI;
			if (angleEdge[i] < -M_PI) angleEdge[i] = angleEdge[i] + 2.0*M_PI;

			tmp = vert1;
			vert1 = vert2;
			vert2 = tmp;
		}

		verticesOnEdge[i].push_back(vert1);
		verticesOnEdge[i].push_back(vert2);

		//}
	}/*}}}*/

	//Order cellsOnVertex and edgesOnVertex, areaTriangle
	for(i = 0; i < ccenters.size(); i++){/*{{{*/

		/*
		 * Since all of the *OnVertex arrays should only have 3 items, it's easy to order CCW
		 *
		 * Then using the triangles, build kitesAreasOnVertex, and areaTrianlge
		 */
		pnt a, b, c;
		pnt edge1, edge2, edge3;;
		int c1, c2, c3;
		int e1, e2, e3;
		int swp_int;

		u_int_itr = cellsOnVertex_u[i].begin();
		c1 = (*u_int_itr);
		u_int_itr++;
		c2 = (*u_int_itr);
		u_int_itr++;
		c3 = (*u_int_itr);

		a = points.at(c1);
		b = points.at(c2);
		c = points.at(c3);

		if(!isCcw(a,b,c)){
			swp_int = c2;
			c2 = c3;
			c3 = swp_int;
		}

		cellsOnVertex[i].clear();
		cellsOnVertex[i].push_back(c1);
		cellsOnVertex[i].push_back(c2);
		cellsOnVertex[i].push_back(c3);

		u_int_itr = edgesOnVertex_u[i].begin();
		e1 = (*u_int_itr);
		u_int_itr++;
		e2 = (*u_int_itr);
		u_int_itr++;
		e3 = (*u_int_itr);

		edge1 = edge_vec.at(e1);
		edge2 = edge_vec.at(e2);
		edge3 = edge_vec.at(e3);

		if(!isCcw(edge1, edge2, edge3)){
			swp_int = e2;
			e2 = e3;
			e3 = swp_int;
		}

		edgesOnVertex[i].clear();
		edgesOnVertex[i].push_back(e1);
		edgesOnVertex[i].push_back(e2);
		edgesOnVertex[i].push_back(e3);

		areaTriangle[i] = triArea(points.at(c1),points.at(c2),points.at(c3));
	}/*}}}*/

	//Order cellsOnCell, edgesOnCell and verticesOnCell
	for(i = 0; i < points.size(); i++){/*{{{*/
		// *
		// * Since we have the *OnEdge arrays ordered correctly, and the *OnVertex arrays ordered correctly,
		// * it should be easy to order the *OnCell arrays
		// *
		pnt cell_center;
		int cur_edge;
		int vert1, vert2;
		int found;
		size_t erased;

		// Choose a starting edge on cell.
		u_int_itr = edgesOnCell_u[i].begin();
		cur_edge = (*u_int_itr);
		cell_center = points.at(i);

		while(!edgesOnCell_u[i].empty()){
			//push this edge into edgesOnCell
			edgesOnCell[i].push_back(cur_edge);
			erased = edgesOnCell_u[i].erase(cur_edge);
			if(erased != 1){
				cout << " Edge " << cur_edge << " not valid" << endl;
				cout << " On cell " << cell_center << endl;
				cout << " Available edges on cell..." << endl;
				for(u_int_itr = edgesOnCell_u[i].begin(); u_int_itr != edgesOnCell_u[i].end(); ++u_int_itr){
					cout << (*u_int_itr) << " ";
				}
				cout << endl;
				assert((int)erased == 1);
			}

			//Add the cell across the edge to cellsOnCell
			if(cellsOnEdge[cur_edge].at(0) == i){
				cellsOnCell[i].push_back(cellsOnEdge[cur_edge].at(1));
			} else {
				cellsOnCell[i].push_back(cellsOnEdge[cur_edge].at(0));
			}

			// Get the correct vertices on current edge, vert1 = starting vertex, vert2 = ending vertex
			// at least in the ccw order here
			vert1 = verticesOnEdge[cur_edge].at(0);
			vert2 = verticesOnEdge[cur_edge].at(1);

			if(!isCcw(cell_center, ccenters.at(vert1), ccenters.at(vert2))){
				vert1 = verticesOnEdge[cur_edge].at(1);
				vert2 = verticesOnEdge[cur_edge].at(0);
			}

			//Push the end vertex back into verticesOnCell (in the correct order).
			verticesOnCell[i].push_back(vert2);

			// Find the next edge by cycling over the edges connceted to the ending vertex. 
			found = 0;
			for(int_itr = edgesOnVertex[vert2].begin(); int_itr != edgesOnVertex[vert2].end(); ++int_itr){
				if((cellsOnEdge[(*int_itr)].at(0) == i || cellsOnEdge[(*int_itr)].at(1) == i) 
						&& ((*int_itr) != cur_edge)){
					cur_edge = (*int_itr);
					found = 1;
					break;
				}
			}

			if(found != 1){
				break;
			}
		}
	}/*}}}*/

	//Order edgesOnEdge
	for(i = 0; i < edge_vec.size(); i++){/*{{{*/
		/*
		 * Edges on edge should be easily built now that all of the other connectivity arrays are built
		 */
		int cell1, cell2;
		int found;

		// Get cells connected to current edge
		cell1 = cellsOnEdge[i].at(0);
		cell2 = cellsOnEdge[i].at(1);

		// Cell 1 loops
		//Find current edge on cell 1, and add all edges from there to the end of edgesOnCell to edgesOnEdge
		//Then, add all edges from the beginning of edgesOnCell to edgesOnEdge
		//Since edgesOnCell should be CCW by now, these will all be CCW as well, and will iterate from cur_edge, around cell 1, 
		//and then back ccw around cell 2
		found = 0;
		for(int_itr = edgesOnCell.at(cell1).begin(); int_itr != edgesOnCell.at(cell1).end(); ++int_itr){
			if((*int_itr) == i){
				found = 1;
			} else if (found && (*int_itr) != i){
				edgesOnEdge[i].push_back((*int_itr));
			}
		}
		assert(found == 1);
		for(int_itr = edgesOnCell.at(cell1).begin(); int_itr != edgesOnCell.at(cell1).end(); ++int_itr){
			if((*int_itr) == i){
				found = 0;
			} else if(found && (*int_itr) != i){
				edgesOnEdge[i].push_back((*int_itr));
			}
		}
		assert(found == 0);
		//Cell 2 loops
		//Do the same thing done in the cell1 loops, just in cell 2 now.
		for(int_itr = edgesOnCell.at(cell2).begin(); int_itr != edgesOnCell.at(cell2).end(); ++int_itr){
			if((*int_itr) == i){
				found = 1;
			} else if (found && (*int_itr) != i){
				edgesOnEdge[i].push_back((*int_itr));
			}
		}
		assert(found == 1);
		for(int_itr = edgesOnCell.at(cell2).begin(); int_itr != edgesOnCell.at(cell2).end(); ++int_itr){
			if((*int_itr) == i){
				found = 0;
			} else if(found && (*int_itr) != i){
				edgesOnEdge[i].push_back((*int_itr));
			}
		}
		assert(found == 0);
	}/*}}}*/

	cellsOnCell_u.clear();
	cellsOnEdge_u.clear();
	cellsOnVertex_u.clear();
	edgesOnCell_u.clear();
	edgesOnEdge_u.clear();
	edgesOnVertex_u.clear();
	verticesOnEdge_u.clear();
	verticesOnCell_u.clear();
}/*}}}*/
void makeWeightsOnEdge(){/*{{{*/
	/************************************************************
	 *
	 * This function computes weightsOnEdge based on kiteAreasOnVertex
	 * and areaCell.
	 *
	 * The weights correspond to edgesOnEdge and allow MPAS to reconstruct the
	 * edge perpendicular (previously defined as v) velocity using the edge
	 * neighbors
	 *
	 * Weight formulation is defined in J. Thurburn, et al. JCP 2009
	 * Numerical representation of geostrophic modes on arbitrarily
	 * structured C-grids
	 *
	 * I'm not entirely sure it's correct yet
	 ************************************************************/
	int i, j, k;

	weightsOnEdge.resize(edge_vec.size());

	for(i = 0; i < edge_vec.size(); i++){
		int jj;
		int cur_edge;
		int prev_edge;
		int nei;
		int cell1, cell2;
		int edge1, edge2;
		int vert1, vert2;
		int neoc1, neoc2;
		double de;
		double area;
		double sum_r;

		cell1 = cellsOnEdge[i].at(0);
		cell2 = cellsOnEdge[i].at(1);

		neoc1 = edgesOnCell[cell1].size();
		neoc2 = edgesOnCell[cell2].size();

		weightsOnEdge[i].resize(edgesOnEdge[i].size());

		sum_r = 0.0;
		prev_edge = i;
		de = dcEdge.at(i);
		jj = 0;
		for(j = 0; j < neoc1-1; j++){
			cur_edge = edgesOnEdge[i].at(jj);

			//Find the vertex that is shared between prev_edge and cur_edge
			if(verticesOnEdge[prev_edge].at(0) == verticesOnEdge[cur_edge].at(0) ||
					verticesOnEdge[prev_edge].at(0) == verticesOnEdge[cur_edge].at(1)){
				vert1 = verticesOnEdge[prev_edge].at(0);
			} else if (verticesOnEdge[prev_edge].at(1) == verticesOnEdge[cur_edge].at(0) ||
					verticesOnEdge[prev_edge].at(1) == verticesOnEdge[cur_edge].at(1)){
				vert1 = verticesOnEdge[prev_edge].at(1);
			} else {
				cout << "Edge " << prev_edge << " doesn't share a vertex with edge " << cur_edge << endl;
				exit(1);
			}
	
			//Using the vertex, cell center, prev_edge and cur_edge compute the kite area using
			//the two sub triangles, CC-PE-V and CC-V-CE
			//
			//This order is CCW
			area = triArea(points.at(cell1),edge_vec.at(prev_edge), ccenters.at(vert1));
			area += triArea(points.at(cell1), ccenters.at(vert1), edge_vec.at(cur_edge));

			for(k = 0; k < 3; k++){
				if(cellsOnVertex[vert1].at(k) == cell1){
					kiteAreasOnVertex[vert1][k] = area;
				}
			}

			//Compute running sum of area ratios for edges (using kites)
			sum_r = sum_r + area/areaCell.at(cell1);

			//Compute indicator function. -1 means inward edge normal, 1 mean outward edge normal.
			//Inward means cell center is end 0 of current edge, where as 
			//Outward mean cell center is end 1 of current edge
			if(cell1 == cellsOnEdge[cur_edge].at(0)){
				nei = 1;
			} else {
				nei = -1;
			}

			//weightsOnEdge as defined in Thuburn paper referenced above.
			//nei is indicator function, 0.5 is alpha (in equation 26)
			//sum_r is running sum of area ratios
			//dvEdge and de are the le and de terms used to scale the weights.
			weightsOnEdge[i][jj] = (0.5 - sum_r)*nei*dvEdge[cur_edge]/de;

			prev_edge = cur_edge;
			jj++;
		}

		sum_r = 0.0;
		prev_edge = i;
		for(j = 0; j < neoc2-1; j++){
			cur_edge = edgesOnEdge[i].at(jj);
			if(verticesOnEdge[prev_edge].at(0) == verticesOnEdge[cur_edge].at(0) ||
					verticesOnEdge[prev_edge].at(0) == verticesOnEdge[cur_edge].at(1)){
				vert1 = verticesOnEdge[prev_edge].at(0);
			} else if (verticesOnEdge[prev_edge].at(1) == verticesOnEdge[cur_edge].at(0) ||
					verticesOnEdge[prev_edge].at(1) == verticesOnEdge[cur_edge].at(1)){
				vert1 = verticesOnEdge[prev_edge].at(1);
			} else {
				cout << "Edge " << prev_edge << " doesn't share a vertex with edge " << cur_edge << endl;
				cout << "Edge " << prev_edge << " has vertices " << verticesOnEdge[prev_edge].at(0) << " " << verticesOnEdge[prev_edge].at(1) << endl;
				cout << "Edge " << cur_edge << " has vertices " << verticesOnEdge[cur_edge].at(0) << " " << verticesOnEdge[cur_edge].at(1) << endl;
				exit(1);
			}


			area = triArea(points.at(cell2),edge_vec.at(prev_edge), ccenters.at(vert1));
			area += triArea(points.at(cell2), ccenters.at(vert1), edge_vec.at(cur_edge));

			for(k = 0; k < 3; k++){
				if(cellsOnVertex[vert1].at(k) == cell2){
					kiteAreasOnVertex[vert1][k] = area;
				}
			}

			sum_r = sum_r + area/areaCell.at(cell2);

			if(cell2 == cellsOnEdge[cur_edge].at(0)){
				nei = -1;
			} else {
				nei = 1;
			}

			weightsOnEdge[i][jj] = (0.5 - sum_r)*nei*dvEdge[cur_edge]/de;
			prev_edge = cur_edge;
			jj++;
		}
	}
}/*}}}*/

int outputGridDimensions( const string outputFilename ){/*{{{*/
	/************************************************************************
	 *
	 * This function writes the grid dimensions to the netcdf file named
	 * outputFilename
	 *
	 * **********************************************************************/
	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	// set error behaviour (matches fortran behaviour)
	NcError err(NcError::verbose_nonfatal);
	
	// open the scvtmesh file
	NcFile grid(outputFilename.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);

	int nCells, maxEdges;
	int junk;

	nCells = points.size();

	maxEdges = 0;

	for(vec_int_itr = edgesOnCell.begin(); vec_int_itr != edgesOnCell.end(); ++vec_int_itr){
		maxEdges = std::max(maxEdges, (int)(*vec_int_itr).size());	
	}
	
	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;
	
	// define dimensions
	NcDim *nCellsDim;
	NcDim *nEdgesDim;
	NcDim *nVerticesDim;
	NcDim *maxEdgesDim;
	NcDim *maxEdges2Dim;
	NcDim *TWODim;
	NcDim *THREEDim;
	NcDim *vertexDegreeDim;
	NcDim *nVertLevelsDim;
	NcDim *nTracersDim;
	NcDim *timeDim;
	
	// write dimensions
	if (!(nCellsDim =		grid.add_dim(	"nCells",		nCells)			)) return NC_ERR;
	if (!(nEdgesDim =		grid.add_dim(	"nEdges",		(3*nCells)-6)	)) return NC_ERR;
	if (!(nVerticesDim =	grid.add_dim(	"nVertices",	(2*nCells)-4)	)) return NC_ERR;
	if (!(maxEdgesDim =		grid.add_dim(	"maxEdges",		maxEdges)		)) return NC_ERR;
	if (!(maxEdges2Dim =	grid.add_dim(	"maxEdges2",	maxEdges*2)		)) return NC_ERR;
	if (!(TWODim =			grid.add_dim(	"TWO",			2)				)) return NC_ERR;
	if (!(vertexDegreeDim = grid.add_dim(   "vertexDegree", 3)				)) return NC_ERR;
	if (!(nVertLevelsDim =  grid.add_dim(   "nVertLevels",  vert_levs)      )) return NC_ERR;
	if (!(nTracersDim =     grid.add_dim(   "nTracers", num_tracers)        )) return NC_ERR;
	if (!(timeDim = 		grid.add_dim(   "Time")					)) return NC_ERR;

	cout << " nCells --- " << nCells << endl;
	cout << " nEdges --- " << (3*nCells)-6 << " " << edge_vec.size() << endl;
	cout << " nVertices --- " << (2*nCells)-4 << " " << triangles.size() << endl;
	cout << " maxEdges --- " << maxEdges << endl;
	cout << " maxEdges2 --- " << maxEdges*2 << endl;
	cout << " nVertLevels --- " << vert_levs << endl;
	cout << " nTracers --- " << num_tracers << endl;
	
	// file closed when file obj goes out of scope
	return 0;
}/*}}}*/
int outputGridAttributes( const string outputFilename ){/*{{{*/
	/************************************************************************
	 *
	 * This function writes the grid dimensions to the netcdf file named
	 * outputFilename
	 *
	 * **********************************************************************/
	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	// set error behaviour (matches fortran behaviour)
	NcError err(NcError::verbose_nonfatal);
	
	// open the scvtmesh file
	NcFile grid(outputFilename.c_str(), NcFile::Write);

	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;
	NcBool sphereAtt, radiusAtt;
	NcBool convAtt;
	
	// write attributes
	if (!(sphereAtt = grid.add_att(   "on_a_sphere", "YES             \0"))) return NC_ERR;
	if (!(radiusAtt = grid.add_att(   "sphere_radius", radius))) return NC_ERR;
	if (!(convAtt   = grid.add_att( "eps", eps ))) return NC_ERR;
	
	// file closed when file obj goes out of scope
	return 0;
}/*}}}*/
int outputGridCoordinates( const string outputFilename) {/*{{{*/
	/************************************************************************
	 *
	 * This function writes the grid coordinates to the netcdf file named
	 * outputFilename
	 * This includes all cell centers, vertices, and edges.
	 * Both cartesian and lat,lon, as well as all of their indices
	 *
	 * **********************************************************************/
	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	// set error behaviour (matches fortran behaviour)
	NcError err(NcError::verbose_nonfatal);
	
	// open the scvtmesh file
	NcFile grid(outputFilename.c_str(), NcFile::Write);
	
	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;
	
	// fetch dimensions
	NcDim *nCellsDim = grid.get_dim( "nCells" );
	NcDim *nEdgesDim = grid.get_dim( "nEdges" );
	NcDim *nVerticesDim = grid.get_dim( "nVertices" );

	int nCells = nCellsDim->size();
	int nEdges = nEdgesDim->size();
	int nVertices = nVerticesDim->size();

	//Define nc variables
	NcVar *xCellVar, *yCellVar, *zCellVar, *xEdgeVar, *yEdgeVar, *zEdgeVar, *xVertexVar, *yVertexVar, *zVertexVar;
	NcVar *lonCellVar, *latCellVar, *lonEdgeVar, *latEdgeVar, *lonVertexVar, *latVertexVar;
	NcVar *idx2cellVar, *idx2edgeVar, *idx2vertexVar;

	int i;
	
	double *x, *y, *z, *lat, *lon;
	int *idxTo;

	// Build and write cell coordinate arrays
	x = new double[nCells];
	y = new double[nCells];
	z = new double[nCells];
	lat = new double[nCells];
	lon = new double[nCells];
	idxTo = new int[nCells];
	i = 0;
	for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
		x[i] = (*point_itr).x;
		y[i] = (*point_itr).y;
		z[i] = (*point_itr).z;
		lat[i] = (*point_itr).getLat();
		lon[i] = (*point_itr).getLon();
		idxTo[i] = (*point_itr).idx+1;

		i++;
	}
	if (!(latCellVar = grid.add_var("latCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!latCellVar->put(lat,nCells)) return NC_ERR;
	if (!(lonCellVar = grid.add_var("lonCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!lonCellVar->put(lon,nCells)) return NC_ERR;
	if (!(xCellVar = grid.add_var("xCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!xCellVar->put(x,nCells)) return NC_ERR;
	if (!(yCellVar = grid.add_var("yCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!yCellVar->put(y,nCells)) return NC_ERR;
	if (!(zCellVar = grid.add_var("zCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!zCellVar->put(z,nCells)) return NC_ERR;
	if (!(idx2cellVar = grid.add_var("indexToCellID", ncInt, nCellsDim))) return NC_ERR;
	if (!idx2cellVar->put(idxTo,nCells)) return NC_ERR;
	free(x);
	free(y);
	free(z);
	free(lat);
	free(lon);
	free(idxTo);
	
	//Build and write edge coordinate arrays
	x = new double[nEdges];
	y = new double[nEdges];
	z = new double[nEdges];
	lat = new double[nEdges];
	lon = new double[nEdges];
	idxTo = new int[nEdges];

	i = 0;
	for(point_itr = edge_vec.begin(); point_itr != edge_vec.end(); ++point_itr){
		x[i] = (*point_itr).x;
		y[i] = (*point_itr).y;
		z[i] = (*point_itr).z;
		lat[i] = (*point_itr).getLat();
		lon[i] = (*point_itr).getLon();
		idxTo[i] = (*point_itr).idx+1;

		i++;
	}
	if (!(latEdgeVar = grid.add_var("latEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!latEdgeVar->put(lat,nEdges)) return NC_ERR;
	if (!(lonEdgeVar = grid.add_var("lonEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!lonEdgeVar->put(lon,nEdges)) return NC_ERR;
	if (!(xEdgeVar = grid.add_var("xEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!xEdgeVar->put(x,nEdges)) return NC_ERR;
	if (!(yEdgeVar = grid.add_var("yEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!yEdgeVar->put(y,nEdges)) return NC_ERR;
	if (!(zEdgeVar = grid.add_var("zEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!zEdgeVar->put(z,nEdges)) return NC_ERR;
	if (!(idx2edgeVar = grid.add_var("indexToEdgeID", ncInt, nEdgesDim))) return NC_ERR;
	if (!idx2edgeVar->put(idxTo, nEdges)) return NC_ERR;
	free(x);
	free(y);
	free(z);
	free(lat);
	free(lon);
	free(idxTo);

	//Build and write vertex coordinate arrays
	x = new double[nVertices];
	y = new double[nVertices];
	z = new double[nVertices];
	lat = new double[nVertices];
	lon = new double[nVertices];
	idxTo = new int[nVertices];

	i = 0;
	for(point_itr = ccenters.begin(); point_itr != ccenters.end(); ++point_itr){
		x[i] = (*point_itr).x;
		y[i] = (*point_itr).y;
		z[i] = (*point_itr).z;
		lat[i] = (*point_itr).getLat();
		lon[i] = (*point_itr).getLon();
		idxTo[i] = (*point_itr).idx+1;

		i++;
	}
	if (!(latVertexVar = grid.add_var("latVertex", ncDouble, nVerticesDim))) return NC_ERR;
	if (!latVertexVar->put(lat,nVertices)) return NC_ERR;
	if (!(lonVertexVar = grid.add_var("lonVertex", ncDouble, nVerticesDim))) return NC_ERR;
	if (!lonVertexVar->put(lon,nVertices)) return NC_ERR;
	if (!(xVertexVar = grid.add_var("xVertex", ncDouble, nVerticesDim))) return NC_ERR;
	if (!xVertexVar->put(x,nVertices)) return NC_ERR;
	if (!(yVertexVar = grid.add_var("yVertex", ncDouble, nVerticesDim))) return NC_ERR;
	if (!yVertexVar->put(y,nVertices)) return NC_ERR;
	if (!(zVertexVar = grid.add_var("zVertex", ncDouble, nVerticesDim))) return NC_ERR;
	if (!zVertexVar->put(z,nVertices)) return NC_ERR;
	if (!(idx2vertexVar = grid.add_var("indexToVertexID", ncInt, nVerticesDim))) return NC_ERR;
	if (!idx2vertexVar->put(idxTo, nVertices)) return NC_ERR;
	free(x);
	free(y);
	free(z);
	free(lat);
	free(lon);
	free(idxTo);

	return 0;
}/*}}}*/
int outputCellConnectivity( const string outputFilename) {/*{{{*/
	/*****************************************************************
	 *
	 * This function writes all of the *OnCell arrays. Including
	 * cellsOnCell
	 * edgesOnCell
	 * verticesOnCell
	 * nEdgesonCell
	 *
	 * ***************************************************************/
	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	// set error behaviour (matches fortran behaviour)
	NcError err(NcError::verbose_nonfatal);
	
	// open the scvtmesh file
	NcFile grid(outputFilename.c_str(), NcFile::Write);
	
	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;
	
	// fetch dimensions
	NcDim *nCellsDim = grid.get_dim( "nCells" );
	NcDim *maxEdgesDim = grid.get_dim( "maxEdges" );

	int nCells = nCellsDim->size();
	int maxEdges = maxEdgesDim->size();
	int i, j;

	// define nc variables
	NcVar *cocVar, *nEocVar, *eocVar, *vocVar;

	int *tmp_arr;
	
	// Build and write COC array
	tmp_arr = new int[nCells*maxEdges];

	for(i = 0; i < nCells; i++){
		for(j = 0; j < maxEdges; j++){
			tmp_arr[i*maxEdges + j] = 0;
		}
	}

	i = 0;
	for(vec_int_itr = cellsOnCell.begin(); vec_int_itr != cellsOnCell.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*maxEdges + j] = (*int_itr) + 1;
			j++;
		}
		i++;
	}
	if (!(cocVar = grid.add_var("cellsOnCell", ncInt, nCellsDim, maxEdgesDim))) return NC_ERR;
	if (!cocVar->put(tmp_arr,nCells,maxEdges)) return NC_ERR;

	// Build and write EOC array
	for(i = 0; i < nCells; i++){
		for(j = 0; j < maxEdges; j++){
			tmp_arr[i*maxEdges + j] = 0;
		}
	}

	i = 0;
	for(vec_int_itr = edgesOnCell.begin(); vec_int_itr != edgesOnCell.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*maxEdges + j] = (*int_itr) + 1;	
			j++;
		}

		i++;
	}

	if (!(eocVar = grid.add_var("edgesOnCell", ncInt, nCellsDim, maxEdgesDim))) return NC_ERR;
	if (!eocVar->put(tmp_arr,nCells,maxEdges)) return NC_ERR;

	// Build and write VOC array 
	for(i = 0; i < nCells; i++){
		for(j = 0; j < maxEdges; j++){
			tmp_arr[i*maxEdges + j] = 0;
		}
	}

	i = 0;
	for(vec_int_itr = verticesOnCell.begin(); vec_int_itr != verticesOnCell.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*maxEdges + j] = (*int_itr) + 1;	
			j++;
		}
		i++;
	}

	if (!(vocVar = grid.add_var("verticesOnCell", ncInt, nCellsDim, maxEdgesDim))) return NC_ERR;
	if (!vocVar->put(tmp_arr,nCells,maxEdges)) return NC_ERR;
	free(tmp_arr);

	//Build and write nEOC array
	tmp_arr = new int[nCells];

	i = 0;
	for(vec_int_itr = edgesOnCell.begin(); vec_int_itr != edgesOnCell.end(); ++vec_int_itr){
		tmp_arr[i] = (*vec_int_itr).size();
		i++;
	}

	if (!(nEocVar = grid.add_var("nEdgesOnCell", ncInt, nCellsDim))) return NC_ERR;
	if (!nEocVar->put(tmp_arr,nCells)) return NC_ERR;
	verticesOnCell.clear();
	edgesOnCell.clear();

	free(tmp_arr);

	return 0;
}/*}}}*/
int outputEdgeConnectivity( const string outputFilename) {/*{{{*/
	/*****************************************************************
	 *
	 * This function writes all of the *OnEdge arrays. Including
	 * cellsOnEdge
	 * edgesOnEdge
	 * verticesOnEdge
	 * nEdgesOnEdge
	 *
	 * ***************************************************************/
	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	// set error behaviour (matches fortran behaviour)
	NcError err(NcError::verbose_nonfatal);
	
	// open the scvtmesh file
	NcFile grid(outputFilename.c_str(), NcFile::Write);
	
	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;
	
	// fetch dimensions
	NcDim *nEdgesDim = grid.get_dim( "nEdges" );
	NcDim *maxEdges2Dim = grid.get_dim( "maxEdges2" );
	NcDim *vertexDegreeDim = grid.get_dim( "vertexDegree" );
	NcDim *twoDim = grid.get_dim( "TWO" );

	// define nc variables
	NcVar *coeVar, *nEoeVar, *eoeVar, *voeVar, *bdryEdgeVar;

	int nEdges = nEdgesDim->size();
	int maxEdges2 = maxEdges2Dim->size();
	int vertexDegree = vertexDegreeDim->size();
	int two = twoDim->size();
	int i, j;

	int *tmp_arr;

	// Build and write EOE array
	tmp_arr = new int[nEdges*maxEdges2];

	for(i = 0; i < nEdges; i++){
		for(j = 0; j < maxEdges2; j++){
			tmp_arr[i*maxEdges2 + j] = 0;
		}
	}

	i = 0;
	for(vec_int_itr = edgesOnEdge.begin(); vec_int_itr != edgesOnEdge.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*maxEdges2 + j] = (*int_itr) + 1;	
			j++;
		}

		i++;
	}

	if (!(eoeVar = grid.add_var("edgesOnEdge", ncInt, nEdgesDim, maxEdges2Dim))) return NC_ERR;
	if (!eoeVar->put(tmp_arr,nEdges,maxEdges2)) return NC_ERR;
	free(tmp_arr);

	// Build and write COE array
	tmp_arr = new int[nEdges*two];
	for(i = 0; i < nEdges; i++){
		for(j = 0; j < two; j++){
			tmp_arr[i*two + j] = 0;
		}
	}
	i = 0;
	for(vec_int_itr = cellsOnEdge.begin(); vec_int_itr != cellsOnEdge.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*two + j] = (*int_itr) + 1;	
			j++;
		}

		i++;
	}

	if (!(coeVar = grid.add_var("cellsOnEdge", ncInt, nEdgesDim, twoDim))) return NC_ERR;
	if (!coeVar->put(tmp_arr,nEdges,two)) return NC_ERR;

	// Build VOE array
	i = 0;
	for(vec_int_itr = verticesOnEdge.begin(); vec_int_itr != verticesOnEdge.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*two + j] = (*int_itr) + 1;	
			j++;
		}

		i++;
	}

	if (!(voeVar = grid.add_var("verticesOnEdge", ncInt, nEdgesDim, twoDim))) return NC_ERR;
	if (!voeVar->put(tmp_arr,nEdges,two)) return NC_ERR;
	free(tmp_arr);

	// Build and write nEoe array
	tmp_arr = new int[nEdges];
	i = 0;
	for(vec_int_itr = edgesOnEdge.begin(); vec_int_itr != edgesOnEdge.end(); ++vec_int_itr){
		tmp_arr[i] = (*vec_int_itr).size();
		i++;
	}

	if (!(nEoeVar = grid.add_var("nEdgesOnEdge", ncInt, nEdgesDim))) return NC_ERR;
	if (!nEoeVar->put(tmp_arr,nEdges)) return NC_ERR;

	// Build and write bdryEdge array
	i = 0;
	for(vec_int_itr = cellsOnEdge.begin(); vec_int_itr != cellsOnEdge.end(); ++vec_int_itr){
		if((*vec_int_itr).size() != 2){
			tmp_arr[i] = 1;
		} else {
			tmp_arr[i] = 0;
		}
	}

	if (!(bdryEdgeVar = grid.add_var("boundaryEdge",ncInt, nEdgesDim))) return NC_ERR;
	if (!bdryEdgeVar->put(tmp_arr,nEdges)) return NC_ERR;
	free(tmp_arr);

	cellsOnEdge.clear();
//	verticesOnEdge.clear(); // Needed for Initial conditions.
	edgesOnEdge.clear();
	
	return 0;
}/*}}}*/
int outputVertexConnectivity( const string outputFilename) {/*{{{*/
	/*****************************************************************
	 *
	 * This function writes all of the *OnVertex arrays. Including
	 * cellsOnVertex
	 * edgesOnVertex
	 *
	 * ***************************************************************/
	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	// set error behaviour (matches fortran behaviour)
	NcError err(NcError::verbose_nonfatal);
	
	// open the scvtmesh file
	NcFile grid(outputFilename.c_str(), NcFile::Write);
	
	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;
	
	// fetch dimensions
	NcDim *nVerticesDim = grid.get_dim( "nVertices" );
	NcDim *vertexDegreeDim = grid.get_dim( "vertexDegree" );

	// define nc variables
	NcVar *covVar, *eovVar, *bdryVertVar;

	int nVertices = nVerticesDim->size();
	int vertexDegree = vertexDegreeDim->size();
	int i, j;

	int *tmp_arr;

	// Build and write COV array
	tmp_arr = new int[nVertices*vertexDegree];
	
	for(i = 0; i < nVertices; i++){
		for(j = 0; j < vertexDegree; j++){
			tmp_arr[i*vertexDegree + j] = 0;
		}
	}

	i = 0;
	for(vec_int_itr = cellsOnVertex.begin(); vec_int_itr != cellsOnVertex.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*vertexDegree + j] = (*int_itr) + 1;	
			j++;
		}
		i++;
	}

	if (!(covVar = grid.add_var("cellsOnVertex", ncInt, nVerticesDim, vertexDegreeDim))) return NC_ERR;
	if (!covVar->put(tmp_arr,nVertices,vertexDegree)) return NC_ERR;

	// Build and write EOV array
	for(i = 0; i < nVertices; i++){
		for(j = 0; j < vertexDegree; j++){
			tmp_arr[i*vertexDegree + j] = 0;
		}
	}
	i = 0;
	for(vec_int_itr = edgesOnVertex.begin(); vec_int_itr != edgesOnVertex.end(); ++vec_int_itr){
		j = 0;
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			tmp_arr[i*vertexDegree + j] = (*int_itr) + 1;	
			j++;
		}

		i++;
	}
	if (!(eovVar = grid.add_var("edgesOnVertex", ncInt, nVerticesDim, vertexDegreeDim))) return NC_ERR;
	if (!eovVar->put(tmp_arr,nVertices,vertexDegree)) return NC_ERR;
	free(tmp_arr);

	// Build and write bdryVert array
	tmp_arr = new int[nVertices];
	
	i = 0;
	for(vec_int_itr = cellsOnVertex.begin(); vec_int_itr != cellsOnVertex.end(); ++vec_int_itr){
		if((*vec_int_itr).size() == vertexDegree){
			tmp_arr[i] = 0;
		} else {
			tmp_arr[i] = 1;
		}
		i++;
	}

	if (!(bdryVertVar = grid.add_var("boundaryVertex", ncInt, nVerticesDim))) return NC_ERR;
	if (!bdryVertVar->put(tmp_arr, nVertices)) return NC_ERR;

	free(tmp_arr);

	cellsOnVertex.clear();
	edgesOnVertex.clear();

	return 0;
}/*}}}*/
int outputCellParameters( const string outputFilename) {/*{{{*/
	/*********************************************************
	 *
	 * This function writes all cell parameters, including
	 * 	areaCell
	 * 	fCell
	 *
	 * *******************************************************/
	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	// set error behaviour (matches fortran behaviour)
	NcError err(NcError::verbose_nonfatal);
	
	// open the scvtmesh file
	NcFile grid(outputFilename.c_str(), NcFile::Write);
	
	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;
	
	// fetch dimensions
	NcDim *nCellsDim = grid.get_dim( "nCells" );

	// define nc variables
	NcVar *fCellVar;
	NcVar *areacVar;

	int nCells = nCellsDim->size();
	int i, j;

	if (!(fCellVar = grid.add_var("fCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!fCellVar->put(&fCell[0],nCells)) return NC_ERR;

	if (!(areacVar = grid.add_var("areaCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!areacVar->put(&areaCell[0],nCells)) return NC_ERR;

	fCell.clear();
	areaCell.clear();

	return 0;
}/*}}}*/
int outputVertexParameters( const string outputFilename) {/*{{{*/
	/*********************************************************
	 *
	 * This function writes all vertex parameters, including
	 * 	areaTriangle
	 * 	kiteAreasOnVertex
	 * 	fVertex
	 *
	 * *******************************************************/
	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	// set error behaviour (matches fortran behaviour)
	NcError err(NcError::verbose_nonfatal);
	
	// open the scvtmesh file
	NcFile grid(outputFilename.c_str(), NcFile::Write);
	
	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;
	
	// fetch dimensions
	NcDim *nVerticesDim = grid.get_dim( "nVertices" );
	NcDim *vertexDegreeDim = grid.get_dim( "vertexDegree" );

	// define nc variables
	NcVar *fVertexVar;
	NcVar *areatVar;
	NcVar *kareaVar;

	int nVertices = nVerticesDim->size();
	int vertexDegree = vertexDegreeDim->size();
	int i, j;

	double *tmp_arr;

	// Build and write fVertex
	if (!(fVertexVar = grid.add_var("fVertex", ncDouble, nVerticesDim))) return NC_ERR;
	if (!fVertexVar->put(&fVertex[0],nVertices)) return NC_ERR;

	// Build and write areaTriangle
	if (!(areatVar = grid.add_var("areaTriangle", ncDouble, nVerticesDim))) return NC_ERR;
	if (!areatVar->put(&areaTriangle[0],nVertices)) return NC_ERR;

	// Build and write kiteAreasOnVertex
	tmp_arr = new double[nVertices*vertexDegree];
	i = 0;
	for(vec_dbl_itr = kiteAreasOnVertex.begin(); vec_dbl_itr != kiteAreasOnVertex.end(); ++vec_dbl_itr){
		j = 0;
		for(dbl_itr = (*vec_dbl_itr).begin(); dbl_itr != (*vec_dbl_itr).end(); ++dbl_itr){
			tmp_arr[i*vertexDegree + j] = (*dbl_itr);
			j++;
		}
		i++;
	}

	if (!(kareaVar = grid.add_var("kiteAreasOnVertex", ncDouble, nVerticesDim, vertexDegreeDim))) return NC_ERR;
	if (!kareaVar->put(tmp_arr,nVertices,vertexDegree)) return NC_ERR;

	free(tmp_arr);

	fVertex.clear();
	areaTriangle.clear();
	kiteAreasOnVertex.clear();

	return 0;
}/*}}}*/
int outputEdgeParameters( const string outputFilename) {/*{{{*/
	/*********************************************************
	 *
	 * This function writes all grid parameters, including
	 * 	fEdge
	 *	angleEdge
	 *	dcEdge
	 *	dvEdge
	 *	weightsOnEdge
	 *
	 * *******************************************************/
	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	// set error behaviour (matches fortran behaviour)
	NcError err(NcError::verbose_nonfatal);
	
	// open the scvtmesh file
	NcFile grid(outputFilename.c_str(), NcFile::Write);
	
	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;
	
	// fetch dimensions
	NcDim *nEdgesDim = grid.get_dim( "nEdges" );
	NcDim *maxEdges2Dim = grid.get_dim( "maxEdges2" );

	// define nc variables
	NcVar *fEdgeVar;
	NcVar *angleVar;
	NcVar *kareaVar, *dcEdgeVar, *dvEdgeVar;
	NcVar *woeVar;

	int nEdges = nEdgesDim->size();
	int maxEdges2 = maxEdges2Dim->size();
	int i, j;

	double *tmp_arr;

	// Build and write fEdges
	if (!(fEdgeVar = grid.add_var("fEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!fEdgeVar->put(&fEdge[0],nEdges)) return NC_ERR;

	//Build and write angleEdge
	if (!(angleVar = grid.add_var("angleEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!angleVar->put(&angleEdge[0],nEdges)) return NC_ERR;

	//Build and write dcEdge
	if (!(dcEdgeVar = grid.add_var("dcEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!dcEdgeVar->put(&dcEdge[0],nEdges)) return NC_ERR;

	//Build and write dvEdge
	if (!(dvEdgeVar = grid.add_var("dvEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!dvEdgeVar->put(&dvEdge[0],nEdges)) return NC_ERR;

	//Build and write weightsOnEdge
	tmp_arr = new double[nEdges*maxEdges2];

	i = 0;
	for(vec_dbl_itr = weightsOnEdge.begin(); vec_dbl_itr != weightsOnEdge.end(); ++vec_dbl_itr){
		for(j = 0; j < maxEdges2; j++){
			tmp_arr[i*maxEdges2 + j] = 0.0;
		}

		j = 0;
		for(dbl_itr = (*vec_dbl_itr).begin(); dbl_itr != (*vec_dbl_itr).end(); ++dbl_itr){
			tmp_arr[i*maxEdges2 + j] = (*dbl_itr);
			j++;
		}
		i++;
	}

	if (!(woeVar = grid.add_var("weightsOnEdge", ncDouble, nEdgesDim, maxEdges2Dim))) return NC_ERR;
	if (!woeVar->put(tmp_arr,nEdges,maxEdges2)) return NC_ERR;

	free(tmp_arr);

	fEdge.clear();
	angleEdge.clear();
	dcEdge.clear();
	weightsOnEdge.clear();

	return 0;
}/*}}}*/
int outputInitialConditions( const string outputFilename) {/*{{{*/
	/***************************************************************************
	 *
	 * This function writes all initial conditions to the grid file, including
	 * (required)
	 * 		h_s
	 * 		u
	 * 		h
	 * 		tracers
	 *
	 * Any extra initial conditions can be written as needed, or added to tracers.
	 *
	 * *************************************************************************/
	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	// set error behaviour (matches fortran behaviour)
	NcError err(NcError::verbose_nonfatal);
	
	// open the scvtmesh file
	NcFile grid(outputFilename.c_str(), NcFile::Write);
	
	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;

	cout << "******************************************* " << endl;
	cout << "* Using default initial conditions." << endl << "* For more fine grained control, edit the function outputInitialConditions." << endl;
	cout << "******************************************* " << endl;
	
	// fetch dimensions
	NcDim *nCellsDim = grid.get_dim( "nCells" );
	NcDim *nEdgesDim = grid.get_dim( "nEdges" );
	NcDim *nVerticesDim = grid.get_dim( "nVertices" );
	NcDim *nVertLevelsDim = grid.get_dim( "nVertLevels" );
	NcDim *nTracersDim = grid.get_dim( "nTracers" );
	NcDim *timeDim = grid.get_dim( "Time" );

	// define nc variables
	NcVar *uVar, *usrcVar, *hVar, *hsVar, *tracerVar;

	int nCells = nCellsDim->size();
	int nEdges = nEdgesDim->size();
	int nVertices = nVerticesDim->size();
	int nVertLevels = nVertLevelsDim->size();
	int nTracers = nTracersDim->size();
	int vert1, vert2;
	int i, j, k;

	double *tmp_arr;
	double *tmp_arr2;
	double h_s_val;

	//Build and write h and h_s
	tmp_arr = new double[nCells]; // h_s
	tmp_arr2 = new double[nCells*nVertLevels]; // h

	for(i = 0; i < nCells; i++){
		h_s_val = 0.0;
		for(j = 0; j < nVertLevels; j++){
			tmp_arr2[i*nVertLevels + j] = 1000.0;
			h_s_val -= tmp_arr2[i*nVertLevels + j];
		}
		tmp_arr[i] = h_s_val;
	}

	if (!(hsVar = grid.add_var("h_s", ncDouble, nCellsDim))) return NC_ERR;
	if (!hsVar->put(tmp_arr,nCells)) return NC_ERR;
	if (!(hVar = grid.add_var("h", ncDouble, timeDim, nCellsDim, nVertLevelsDim))) return NC_ERR;
	if (!hVar->put(tmp_arr2,1,nCells,nVertLevels)) return NC_ERR;
	free(tmp_arr);
	free(tmp_arr2);

	//Build and write u and u_src
	tmp_arr = new double[nEdges*nVertLevels]; // u
	tmp_arr2 = new double[nEdges*nVertLevels]; // u_src

	for(i = 0; i < nEdges; i++){
		for(j = 0; j < nVertLevels; j++){
			vert1 = verticesOnEdge.at(i).at(0);
			vert2 = verticesOnEdge.at(i).at(1);
			if(j == 0){
				// Setup solid body rotation
				tmp_arr[i*nVertLevels + j] = (ccenters.at(vert2).z - ccenters.at(vert1).z) / dvEdge.at(i);
				tmp_arr2[i*nVertLevels + j] = (ccenters.at(vert2).z - ccenters.at(vert1).z) / dvEdge.at(i);
			} else {
				tmp_arr2[i*nVertLevels + j] = 0.0;
			}
		}
	}

	if (!(uVar = grid.add_var("u", ncDouble, timeDim, nEdgesDim, nVertLevelsDim))) return NC_ERR;
	if (!uVar->put(tmp_arr,1,nEdges,nVertLevels)) return NC_ERR;
	if (!(usrcVar = grid.add_var("u_src", ncDouble, nEdgesDim, nVertLevelsDim))) return NC_ERR;
	if (!usrcVar->put(tmp_arr2,nEdges,nVertLevels)) return NC_ERR;
	free(tmp_arr);
	free(tmp_arr2);

	//Build and write tracers
	tmp_arr = new double[nCells*nVertLevels*nTracers];
	for(i = 0; i < nCells; i++){
		for(j = 0; j < nVertLevels; j++){
			for(k = 0; k < nTracers; k++){
				tmp_arr[i*nVertLevels + j*nTracers + k] = 0.0;
			}
		}
	}

	if (!(tracerVar = grid.add_var("tracers", ncDouble, timeDim, nCellsDim, nVertLevelsDim, nTracersDim))) return NC_ERR;
	if (!tracerVar->put(tmp_arr,1,nCells,nVertLevels,nTracers)) return NC_ERR;
	free(tmp_arr);

	return 0;
}/*}}}*/
int outputMeshDensity( const string outputFilename) {/*{{{*/
	/***************************************************************************
	 *
	 * This function writes the meshDensity variable. Read in from the file SaveDensity
	 *
	 * *************************************************************************/
	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	// set error behaviour (matches fortran behaviour)
	NcError err(NcError::verbose_nonfatal);
	
	// open the scvtmesh file
	NcFile grid(outputFilename.c_str(), NcFile::Write);
	
	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;

	// fetch dimensions
	NcDim *nCellsDim = grid.get_dim( "nCells" );

	NcVar *cDensVar;

	int nCells = nCellsDim->size();
	int i, j, k;
	int junk_int;
	double junk_dbl;

	vector<double> dbl_tmp_arr;


	//Build and write meshDensity
	dbl_tmp_arr.resize(nCells);
	ifstream celldens_in("SaveDensity");

	if(!celldens_in){
		for(i = 0 ; i < nCells; i++){
			dbl_tmp_arr.at(i) = 1;
		}
	} else {
		for(i = 0; i < nCells; i++){
			celldens_in >> dbl_tmp_arr.at(i);
			
		}
	}

	celldens_in.close();

	if (!(cDensVar = grid.add_var("meshDensity", ncDouble, nCellsDim))) return NC_ERR;
	if (!cDensVar->put(&dbl_tmp_arr.at(0),nCells)) return NC_ERR;

	return 0;
}/*}}}*/
int outputVordrawArrays( const string outputFilename) {/*{{{*/
	/***************************************************************************
	 *
	 * This function writes all of the arrays for Vordraw compatibilty to the grid file, including
	 * 	cellProxy
	 * 	cellDensity
	 * 	vertexProxy
	 * 	vertexDensity
	 * 	nBorder
	 * 	xBorder
	 * 	yBorder
	 * 	zBorder
	 *
	 * *************************************************************************/
	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	// set error behaviour (matches fortran behaviour)
	NcError err(NcError::verbose_nonfatal);
	
	// open the scvtmesh file
	NcFile grid(outputFilename.c_str(), NcFile::Write);
	
	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;

	// fetch dimensions
	NcDim *nCellsDim = grid.get_dim( "nCells" );
	NcDim *nVerticesDim = grid.get_dim( "nVertices" );

	NcDim *nBorderDim;
	if (!(nBorderDim =		grid.add_dim(	"nBorder",		1)			)) return NC_ERR;

	// define nc variables
	NcVar *cDensVar, *cProxyVar, *vDensVar, *vProxyVar;
	NcVar *xBordVar, *yBordVar, *zBordVar;

	int nCells = nCellsDim->size();
	int nVertices = nVerticesDim->size();
	int nBorder = nBorderDim->size();
	int i, j, k;
	int junk_int;
	double junk_dbl;

	vector<double> dbl_tmp_arr;
	vector<int> int_tmp_arr;
	vector<double> xBorder, yBorder, zBorder;


	//Build and write cellProxy and cellDens
	int_tmp_arr.resize(nCells);
	dbl_tmp_arr.resize(nCells);
	ifstream cellproxy_in("SaveCellProxy");
	ifstream celldens_in("SaveCellDensity");
	if(!cellproxy_in){
		for(i = 0 ; i < nCells; i++){
			int_tmp_arr.at(i) = 1;
		}
	} else {
		for(i = 0; i < nCells; i++){
			cellproxy_in >> int_tmp_arr.at(i);
			
		}
	}

	if(!celldens_in){
		for(i = 0 ; i < nCells; i++){
			dbl_tmp_arr.at(i) = 1;
		}
	} else {
		for(i = 0; i < nCells; i++){
			cellproxy_in >> dbl_tmp_arr.at(i);
			
		}
	}

	cellproxy_in.close();
	celldens_in.close();

	if (!(cDensVar = grid.add_var("cellDensity", ncDouble, nCellsDim))) return NC_ERR;
	if (!cDensVar->put(&dbl_tmp_arr.at(0),nCells)) return NC_ERR;
	if (!(cProxyVar = grid.add_var("cellProxy", ncInt, nCellsDim))) return NC_ERR;
	if (!cProxyVar->put(&int_tmp_arr.at(0),nCells)) return NC_ERR;

	//Build and write vertProxy and vertDens
	int_tmp_arr.clear();
	dbl_tmp_arr.clear();
	int_tmp_arr.resize(nVertices);
	dbl_tmp_arr.resize(nVertices);
	ifstream vertproxy_in("SaveVertexProxy");
	ifstream vertdens_in("SaveVertexDensity");

	if(!vertproxy_in){
		for(i = 0 ; i < nVertices; i++){
			int_tmp_arr.at(i) = 1;
		}
	} else {
		for(i = 0; i < nVertices; i++){
			vertproxy_in >> int_tmp_arr.at(i);
			
		}
	}

	if(!vertdens_in){
		for(i = 0 ; i < nVertices; i++){
			dbl_tmp_arr.at(i) = 1;
		}
	} else {
		for(i = 0; i < nVertices; i++){
			vertproxy_in >> dbl_tmp_arr.at(i);
			
		}
	}

	vertproxy_in.close();
	vertdens_in.close();
	if (!(vDensVar = grid.add_var("vertexDensity", ncDouble, nVerticesDim))) return NC_ERR;
	if (!vDensVar->put(&dbl_tmp_arr.at(0),nVertices)) return NC_ERR;
	if (!(vProxyVar = grid.add_var("vertexProxy", ncInt, nVerticesDim))) return NC_ERR;
	if (!vProxyVar->put(&int_tmp_arr.at(0),nVertices)) return NC_ERR;

	int_tmp_arr.clear();
	dbl_tmp_arr.clear();
	xBorder.push_back(0.0);
	yBorder.push_back(0.0);
	zBorder.push_back(0.0);

	if (!(xBordVar = grid.add_var("xBorder", ncDouble, nBorderDim))) return NC_ERR;
	if (!xBordVar->put(&xBorder.at(0),nBorder)) return NC_ERR;

	if (!(yBordVar = grid.add_var("yBorder", ncDouble, nBorderDim))) return NC_ERR;
	if (!yBordVar->put(&yBorder.at(0),nBorder)) return NC_ERR;

	if (!(zBordVar = grid.add_var("zBorder", ncDouble, nBorderDim))) return NC_ERR;
	if (!zBordVar->put(&zBorder.at(0),nBorder)) return NC_ERR;

	return 0;
}/*}}}*/

int writeGraphFile() {/*{{{*/
	//This function writes out the graph.info file, for use with metis domain decomposition software.
	ofstream graph("graph.info");

	graph << points.size() << " " << edge_vec.size() << endl;

	for(vec_int_itr = cellsOnCell.begin(); vec_int_itr != cellsOnCell.end(); ++vec_int_itr){
		for(int_itr = (*vec_int_itr).begin(); int_itr != (*vec_int_itr).end(); ++int_itr){
			graph << (*int_itr)+1 << " ";
		}
		graph << endl;
	}

	graph.close();
	cellsOnCell.clear();
	return 0;
}/*}}}*/

double coriolisParameter(const pnt &p){/*{{{*/
	/******************************************************
	 * This function returns the Coriolis parameter at a point.
	 ******************************************************/
	return 2.0*7.292E-5*sin(p.getLat());
}/*}}}*/

