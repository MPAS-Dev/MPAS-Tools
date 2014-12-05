#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <netcdfcpp.h>
#include <float.h>
#include <sstream>
#include <limits>

#define ID_LEN 10

using namespace std;
//using namespace tr1;

int nCells, nVertices, vertexDegree;
bool spherical=false;
double sphereRadius=1.0;
int connectivityBase;
string in_history = "";
string in_file_id = "";

// Connectivity and location information {{{

vector<double> xCell, yCell, zCell;
vector<double> xVertex, yVertex, zVertex;
vector< vector<int> > cellsOnVertex;
vector<double> meshDensity;

// }}}

// Iterators {{{
vector<int>::iterator int_itr;
vector< vector<int> >::iterator vec_int_itr;
vector< vector<double> >::iterator vec_dbl_itr;
vector<double>::iterator dbl_itr;
// }}}

/* Building/Ordering functions {{{ */
int readGridInput(const string inputFilename);
int buildVertices();
/*}}}*/

/* Output functions {{{*/
int outputGridDimensions(const string outputFilename);
int outputGridAttributes(const string outputFilename, const string inputFilename);
int outputGridCoordinates(const string outputFilename);
int outputVertexConnectivity(const string outputFilename);
int outputMeshDensity(const string outputFilename);
/*}}}*/

/* Utility functions {{{*/
int circumcenter(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double *cx, double *cy, double *cz);
int isCCW(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3);
/*}}}*/

string gen_random(const int len);

int main ( int argc, char *argv[] ) {
	int error;
	ostringstream out_name_stream;
	string out_name;
	string in_name = "grid.nc";

	cout << endl << endl;
	cout << "************************************************************" << endl;
	cout << "ASCII_TO_NETCDF_PACKAGER:\n";
	cout << "  C++ version\n";
	cout << "  Convert a set of ascii files describing a grid into a NetCDF file describing the same grid.\n";
	cout << "  Requires cell information, and connectivity of the dual grid. Along with density values of each cell.\n";
	cout << endl << endl;
	cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
	cout << "************************************************************" << endl;
	cout << endl << endl;

	srand(time(NULL));

	cout << "Reading input grid." << endl;
	error = readGridInput(in_name);
	if(error) return 1;


	if ( argc > 1 ) { 
		out_name_stream << "grid." << argv[1] << "." << nCells << ".nc";
	} else {
		out_name_stream << "grid." << nCells << ".nc";
	}
	out_name = out_name_stream.str();

	cout << "Building veritces." << endl;
	error = buildVertices();
	if(error) return 1;

	cout << endl << "Writing file: " << out_name << endl << endl;

	cout << "Writing grid dimensions" << endl;
	if(error = outputGridDimensions(out_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Writing grid attributes" << endl;
	if(error = outputGridAttributes(out_name, in_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Writing grid coordinates" << endl;
	if(error = outputGridCoordinates(out_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Writing vertex connectivity" << endl;
	if(error = outputVertexConnectivity(out_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Reading and writing meshDensity" << endl;
	if(error = outputMeshDensity(out_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}

	return 0;
}

/* Building/Ordering functions {{{ */
int readGridInput(const string inputFilename){/*{{{*/
	double x, y, z;
	ifstream cells("end_points.dat");
	ifstream dual_cells("triangles.dat");
	ifstream density("point_density.dat");
	string line;
	vector<int> *dual_cell;
	int iVtx;

	double xRange[2], yRange[2], zRange[2];

#ifdef _DEBUG
	cout << endl << endl << "Begin function: readGridInput" << endl << endl;
#endif

	xCell.clear();
	yCell.clear();
	zCell.clear();

	xRange[0] = DBL_MAX;
	xRange[1] = DBL_MIN;
	yRange[0] = DBL_MAX;
	yRange[1] = DBL_MIN;
	zRange[0] = DBL_MAX;
	zRange[1] = DBL_MIN;

	while(!cells.eof()){
		cells >> x >> y >> z;

		if(cells.good()){
			xRange[0] = min(xRange[0], x);
			xRange[1] = max(xRange[1], x);
			yRange[0] = min(yRange[0], y);
			yRange[1] = max(yRange[1], y);
			zRange[0] = min(zRange[0], z);
			zRange[1] = max(zRange[1], z);
			xCell.push_back(x);
			yCell.push_back(y);
			zCell.push_back(z);
		}
	}
	cells.close();

	if( fabs(xRange[1] - xRange[0]) > FLT_EPSILON && fabs(yRange[1] - yRange[0]) > FLT_EPSILON && fabs(zRange[1] - zRange[0]) > FLT_EPSILON ){
		spherical = true;
	}

	if (spherical) {
		sphereRadius = sqrt(xCell[0]*xCell[0] + yCell[0]*yCell[0] + zCell[0]*zCell[0]);
	}
	
	cellsOnVertex.clear();
	connectivityBase = INT_MAX;
	
	nVertices = 0;

	for(std::string line; getline(dual_cells, line); ){
		nVertices++;
	}
	cellsOnVertex.resize(nVertices);

	dual_cells.close();
	dual_cells.open("triangles.dat");

	iVtx = 0;
	for(std::string line; getline(dual_cells, line); ){
		int start_idx = 0;
		int count = 0;
		for(int i = 0; i < line.length(); i++){
			count++;
			if(line[i] == ' '){
				std::string idx = line.substr(start_idx, count);

				cellsOnVertex.at(iVtx).push_back( atoi(idx.c_str()) );

				if (atoi(idx.c_str()) >= 0){
					connectivityBase = min(connectivityBase, atoi(idx.c_str()));
				}

				count = 0;
				start_idx = i;
			}
		}

		std::string last_idx = line.substr(start_idx);
		cellsOnVertex.at(iVtx).push_back( atoi(last_idx.c_str()) );

		if (atoi(last_idx.c_str()) >= 0){
			connectivityBase = min(connectivityBase, atoi(last_idx.c_str()));
		}

		vertexDegree = cellsOnVertex.at(iVtx).size();
		iVtx++;
	}
	dual_cells.close();

	meshDensity.clear();
	while(!density.eof()){
		double dens;

		density >> dens;
		meshDensity.push_back(dens);
	}
	density.close();

	nCells = xCell.size();
	nVertices = cellsOnVertex.size();

	cout << "Read dimensions:" << endl;
	cout << "    nCells = " << xCell.size() << endl;
	cout << "    nVertices = " << cellsOnVertex.size() << endl;
	cout << "    vertexDegree = " << vertexDegree << endl;
	cout << "    Spherical? = " << spherical << endl;
	cout << "    Sphere Radius = " << sphereRadius << endl;
	cout << "    Connectivity base = " << connectivityBase << endl;

	cout << "" << endl;
	cout << "X range: " << xRange[0] << " " << xRange[1] << endl;
	cout << "Y range: " << yRange[0] << " " << yRange[1] << endl;
	cout << "Z range: " << zRange[0] << " " << zRange[1] << endl;

	return 0;
}/*}}}*/

int buildVertices(){/*{{{*/
	double x, y, z, norm;
	int v1, v2, v3;

	xVertex.clear();
	yVertex.clear();
	zVertex.clear();

	for(int i = 0; i < cellsOnVertex.size(); i++){
		v1 = cellsOnVertex.at(i).at(0) - connectivityBase;
		v2 = cellsOnVertex.at(i).at(1) - connectivityBase;
		v3 = cellsOnVertex.at(i).at(2) - connectivityBase;

		if(!isCCW(xCell[v1], yCell[v1], zCell[v1], xCell[v2], yCell[v2], zCell[v2], xCell[v3], yCell[v3], zCell[v3])){
			v2 = cellsOnVertex.at(i).at(2) - connectivityBase;
			v3 = cellsOnVertex.at(i).at(1) - connectivityBase;
		}

		/*
		cout << "Circumcenter of : " << v1 << " " << v2 << " " << v3 << endl;
		cout << "    1 - " << xCell[v1] << " " << yCell[v1] << " " << zCell[v1] << endl;
		cout << "    2 - " << xCell[v2] << " " << yCell[v2] << " " << zCell[v2] << endl;
		cout << "    3 - " << xCell[v3] << " " << yCell[v3] << " " << zCell[v3] << endl;
		// */

		circumcenter(xCell[v1], yCell[v1], zCell[v1],
					 xCell[v2], yCell[v2], zCell[v2],
					 xCell[v3], yCell[v3], zCell[v3],
					 &x, &y, &z);

		if (spherical){
			norm = sqrt( x*x + y*y + z*z );
			x = (x / norm) * sphereRadius;
			y = (y / norm) * sphereRadius;
			z = (z / norm) * sphereRadius;
		}

		xVertex.push_back(x);
		yVertex.push_back(y);
		zVertex.push_back(z);
	}

	return 0;

}/*}}}*/
/*}}}*/

/* Output functions {{{*/
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

	int junk;

	nCells = xCell.size();

	/*
	for(vec_int_itr = edgesOnCell.begin(); vec_int_itr != edgesOnCell.end(); ++vec_int_itr){
		maxEdges = std::max(maxEdges, (int)(*vec_int_itr).size());	
	}*/
	
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
	NcDim *timeDim;
	
	// write dimensions
	if (!(nCellsDim =		grid.add_dim(	"nCells",		xCell.size())		)) return NC_ERR;
	if (!(nVerticesDim =	grid.add_dim(	"nVertices",	xVertex.size())	)) return NC_ERR;
	if (!(TWODim =			grid.add_dim(	"TWO",			2)					)) return NC_ERR;
	if (!(vertexDegreeDim = grid.add_dim(   "vertexDegree", vertexDegree)		)) return NC_ERR;
	if (!(timeDim = 		grid.add_dim(   "Time")								)) return NC_ERR;

	grid.close();
	
	// file closed when file obj goes out of scope
	return 0;
}/*}}}*/
int outputGridAttributes( const string outputFilename, const string inputFilename ){/*{{{*/
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
	NcBool history, id, spec, conventions, source, periodic;
	string history_str = "";
	string id_str = "";
	string parent_str ="";
	
	// write attributes
	if(!spherical){
		if (!(sphereAtt = grid.add_att(   "on_a_sphere", "NO\0"))) return NC_ERR;
		if (!(radiusAtt = grid.add_att(   "sphere_radius", 1.0))) return NC_ERR;
	} else {
		if (!(sphereAtt = grid.add_att(   "on_a_sphere", "YES\0"))) return NC_ERR;
		if (!(radiusAtt = grid.add_att(   "sphere_radius", sphereRadius))) return NC_ERR;
	}

	history_str += "AsciiToNetCDFPackager.x ";
	if(in_history != ""){
		history_str += "\n";
		history_str += in_history;
	}

	id_str = gen_random(ID_LEN);

	if (!(history = grid.add_att(   "history", history_str.c_str() ))) return NC_ERR;
	if (!(conventions = grid.add_att(   "Conventions", "MPAS" ))) return NC_ERR;
	if (!(source = grid.add_att(   "source", "MpasMeshConverter.x" ))) return NC_ERR;
	if (!(id = grid.add_att(   "file_id", id_str.c_str() ))) return NC_ERR;
	if (!(periodic = grid.add_att(   "is_periodic", "NO\0" ))) return NC_ERR;

	grid.close();
	
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
	NcDim *nVerticesDim = grid.get_dim( "nVertices" );

	int nCells = nCellsDim->size();
	int nVertices = nVerticesDim->size();

	//Define nc variables
	NcVar *xCellVar, *yCellVar, *zCellVar, *xVertexVar, *yVertexVar, *zVertexVar;

	int i;

	// Build and write cell coordinate arrays
	cout << "Writing xcell" << endl;
	if (!(xCellVar = grid.add_var("xCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!xCellVar->put(&xCell[0],nCells)) return NC_ERR;
	cout << "Writing ycell" << endl;
	if (!(yCellVar = grid.add_var("yCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!yCellVar->put(&yCell[0],nCells)) return NC_ERR;
	cout << "Writing zcell" << endl;
	if (!(zCellVar = grid.add_var("zCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!zCellVar->put(&zCell[0],nCells)) return NC_ERR;
	
	//Build and write vertex coordinate arrays
	cout << "Writing xvertex" << endl;
	if (!(xVertexVar = grid.add_var("xVertex", ncDouble, nVerticesDim))) return NC_ERR;
	if (!xVertexVar->put(&xVertex[0],nVertices)) return NC_ERR;
	cout << "Writing yvertex" << endl;
	if (!(yVertexVar = grid.add_var("yVertex", ncDouble, nVerticesDim))) return NC_ERR;
	if (!yVertexVar->put(&yVertex[0],nVertices)) return NC_ERR;
	cout << "Writing zvertex" << endl;
	if (!(zVertexVar = grid.add_var("zVertex", ncDouble, nVerticesDim))) return NC_ERR;
	if (!zVertexVar->put(&zVertex[0],nVertices)) return NC_ERR;

	grid.close();

	return 0;
}/*}}}*/
int outputVertexConnectivity( const string outputFilename) {/*{{{*/
	/*****************************************************************
	 *
	 * This function writes all of the *OnVertex arrays. Including
	 * cellsOnVertex
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
			tmp_arr[i*vertexDegree + j] = (*int_itr) - connectivityBase + 1;	
			j++;
		}
		i++;
	}

	if (!(covVar = grid.add_var("cellsOnVertex", ncInt, nVerticesDim, vertexDegreeDim))) return NC_ERR;
	if (!covVar->put(tmp_arr,nVertices,vertexDegree)) return NC_ERR;

	cellsOnVertex.clear();

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

	//Write meshDensity
	if (!(cDensVar = grid.add_var("meshDensity", ncDouble, nCellsDim))) return NC_ERR;
	if (!cDensVar->put(&meshDensity.at(0),nCells)) return NC_ERR;

	return 0;
}/*}}}*/
/*}}}*/

string gen_random(const int len) {/*{{{*/
	static const char alphanum[] =
		"0123456789"
//		"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		"abcdefghijklmnopqrstuvwxyz";

	string rand_str = "";

	for (int i = 0; i < len; ++i) {
		rand_str += alphanum[rand() % (sizeof(alphanum) - 1)];
	}

	return rand_str;
}/*}}}*/

int circumcenter(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double *cx, double *cy, double *cz){/*{{{*/

	if(spherical){
		double a, b, c, pbc, apc, abp;
		double bottom;
		double x23, y23, z23;
		double x31, y31, z31;
		double x12, y12, z12;

		x23 = x2 - x3;
		y23 = y2 - y3;
		z23 = z2 - z3;

		x31 = x3 - x1;
		y31 = y3 - y1;
		z31 = z3 - z1;

		x12 = x1 - x2;
		y12 = y1 - y2;
		z12 = z1 - z2;

		a = pow(x23, 2) + pow(y23, 2) + pow(z23, 2);
		b = pow(x31, 2) + pow(y31, 2) + pow(z31, 2);
		c = pow(x12, 2) + pow(y12, 2) + pow(z12, 2);
//		cout << " ABC: " << a << " " << b << " " << c << endl;

		pbc = a*(-a + b + c);
		apc = b*( a - b + c);
		abp = c*( a + b - c);

		bottom = pbc + apc + abp;

		*cx = (pbc * x1 + apc * x2 + abp * x3) / bottom;
		*cy = (pbc * y1 + apc * y2 + abp * y3) / bottom;
		*cz = (pbc * z1 + apc * z2 + abp * z3) / bottom;
	} else {
		double d;

		d = 2.0 * ( x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));

		*cx = (( powf(x1, 2) + powf(y1, 2) ) * (y2 - y3) + ( powf(x2, 2) + powf(y2, 2) ) * (y3 - y1) + ( powf(x3, 2) + powf(y3, 2) ) * (y1 - y2)) / d;
		*cy = (( powf(x1, 2) + powf(y1, 2) ) * (x3 - x2) + ( powf(x2, 2) + powf(y2, 2) ) * (x1 - x3) + ( powf(x3, 2) + powf(y3, 2) ) * (x2 - x1)) / d;
		*cz = 0.0;
	}

	return 0;
}/*}}}*/

int isCCW(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3){/*{{{*/
	double nx, ny, nz;
	double ux, uy, uz;
	double vx, vy, vz;
	double cx, cy, cz;
	double dot;

	if (spherical){
		nx = x1;
		ny = y1;
		nz = z1;
	} else {
		nx = 0.0;
		ny = 0.0;
		nz = 1.0;
	}

	ux = x2 - x1;
	uy = y2 - y1;
	uz = z2 - z1;
	vx = x3 - x1;
	vy = y3 - y1;
	vz = z3 - z1;

	cx = uy * vz - uz * vy;
	cy = uz * vx - ux * vz;
	cz = ux * vy - uy * vx;

	dot = cx * nx + cy * ny + cz * nz;

	if (dot > 0.0) {
		return 1;
	} else {
		return 0;
	}
}/*}}}*/
