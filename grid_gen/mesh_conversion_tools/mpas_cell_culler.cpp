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
#include <string.h>

#include "netcdf_utils.h"

#define ID_LEN 10

using namespace std;

enum { merge, invert, preserve };

int nCells, nVertices, nEdges, vertexDegree, maxEdges;
bool spherical, periodic;
bool cullMasks = false;
double sphere_radius, xPeriod, yPeriod;
string in_history = "";
string in_file_id = "";
string in_parent_id = "";
double in_mesh_spec = 1.0;
bool outputMap = false;

// Connectivity and location information {{{

vector<int> cullCell;
vector<int> cellMap;
vector<int> vertexMap;
vector<int> edgeMap;
vector< vector<int> > verticesOnEdge;
vector< vector<int> > cellsOnEdge;
vector< vector<int> > cellsOnVertex;
vector<double> areaCell;

// }}}

/* Input/Marking Functions {{{ */
int readGridInput(const string inputFilename);
int mergeCellMasks(const string masksFilename, const int maskOp);
int markCells();
int markVertices();
int markEdges();
/*}}}*/

/* Mapping/Output functions {{{*/
int outputGridDimensions(const string outputFilename);
int outputGridAttributes(const string inputFilename, const string outputFilename);
int mapAndOutputGridCoordinates(const string inputFilename, const string outputFilename);
int mapAndOutputCellFields(const string inputFilename, const string outputFilename);
int mapAndOutputEdgeFields(const string inputFilename, const string outputFilename);
int mapAndOutputVertexFields(const string inputFilename, const string outputFilename);
int outputCellMap();
/*}}}*/

void print_usage(){/*{{{*/
	cout << endl << endl;
	cout << "Usage:" << endl;
	cout << "\tMpasCellCuller.x [input_name] [output_name] [[-m/-i/-p] masks_name] [-c]" << endl;
	cout << endl;
	cout << "\t\tinput_name:" << endl;
	cout << "\t\t\tThis argument specifies the input MPAS mesh." << endl;
	cout << "\t\toutput_name:" << endl;
	cout << "\t\t\tThis argument specifies the output culled MPAS mesh." << endl;
	cout << "\t\t\tIf not specified, it defaults to culled_mesh.nc, but" << endl;
	cout << "\t\t\tit is required if additional arguments are specified." << endl;
	cout << "\t\t-m/-i/-p:" << endl;
	cout << "\t\t\tThese arguments control how a set of masks is used when" << endl;
        cout << "\t\t\tculling a mesh." << endl;
	cout << "\t\t\tThe -m argument applies a mask to cull based on (i.e." << endl;
        cout << "\t\t\twhere the mask is 1, the mesh will be culled)." << endl;
	cout << "\t\t\tThe -i argument applies the inverse mask to cull based" << endl;
        cout << "\t\t\ton (i.e. where the mask is 0, the mesh will be" << endl;
        cout << "\t\t\tculled)." << endl;
	cout << "\t\t\tThe -p argument forces any marked cells to not be" << endl;
        cout << "\t\t\tculled." << endl;
	cout << "\t\t\tIf this argument is specified, the masks_name argument" << endl;
        cout << "\t\t\tis required" << endl;
	cout << "\t\t-c:" << endl;
	cout << "\t\t\tOutput the mapping from old to new mesh (cellMap) in" << endl;
        cout << "\t\t\t\tcellMapForward.txt, " << endl;
	cout << "\t\t\tand output the reverse mapping from new to old mesh in" << endl;
        cout << "\t\t\t\tcellMapBackward.txt." << endl;
}/*}}}*/

string gen_random(const int len);

int main ( int argc, char *argv[] ) {
	int error;
	string out_name = "culled_mesh.nc";
	string in_name = "mesh.nc";
	vector<string> mask_names;
	vector<int> mask_ops;

	cout << endl << endl;
	cout << "************************************************************" << endl;
	cout << "MPAS_CELL_CULLER:\n";
	cout << "  C++ version\n";
	cout << "  Remove cells/edges/vertices from a NetCDF MPAS Mesh file. \n";
	cout << endl << endl;
	cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
	cout << "************************************************************" << endl;
	cout << endl << endl;
	//
	//  If the input file was not specified, get it now.
	//
	if ( argc <= 1 )
	{
		cout << "\n";
		cout << "MPAS_CELL_CULLER:\n";
		cout << "  Please enter the NetCDF input filename.\n";

		cin >> in_name;

		cout << "\n";
		cout << "MPAS_CELL_CULLER:\n";
		cout << "  Please enter the output NetCDF MPAS Mesh filename.\n";

		cin >> out_name;
	}
	else if (argc == 2)
	{
		in_name = argv[1];

		cout << "\n";
		cout << "MPAS_CELL_CULLER:\n";
		cout << "  Output name not specified. Using default of culled_mesh.nc\n";
	}
	else if (argc == 3)
	{
		in_name = argv[1];
		out_name = argv[2];
	}
	else if (argc >= 10)
	{
		cout << "\n";
		cout << " ERROR: Incorrect number of arguments specified. See usage statement" << endl;
		print_usage();
		exit(1);
	}
	else
	{

		cullMasks = true;
		in_name = argv[1];
		out_name = argv[2];
		bool foundOperation;

		for ( int i = 3; i < argc; i+=2 ) {
			foundOperation = false;
			if (strcmp(argv[i], "-m") == 0 ) {
				mask_ops.push_back(static_cast<int>(merge));
				foundOperation = true;
			} else if ( strcmp(argv[i], "-i") == 0 ){
				mask_ops.push_back(static_cast<int>(invert));
				foundOperation = true;
			} else if ( strcmp(argv[i], "-p") == 0 ){
				mask_ops.push_back(static_cast<int>(preserve));
				foundOperation = true;
			} else if ( strcmp(argv[i], "-c") == 0 ){
				outputMap = true;
			} else {
				cout << " ERROR: Invalid option passed on the command line " << argv[i] << ". Exiting..." << endl;
				print_usage();
				exit(1);
			}

			if (foundOperation) {
				mask_names.push_back( argv[i+1] );
			}
		}
	}


	if(out_name == in_name){
		cout << "   ERROR: Input and Output names are the same." << endl;
		return 1;
	}

	srand(time(NULL));

	cout << "Reading input grid." << endl;
	error = readGridInput(in_name);
	if(error) return 1;

	if ( cullMasks ) {
		cout << "Reading in mask information." << endl;
		for ( int i = 0; i < mask_names.size(); i++ ) {
			error = mergeCellMasks(mask_names[i], mask_ops[i]);
			if(error) return 1;
		}
	}

	cout << "Marking cells for removal." << endl;
	error = markCells();
	if(error) return 1;

	cout << "Marking vertices for removal." << endl;
	error = markVertices();
	if(error) return 1;

	cout << "Marking edges for removal." << endl;
	error = markEdges();
	if(error) return 1;

	cout << "Writing grid dimensions" << endl;
	if(error = outputGridDimensions(out_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Writing grid attributes" << endl;
	if(error = outputGridAttributes(in_name, out_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Writing grid coordinates" << endl;
	if(error = mapAndOutputGridCoordinates(in_name, out_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Mapping and writing cell fields and culled_graph.info" << endl;
	if(error = mapAndOutputCellFields(in_name, out_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Mapping and writing edge fields" << endl;
	if(error = mapAndOutputEdgeFields(in_name, out_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Mapping and writing vertex fields" << endl;
	if(error = mapAndOutputVertexFields(in_name, out_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Outputting cell map" << endl;
	if (outputMap) {
		if(error = outputCellMap()){
			cout << "Error - " << error << endl;
			exit(error);
		}
	}

	return 0;
}

int outputCellMap(){/*{{{*/

	int iCell;
	ofstream outputfileForward, outputfileBackward;

	// forwards mapping
	outputfileForward.open("cellMapForward.txt");

	for (iCell=0 ; iCell < nCells ; iCell++) {

		outputfileForward << cellMap.at(iCell) << endl;

	}

	outputfileForward.close();

	// backwards mapping
	int nCellsNew = 0;
	vector<int> cellMapBackward;

	cellMapBackward.clear();
	cellMapBackward.resize(nCells);

	for (iCell=0 ; iCell < nCells ; iCell++) {

		if (cellMap.at(iCell) >= 0) {

			cellMapBackward.at(cellMap.at(iCell)) = iCell;
			nCellsNew++;

		}

	}

	outputfileBackward.open("cellMapBackward.txt");

	for (iCell=0 ; iCell < nCellsNew ; iCell++) {

		outputfileBackward << cellMapBackward.at(iCell) << endl;

	}

	outputfileBackward.close();

	cellMapBackward.clear();

	return 0;
}/*}}}*/

/* Input/Marking Functions {{{ */
int readGridInput(const string inputFilename){/*{{{*/
	double *xcell, *ycell,*zcell;
	double *xvertex, *yvertex,*zvertex;
	int *cellsonvertex_list, *verticesonedge_list, *cellsonedge_list;

#ifdef _DEBUG
	cout << endl << endl << "Begin function: readGridInput" << endl << endl;
#endif

	nCells = netcdf_mpas_read_dim(inputFilename, "nCells");
	nVertices = netcdf_mpas_read_dim(inputFilename, "nVertices");
	nEdges = netcdf_mpas_read_dim(inputFilename, "nEdges");
	vertexDegree = netcdf_mpas_read_dim(inputFilename, "vertexDegree");
	maxEdges = netcdf_mpas_read_dim(inputFilename, "maxEdges");
#ifdef _DEBUG
	cout << "   Reading on_a_sphere" << endl;
#endif
	spherical = netcdf_mpas_read_onsphere(inputFilename);
#ifdef _DEBUG
	cout << "   Reading sphere_radius" << endl;
#endif
	sphere_radius = netcdf_mpas_read_sphereradius(inputFilename);
#ifdef _DEBUG
	cout << "   Reading history" << endl;
#endif
	in_history = netcdf_mpas_read_history(inputFilename);
#ifdef _DEBUG
	cout << "   Reading file_id" << endl;
#endif
	in_file_id = netcdf_mpas_read_fileid(inputFilename);
#ifdef _DEBUG
	cout << "   Reading parent_id" << endl;
#endif
	in_parent_id = netcdf_mpas_read_parentid(inputFilename);
#ifdef _DEBUG
	cout << "   Reading mesh_spec" << endl;
#endif
	in_mesh_spec = netcdf_mpas_read_meshspec(inputFilename);

#ifdef _DEBUG
	cout << "   Reading is_periodic" << endl;
#endif
	periodic = netcdf_mpas_read_isperiodic(inputFilename);

#ifdef _DEBUG
	cout << "   Reading x_period" << endl;
#endif
	xPeriod = netcdf_mpas_read_xperiod(inputFilename);

#ifdef _DEBUG
	cout << "   Reading y_period" << endl;
#endif
	yPeriod = netcdf_mpas_read_yperiod(inputFilename);



	cout << "Read dimensions:" << endl;
	cout << "    nCells = " << nCells << endl;
	cout << "    nVertices = " << nVertices << endl;
	cout << "    nEdges = " << nEdges << endl;
	cout << "    vertexDegree = " << vertexDegree << endl;
	cout << "    maxEdges = " << maxEdges << endl;
	cout << "    Spherical? = " << spherical << endl;
	cout << "    Periodic? = " << periodic << endl;
	if ( periodic ) {
		cout << "    x_period = " << xPeriod << endl;
		cout << "    y_period = " << yPeriod << endl;
	}

#ifdef _DEBUG
	cout << " Read areaCell" << endl;
#endif
	areaCell.clear();
	areaCell.resize(nCells);
	netcdf_mpas_read_areacell ( inputFilename, nCells, &areaCell[0] );

#ifdef _DEBUG
	cout << " Read cullCell" << endl;
#endif
	cullCell.clear();
	cullCell.resize(nCells);
	netcdf_mpas_read_cullcell ( inputFilename, nCells, &cullCell[0] );

	// Build cellsOnVertex information
	cellsonvertex_list = new int[nVertices * vertexDegree];

#ifdef _DEBUG
	cout << " Read cellsOnVertex" << endl;
#endif
	netcdf_mpas_read_cellsonvertex ( inputFilename, nVertices, vertexDegree, cellsonvertex_list );
	cellsOnVertex.resize(nVertices);

	for(int i = 0; i < nVertices; i++){
		for(int j = 0; j < vertexDegree; j++){
			// Subtract 1 to convert into base 0 (c index space).
			cellsOnVertex.at(i).push_back(cellsonvertex_list[i*vertexDegree + j] - 1);
		}
	}

	delete[] cellsonvertex_list;

	// Build verticesOnEdge information
	verticesonedge_list = new int[nEdges * 2];
	verticesOnEdge.clear();
	verticesOnEdge.resize(nEdges);

#ifdef _DEBUG
	cout << " Read verticesOnEdge" << endl;
#endif
	netcdf_mpas_read_verticesonedge ( inputFilename, nEdges, verticesonedge_list );
	for(int i = 0; i < nEdges; i++){
		// Subtract 1 to convert into base 0 (c index space).
		verticesOnEdge.at(i).push_back(verticesonedge_list[i*2] - 1);
		verticesOnEdge.at(i).push_back(verticesonedge_list[i*2+1] - 1);
	}

	delete[] verticesonedge_list;

	// Build cellsOnEdge information
	cellsonedge_list = new int[nEdges * 2];
	cellsOnEdge.clear();
	cellsOnEdge.resize(nEdges);

#ifdef _DEBUG
	cout << " Read cellsOnEdge" << endl;
#endif
	netcdf_mpas_read_cellsonedge ( inputFilename, nEdges, cellsonedge_list );
	for(int i = 0; i < nEdges; i++){
		// Subtract 1 to convert into base 0 (c index space).
		cellsOnEdge.at(i).push_back(cellsonedge_list[i*2] - 1);
		cellsOnEdge.at(i).push_back(cellsonedge_list[i*2+1] - 1);
	}

	delete[] cellsonedge_list;

	return 0;
}/*}}}*/
int mergeCellMasks(const string masksFilename, const int maskOp){/*{{{*/
	int nRegions, nTransects;
	int *regionCellMasks, *transectCellMasks, *cellSeedMask, *flattenedMask;
	int i, j;

	nRegions = netcdf_mpas_read_dim(masksFilename, "nRegions");
	nTransects = netcdf_mpas_read_dim(masksFilename, "nTransects");

	regionCellMasks = new int[nCells*nRegions];
	transectCellMasks = new int[nCells*nTransects];
	cellSeedMask = new int[nCells];
	flattenedMask = new int[nCells];

	netcdf_mpas_read_regioncellmasks(masksFilename, nCells, nRegions, regionCellMasks);
	netcdf_mpas_read_transectcellmasks(masksFilename, nCells, nTransects, transectCellMasks);
	netcdf_mpas_read_cellseedmask(masksFilename, nCells, cellSeedMask);

	for ( i = 0; i < nCells; i++){
		flattenedMask[i] = cellSeedMask[i];
		for ( j = 0; j < nRegions; j++){
			flattenedMask[i] = max(flattenedMask[i], regionCellMasks[i * nRegions + j]);
		}

		for ( j = 0; j < nTransects; j++ ) {
			flattenedMask[i] = max(flattenedMask[i], transectCellMasks[i * nTransects + j]);
		}
	}

	if ( maskOp == invert || maskOp == merge ) {
		if ( maskOp == invert ) {
			for (i = 0; i < nCells; i++){
				flattenedMask[i] = (flattenedMask[i] + 1) % 2;
			}
		}

		for ( i = 0; i < nCells; i++ ){
			cullCell[i] = max(cullCell[i], flattenedMask[i]);
		}
	} else if ( maskOp == preserve ) {
		for ( i = 0; i < nCells; i++ ) {
			if ( flattenedMask[i] && cullCell[i] ) {
				cullCell[i] = 0;
			}
		}
	}

	delete[] cellSeedMask;
	delete[] transectCellMasks;
	delete[] regionCellMasks;
	delete[] flattenedMask;

	return 0;
}/*}}}*/
int markCells(){/*{{{*/
	int new_idx;
	int cells_removed;
	cellMap.clear();
	cellMap.resize(nCells);

	new_idx = 0;
	cells_removed = 0;
	for(int iCell = 0; iCell < nCells; iCell++){
		// Remove all cells with negative area, and cells that shouldn't be in the grid.
		if(areaCell.at(iCell) < 0 || cullCell.at(iCell) == 1){
			cellMap.at(iCell) = -1;
			cells_removed++;
		} else {
			cellMap.at(iCell) = new_idx;
			new_idx++;
		}
	}

	cout << "Removing " << cells_removed << " cells." << endl;

	return 0;
}/*}}}*/
int markVertices(){/*{{{*/
	bool keep_vertex;
	int iCell, new_idx;
	int vertices_removed;

	vertexMap.clear();
	vertexMap.resize(nVertices);

	new_idx = 0;
	vertices_removed = 0;
	for(int iVertex = 0; iVertex < nVertices; iVertex++){

		// Only keep vertices that have at least one cell connected to them
		// after cell removal.
		keep_vertex = false;
		for(int j = 0; j < cellsOnVertex.at(iVertex).size(); j++){
			iCell = cellsOnVertex.at(iVertex).at(j);
			if(iCell != -1) {
				keep_vertex = keep_vertex || (cellMap.at(iCell) != -1);
			}
		}

		if(keep_vertex){
			vertexMap.at(iVertex) = new_idx;
			new_idx++;
		} else {
			vertexMap.at(iVertex) = -1;
			vertices_removed++;
		}
	}

	cout << "Removing " << vertices_removed << " vertices." << endl;

	return 0;
}/*}}}*/
int markEdges(){/*{{{*/
	bool keep_edge;
	int vertex1, vertex2;
	int cell1, cell2;
	int new_idx;
	int edges_removed;

	edgeMap.clear();
	edgeMap.resize(nEdges);

	new_idx = 0;
	edges_removed = 0;
	for(int iEdge = 0; iEdge < nEdges; iEdge++){
		vertex1 = verticesOnEdge.at(iEdge).at(0);
		vertex2 = verticesOnEdge.at(iEdge).at(1);
		cell1 = cellsOnEdge.at(iEdge).at(0);
		cell2 = cellsOnEdge.at(iEdge).at(1);

		// Only keep an edge if it has two vertices
		// after vertex removal and at least one cell
		// after cell removal.
		if(vertex2 != -1){
			keep_edge = (vertexMap.at(vertex1) != -1 && vertexMap.at(vertex2) != -1);
		} else {
			keep_edge = false;
		}

		if(cell2 != -1){
			keep_edge = keep_edge && (cellMap.at(cell1) != -1 || cellMap.at(cell2) != -1);
		} else {
			keep_edge = keep_edge && (cellMap.at(cell1) != -1);
		}

		if(keep_edge){
			edgeMap.at(iEdge) = new_idx;
			new_idx++;
		} else {
			edgeMap.at(iEdge) = -1;
			edges_removed++;
		}
	}

	cout << "Removing " << edges_removed << " edges." << endl;

	return 0;
}/*}}}*/
/*}}}*/

/* Mapping/Output functions {{{*/
int outputGridDimensions( const string outputFilename ){/*{{{*/
	/************************************************************************
	 *
	 * This function writes the grid dimensions to the netcdf file named
	 * outputFilename
	 *
	 * **********************************************************************/

	int nCellsNew, nEdgesNew, nVerticesNew;
	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	
	// set error behaviour (matches fortran behaviour)
	NcError err(NcError::verbose_nonfatal);
	
	// open the scvtmesh file
	NcFile grid(outputFilename.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);

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
	NcDim *TWODim;
	NcDim *THREEDim;
	NcDim *vertexDegreeDim;
	NcDim *timeDim;

	nCellsNew = 0;
	for(int iCell = 0; iCell < nCells; iCell++){
		nCellsNew += (cellMap.at(iCell) != -1);
	}

	nVerticesNew = 0;
	for(int iVertex = 0; iVertex < nVertices; iVertex++){
		nVerticesNew += (vertexMap.at(iVertex) != -1);
	}

	nEdgesNew = 0;
	for(int iEdge = 0; iEdge < nEdges; iEdge++){
		nEdgesNew += (edgeMap.at(iEdge) != -1);
	}
	
	// write dimensions
	if (!(nCellsDim =		grid.add_dim(	"nCells",		nCellsNew)		)) return NC_ERR;
	if (!(nEdgesDim =		grid.add_dim(	"nEdges",		nEdgesNew)		)) return NC_ERR;
	if (!(nVerticesDim =	grid.add_dim(	"nVertices",	nVerticesNew)	)) return NC_ERR;
	if (!(TWODim =			grid.add_dim(	"TWO",			2)				)) return NC_ERR;
	if (!(vertexDegreeDim = grid.add_dim(   "vertexDegree", vertexDegree)	)) return NC_ERR;
	if (!(timeDim = 		grid.add_dim(   "Time")							)) return NC_ERR;

	grid.close();
	
	// file closed when file obj goes out of scope
	return 0;
}/*}}}*/
int outputGridAttributes( const string inputFilename, const string outputFilename ){/*{{{*/
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
	NcBool sphereAtt, radiusAtt, periodicAtt, xPeriodAtt, yPeriodAtt;
	NcBool history, id, spec, conventions, source, parent_id;
	string history_str = "";
	string id_str = "";
	string parent_str = "";
	
	// write attributes
	if(!spherical){
		if (!(sphereAtt = grid.add_att(   "on_a_sphere", "NO\0"))) return NC_ERR;
		if (!(radiusAtt = grid.add_att(   "sphere_radius", 0.0))) return NC_ERR;
	} else {
		if (!(sphereAtt = grid.add_att(   "on_a_sphere", "YES\0"))) return NC_ERR;
		if (!(radiusAtt = grid.add_att(   "sphere_radius", sphere_radius))) return NC_ERR;
	}

	if(!periodic){
		if (!(periodicAtt = grid.add_att(   "is_periodic", "NO\0"))) return NC_ERR;
	} else {
		if (!(periodicAtt = grid.add_att(   "is_periodic", "YES\0"))) return NC_ERR;
		if (!(xPeriodAtt = grid.add_att(   "x_period", xPeriod))) return NC_ERR;
		if (!(xPeriodAtt = grid.add_att(   "y_period", yPeriod))) return NC_ERR;
	}

	history_str += "MpasCellCuller.x ";
	history_str += inputFilename;
	history_str += " ";
	history_str += outputFilename;
	if(in_history != ""){
		history_str += "\n";
		history_str += in_history;
	}

	if(in_file_id != "" ){
		parent_str = in_file_id;
		if(in_parent_id != ""){
			parent_str += "\n";
			parent_str += in_parent_id;
		}
		if (!(id = grid.add_att(   "parent_id", parent_str.c_str() ))) return NC_ERR;
	}

	id_str = gen_random(ID_LEN);

	if (!(history = grid.add_att(   "history", history_str.c_str() ))) return NC_ERR;
	if (!(spec = grid.add_att(   "mesh_spec", in_mesh_spec ))) return NC_ERR;
	if (!(conventions = grid.add_att(   "Conventions", "MPAS" ))) return NC_ERR;
	if (!(source = grid.add_att(   "source", "MpasCellCuller.x" ))) return NC_ERR;
	if (!(id = grid.add_att(   "file_id", id_str.c_str() ))) return NC_ERR;

	grid.close();
	
	// file closed when file obj goes out of scope
	return 0;
}/*}}}*/
int mapAndOutputGridCoordinates( const string inputFilename, const string outputFilename) {/*{{{*/
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

	int nCellsNew = nCellsDim->size();
	int nEdgesNew = nEdgesDim->size();
	int nVerticesNew = nVerticesDim->size();

	//Define nc variables
	NcVar *xCellVar, *yCellVar, *zCellVar, *xEdgeVar, *yEdgeVar, *zEdgeVar, *xVertexVar, *yVertexVar, *zVertexVar;
	NcVar *lonCellVar, *latCellVar, *lonEdgeVar, *latEdgeVar, *lonVertexVar, *latVertexVar;
	NcVar *idx2cellVar, *idx2edgeVar, *idx2vertexVar;

	int i, idx_map;
	
	double *xOld, *yOld, *zOld, *latOld, *lonOld;
	double *xNew, *yNew, *zNew, *latNew, *lonNew;
	int *idxToNew;

	// Build and write cell coordinate arrays
	xOld = new double[nCells];
	yOld = new double[nCells];
	zOld = new double[nCells];
	latOld = new double[nCells];
	lonOld = new double[nCells];

	xNew = new double[nCellsNew];
	yNew = new double[nCellsNew];
	zNew = new double[nCellsNew];
	latNew = new double[nCellsNew];
	lonNew = new double[nCellsNew];
	idxToNew = new int[nCellsNew];

	netcdf_mpas_read_xyzcell ( inputFilename, nCells, xOld, yOld, zOld );
	netcdf_mpas_read_latloncell ( inputFilename, nCells, latOld, lonOld );

	idx_map = 0;
	for(int iCell = 0; iCell < nCells; iCell++){
		if(cellMap.at(iCell) != -1){
			xNew[idx_map] = xOld[iCell];
			yNew[idx_map] = yOld[iCell];
			zNew[idx_map] = zOld[iCell];
			latNew[idx_map] = latOld[iCell];
			lonNew[idx_map] = lonOld[iCell];
			idxToNew[idx_map] = idx_map+1;
			idx_map++;
		}
	}

	if (!(latCellVar = grid.add_var("latCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!latCellVar->put(latNew,nCellsNew)) return NC_ERR;
	if (!(lonCellVar = grid.add_var("lonCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!lonCellVar->put(lonNew,nCellsNew)) return NC_ERR;
	if (!(xCellVar = grid.add_var("xCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!xCellVar->put(xNew,nCellsNew)) return NC_ERR;
	if (!(yCellVar = grid.add_var("yCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!yCellVar->put(yNew,nCellsNew)) return NC_ERR;
	if (!(zCellVar = grid.add_var("zCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!zCellVar->put(zNew,nCellsNew)) return NC_ERR;
	if (!(idx2cellVar = grid.add_var("indexToCellID", ncInt, nCellsDim))) return NC_ERR;
	if (!idx2cellVar->put(idxToNew,nCellsNew)) return NC_ERR;
	delete[] xOld;
	delete[] yOld;
	delete[] zOld;
	delete[] latOld;
	delete[] lonOld;

	delete[] xNew;
	delete[] yNew;
	delete[] zNew;
	delete[] latNew;
	delete[] lonNew;
	delete[] idxToNew;
	
	//Build and write edge coordinate arrays
	xOld = new double[nEdges];
	yOld = new double[nEdges];
	zOld = new double[nEdges];
	latOld = new double[nEdges];
	lonOld = new double[nEdges];

	xNew = new double[nEdgesNew];
	yNew = new double[nEdgesNew];
	zNew = new double[nEdgesNew];
	latNew = new double[nEdgesNew];
	lonNew = new double[nEdgesNew];
	idxToNew = new int[nEdgesNew];

	netcdf_mpas_read_xyzedge ( inputFilename, nEdges, xOld, yOld, zOld );
	netcdf_mpas_read_latlonedge ( inputFilename, nEdges, latOld, lonOld );

	idx_map = 0;
	for(int iEdge = 0; iEdge < nEdges; iEdge++){
		if(edgeMap.at(iEdge) != -1){
			xNew[idx_map] = xOld[iEdge];
			yNew[idx_map] = yOld[iEdge];
			zNew[idx_map] = zOld[iEdge];
			latNew[idx_map] = latOld[iEdge];
			lonNew[idx_map] = lonOld[iEdge];
			idxToNew[idx_map] = idx_map+1;
			idx_map++;
		}
	}

	if (!(latEdgeVar = grid.add_var("latEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!latEdgeVar->put(latNew,nEdgesNew)) return NC_ERR;
	if (!(lonEdgeVar = grid.add_var("lonEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!lonEdgeVar->put(lonNew,nEdgesNew)) return NC_ERR;
	if (!(xEdgeVar = grid.add_var("xEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!xEdgeVar->put(xNew,nEdgesNew)) return NC_ERR;
	if (!(yEdgeVar = grid.add_var("yEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!yEdgeVar->put(yNew,nEdgesNew)) return NC_ERR;
	if (!(zEdgeVar = grid.add_var("zEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!zEdgeVar->put(zNew,nEdgesNew)) return NC_ERR;
	if (!(idx2edgeVar = grid.add_var("indexToEdgeID", ncInt, nEdgesDim))) return NC_ERR;
	if (!idx2edgeVar->put(idxToNew, nEdgesNew)) return NC_ERR;
	delete[] xOld;
	delete[] yOld;
	delete[] zOld;
	delete[] latOld;
	delete[] lonOld;

	delete[] xNew;
	delete[] yNew;
	delete[] zNew;
	delete[] latNew;
	delete[] lonNew;
	delete[] idxToNew;

	//Build and write vertex coordinate arrays
	xOld = new double[nVertices];
	yOld = new double[nVertices];
	zOld = new double[nVertices];
	latOld = new double[nVertices];
	lonOld = new double[nVertices];

	xNew = new double[nVerticesNew];
	yNew = new double[nVerticesNew];
	zNew = new double[nVerticesNew];
	latNew = new double[nVerticesNew];
	lonNew = new double[nVerticesNew];
	idxToNew = new int[nVerticesNew];

	netcdf_mpas_read_xyzvertex ( inputFilename, nVertices, xOld, yOld, zOld );
	netcdf_mpas_read_latlonvertex ( inputFilename, nVertices, latOld, lonOld );

	idx_map = 0;
	for(int iVertex = 0; iVertex < nVertices; iVertex++){
		if(vertexMap.at(iVertex) != -1){
			xNew[idx_map] = xOld[iVertex];
			yNew[idx_map] = yOld[iVertex];
			zNew[idx_map] = zOld[iVertex];
			latNew[idx_map] = latOld[iVertex];
			lonNew[idx_map] = lonOld[iVertex];
			idxToNew[idx_map] = idx_map+1;
			idx_map++;
		}
	}

	if (!(latVertexVar = grid.add_var("latVertex", ncDouble, nVerticesDim))) return NC_ERR;
	if (!latVertexVar->put(latNew,nVerticesNew)) return NC_ERR;
	if (!(lonVertexVar = grid.add_var("lonVertex", ncDouble, nVerticesDim))) return NC_ERR;
	if (!lonVertexVar->put(lonNew,nVerticesNew)) return NC_ERR;
	if (!(xVertexVar = grid.add_var("xVertex", ncDouble, nVerticesDim))) return NC_ERR;
	if (!xVertexVar->put(xNew,nVerticesNew)) return NC_ERR;
	if (!(yVertexVar = grid.add_var("yVertex", ncDouble, nVerticesDim))) return NC_ERR;
	if (!yVertexVar->put(yNew,nVerticesNew)) return NC_ERR;
	if (!(zVertexVar = grid.add_var("zVertex", ncDouble, nVerticesDim))) return NC_ERR;
	if (!zVertexVar->put(zNew,nVerticesNew)) return NC_ERR;
	if (!(idx2vertexVar = grid.add_var("indexToVertexID", ncInt, nVerticesDim))) return NC_ERR;
	if (!idx2vertexVar->put(idxToNew, nVerticesNew)) return NC_ERR;
	delete[] xOld;
	delete[] yOld;
	delete[] zOld;
	delete[] latOld;
	delete[] lonOld;

	delete[] xNew;
	delete[] yNew;
	delete[] zNew;
	delete[] latNew;
	delete[] lonNew;
	delete[] idxToNew;

	grid.close();

	return 0;
}/*}}}*/
int mapAndOutputCellFields( const string inputFilename, const string outputFilename) {/*{{{*/
	/*****************************************************************
	 *
	 * This function maps and writes all of the cell related fields. Including
	 * cellsOnCell
	 * edgesOnCell
	 * verticesOnCell
	 * nEdgesonCell
	 * areaCell
	 * meshDensity
	 *
	 * It also writes the graph.info file which can be used to decompose the mesh.
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
	NcDim *nEdgesDim = grid.get_dim( "nEdges" );
	NcDim *maxEdgesDim;
	NcDim *maxEdges2Dim;

	int nCellsNew = nCellsDim->size();
	int nEdgesNew = nEdgesDim->size();
	int maxEdgesNew;
	int edgeCount;

	// define nc variables
	NcVar *cocVar, *nEocVar, *eocVar, *vocVar, *areacVar;
	NcVar *cDensVar;

	double *meshDensityOld, *meshDensityNew;
	double *areaCellNew;
	int *tmp_arr_old, *nEdgesOnCellOld, *nEdgesOnCellNew;
	int *tmp_arr_new;
	
	tmp_arr_old = new int[nCells*maxEdges];
	nEdgesOnCellOld = new int[nCells];
	nEdgesOnCellNew = new int[nCellsNew];

	netcdf_mpas_read_edgesoncell ( inputFilename, nCells, maxEdges, tmp_arr_old );
	netcdf_mpas_read_nedgesoncell ( inputFilename, nCells, nEdgesOnCellOld );

	// Need to map nEdgesOnCell to get maxEdges
	maxEdgesNew = 0;
	for(int iCell = 0; iCell < nCells; iCell++){
		if(cellMap.at(iCell) != -1){
			nEdgesOnCellNew[cellMap.at(iCell)] = nEdgesOnCellOld[iCell];
			maxEdgesNew = max(maxEdgesNew, nEdgesOnCellNew[cellMap.at(iCell)]);
		}
	}
	tmp_arr_new = new int[nCells * maxEdgesNew];

	// Write maxEdges and maxEdges2 to output file.
	if (!(maxEdgesDim =		grid.add_dim(	"maxEdges",		maxEdgesNew)			)) return NC_ERR;
	if (!(maxEdges2Dim =	grid.add_dim(	"maxEdges2",	maxEdgesNew*2)			)) return NC_ERR;

	// Write nEdgesOncell to output file
	if (!(nEocVar = grid.add_var("nEdgesOnCell", ncInt, nCellsDim))) return NC_ERR;
	if (!nEocVar->put(nEdgesOnCellNew,nCellsNew)) return NC_ERR;

	// Map edgesOnCell
	for(int iCell = 0; iCell < nCells; iCell++){
#ifdef _DEBUG
		cout << "On cell: " << iCell << endl;
#endif
		if(cellMap.at(iCell) != -1){
			for(int j = 0; j < maxEdgesNew; j++){
				int iEdge = tmp_arr_old[iCell*maxEdges + j] - 1;

				if(iEdge != -1 && iEdge < edgeMap.size()){
					tmp_arr_new[cellMap.at(iCell)*maxEdgesNew + j] = edgeMap.at(iEdge)+1;
				} else {
					tmp_arr_new[cellMap.at(iCell)*maxEdgesNew + j] = 0;
				}

#ifdef _DEBUG
				cout << "    Mapping edge: " << iEdge << " to " << tmp_arr_new[cellMap.at(iCell)*maxEdgesNew + j] << " dbg info: " << tmp_arr_old[iCell*maxEdges + j] << " " << j << endl;
#endif
			}
		}
	}

	if (!(eocVar = grid.add_var("edgesOnCell", ncInt, nCellsDim, maxEdgesDim))) return NC_ERR;
	if (!eocVar->put(tmp_arr_new,nCellsNew,maxEdgesNew)) return NC_ERR;

	netcdf_mpas_read_cellsoncell ( inputFilename, nCells, maxEdges, tmp_arr_old );

	// Map cellsOnCell, and determine number of edges in graph.
	edgeCount = 0;
	for(int iCell = 0; iCell < nCells; iCell++){
		if(cellMap.at(iCell) != -1){
			for(int j = 0; j < nEdgesOnCellOld[iCell]; j++){
				int coc = tmp_arr_old[iCell*maxEdges + j] - 1;

				if(coc != -1 && coc < nCells && cellMap.at(coc) < nCellsNew && cellMap.at(coc) != -1){
					tmp_arr_new[cellMap.at(iCell)*maxEdgesNew + j] = cellMap.at(coc)+1;
					edgeCount++;
				} else {
					tmp_arr_new[cellMap.at(iCell)*maxEdgesNew + j] = 0;
				}
			}
		}
	}
	edgeCount = edgeCount / 2;

	// Build graph.info file
	ofstream graph("culled_graph.info");
	graph << nCellsNew << " " << edgeCount << endl;
	for(int iCell = 0; iCell < nCellsNew; iCell++){
		for(int j = 0; j < nEdgesOnCellNew[iCell]; j++){
			if (tmp_arr_new[iCell * maxEdgesNew + j] != 0) {
				graph << tmp_arr_new[iCell * maxEdgesNew + j] << " ";
			}
		}
		graph << endl;
	}
	graph.close();

	if (!(cocVar = grid.add_var("cellsOnCell", ncInt, nCellsDim, maxEdgesDim))) return NC_ERR;
	if (!cocVar->put(tmp_arr_new,nCellsNew,maxEdgesNew)) return NC_ERR;
	delete[] nEdgesOnCellNew;
	delete[] nEdgesOnCellOld;

	netcdf_mpas_read_verticesoncell ( inputFilename, nCells, maxEdges, tmp_arr_old );

	// Map verticesOnCell
	for(int iCell = 0; iCell < nCells; iCell++){
		if(cellMap.at(iCell) != -1){
			for(int j = 0; j < maxEdgesNew; j++){
				int iVertex = tmp_arr_old[iCell*maxEdges + j] - 1;

				if(iVertex != -1 && iVertex < vertexMap.size()){
					tmp_arr_new[cellMap.at(iCell)*maxEdgesNew + j] = vertexMap.at(iVertex)+1;
				} else {
					tmp_arr_new[cellMap.at(iCell)*maxEdgesNew + j] = 0;
				}
			}
		}
	}

	if (!(vocVar = grid.add_var("verticesOnCell", ncInt, nCellsDim, maxEdgesDim))) return NC_ERR;
	if (!vocVar->put(tmp_arr_new,nCellsNew,maxEdgesNew)) return NC_ERR;

	delete[] tmp_arr_old;
	delete[] tmp_arr_new;

	// Map areaCell
	areaCellNew = new double[nCellsNew]; 

	for(int iCell = 0; iCell < nCells; iCell++){
		if(cellMap.at(iCell) != -1){
			areaCellNew[cellMap.at(iCell)] = areaCell.at(iCell);
		}
	}

	if (!(areacVar = grid.add_var("areaCell", ncDouble, nCellsDim))) return NC_ERR;
	if (!areacVar->put(areaCellNew,nCellsNew)) return NC_ERR;

	delete[] areaCellNew;

	// Map meshDensity
	meshDensityOld = new double[nCells];
	meshDensityNew = new double[nCellsNew];

	netcdf_mpas_read_mesh_density ( inputFilename, nCells, meshDensityOld );

	for(int iCell = 0; iCell < nCells; iCell++){
		if(cellMap.at(iCell) != -1){
			meshDensityNew[cellMap.at(iCell)] = meshDensityOld[iCell];
		}
	}


	//Write meshDensity
	if (!(cDensVar = grid.add_var("meshDensity", ncDouble, nCellsDim))) return NC_ERR;
	if (!cDensVar->put(meshDensityNew,nCellsNew)) return NC_ERR;

	return 0;
}/*}}}*/
int mapAndOutputEdgeFields( const string inputFilename, const string outputFilename) {/*{{{*/
	/*****************************************************************
	 *
	 * This function maps and writes all of the edge related fields. Including
	 * cellsOnEdge
	 * edgesOnEdge
	 * verticesOnEdge
	 * nEdgesOnEdge
	 * weightsOnEdge
	 * dvEdge
	 * dcEdge
	 * angleEdge
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
	NcVar *coeVar, *nEoeVar, *eoeVar, *voeVar, *woeVar;
	NcVar *angleVar;
	NcVar *dcEdgeVar, *dvEdgeVar;

	int nEdgesNew = nEdgesDim->size();
	int maxEdges2New = maxEdges2Dim->size();
	int vertexDegree = vertexDegreeDim->size();
	int two = twoDim->size();

	int *nEdgesOnEdgeOld, *nEdgesOnEdgeNew;
	int *edgesOnEdgeOld, *edgesOnEdgeNew;
	int *cellsOnEdgeOld, *cellsOnEdgeNew;
	int *verticesOnEdgeOld, *verticesOnEdgeNew;
	double *weightsOnEdgeOld, *weightsOnEdgeNew;
	double *dvEdgeOld, *dvEdgeNew;
	double *dcEdgeOld, *dcEdgeNew;
	double *angleEdgeOld, *angleEdgeNew;

	// Need to map cellsOnEdge and verticesOnEdge
	cellsOnEdgeOld = new int[nEdges*2];
	cellsOnEdgeNew = new int[nEdgesNew*2];
	verticesOnEdgeOld = new int[nEdges*2];
	verticesOnEdgeNew = new int[nEdgesNew*2];

	netcdf_mpas_read_cellsonedge( inputFilename, nEdges, cellsOnEdgeOld);
	netcdf_mpas_read_verticesonedge( inputFilename, nEdges, verticesOnEdgeOld);

	// Map cellsOnEdge and verticesOnEdge
	for(int iEdge = 0; iEdge < nEdges; iEdge++){
		if(edgeMap.at(iEdge) != -1){
			int cell1, cell2;
			int vertex1, vertex2;

			cell1 = cellsOnEdgeOld[iEdge * 2] - 1;
			cell2 = cellsOnEdgeOld[iEdge * 2 + 1] - 1;
			vertex1 = verticesOnEdgeOld[iEdge * 2] - 1;
			vertex2 = verticesOnEdgeOld[iEdge * 2 + 1] - 1;

#ifdef _DEBUG
			cout << "Defining edge: " << endl;
			cout << "   Old cell1: " << cell1 << endl;
			cout << "   Old cell2: " << cell2 << endl;
			cout << "   Old vertex1: " << vertex1 << endl;
			cout << "   Old vertex2: " << vertex2 << endl;
#endif

			if(cell1 != -1 && cell2 != -1){
				if(cellMap.at(cell1) != -1 && cellMap.at(cell2) != -1){
					cellsOnEdgeNew[edgeMap.at(iEdge)*2] = cellMap.at(cell1) + 1;
					cellsOnEdgeNew[edgeMap.at(iEdge)*2+1] = cellMap.at(cell2) + 1;

					verticesOnEdgeNew[edgeMap.at(iEdge)*2] = vertexMap.at(vertex1) + 1;
					verticesOnEdgeNew[edgeMap.at(iEdge)*2+1] = vertexMap.at(vertex2) + 1;
				} else if (cellMap.at(cell2) == -1){
					cellsOnEdgeNew[edgeMap.at(iEdge)*2] = cellMap.at(cell1) + 1;
					cellsOnEdgeNew[edgeMap.at(iEdge)*2+1] = 0;

					verticesOnEdgeNew[edgeMap.at(iEdge)*2] = vertexMap.at(vertex1) + 1;
					verticesOnEdgeNew[edgeMap.at(iEdge)*2+1] = vertexMap.at(vertex2) + 1;
				} else if (cellMap.at(cell1) == -1){
					cellsOnEdgeNew[edgeMap.at(iEdge)*2] = cellMap.at(cell2) + 1;
					cellsOnEdgeNew[edgeMap.at(iEdge)*2+1] = 0;

					verticesOnEdgeNew[edgeMap.at(iEdge)*2] = vertexMap.at(vertex2) + 1;
					verticesOnEdgeNew[edgeMap.at(iEdge)*2+1] = vertexMap.at(vertex1) + 1;
				}
			} else if(cell2 == -1){
				if(cellMap.at(cell1) != -1){
					cellsOnEdgeNew[edgeMap.at(iEdge)*2] = cellMap.at(cell1) + 1;
					cellsOnEdgeNew[edgeMap.at(iEdge)*2+1] = 0;

					verticesOnEdgeNew[edgeMap.at(iEdge)*2] = vertexMap.at(vertex1) + 1;
					verticesOnEdgeNew[edgeMap.at(iEdge)*2+1] = vertexMap.at(vertex2) + 1;
				} else {
					cout << "ERROR: Edge mask is 1, but has no cells." << endl;
				}
			} else if(cell1 == -1){
				cellsOnEdgeNew[edgeMap.at(iEdge)*2] = cellMap.at(cell2) + 1;
				cellsOnEdgeNew[edgeMap.at(iEdge)*2+1] = 0;

				verticesOnEdgeNew[edgeMap.at(iEdge)*2] = vertexMap.at(vertex2) + 1;
				verticesOnEdgeNew[edgeMap.at(iEdge)*2+1] = vertexMap.at(vertex1) + 1;
			} else {
				cout << "ERROR: Edge mask is 1, but has no cells." << endl;
			}

#ifdef _DEBUG
			cout << "   New cell1: " << cellsOnEdgeNew[edgeMap.at(iEdge)*2] << endl;
			cout << "   New cell2: " << cellsOnEdgeNew[edgeMap.at(iEdge)*2+1] << endl;
			cout << "   New vertex1: " << verticesOnEdgeNew[edgeMap.at(iEdge)*2] << endl;
			cout << "   New vertex2: " << verticesOnEdgeNew[edgeMap.at(iEdge)*2+1] << endl;
#endif
		}
	}
	
	if (!(voeVar = grid.add_var("verticesOnEdge", ncInt, nEdgesDim, twoDim))) return NC_ERR;
	if (!voeVar->put(verticesOnEdgeNew,nEdgesNew,2)) return NC_ERR;
	if (!(coeVar = grid.add_var("cellsOnEdge", ncInt, nEdgesDim, twoDim))) return NC_ERR;
	if (!coeVar->put(cellsOnEdgeNew,nEdgesNew,2)) return NC_ERR;

	// Don't delete cellsOnEdgeOld yet. It's needed to map edgesOnEdge and weightsOnEdge
	delete[] cellsOnEdgeNew;
	delete[] verticesOnEdgeOld;
	delete[] verticesOnEdgeNew;

	// Map edgesOnEdge, nEdgesOnEdge, and weightsOnEdge
	nEdgesOnEdgeOld = new int[nEdges];
	nEdgesOnEdgeNew = new int[nEdges];
	edgesOnEdgeOld = new int[nEdges*maxEdges*2];
	edgesOnEdgeNew = new int[nEdgesNew*maxEdges2New];
	weightsOnEdgeOld = new double[nEdges*maxEdges*2];
	weightsOnEdgeNew = new double[nEdgesNew*maxEdges2New];

	netcdf_mpas_read_nedgesonedge( inputFilename, nEdges, nEdgesOnEdgeOld );
	netcdf_mpas_read_edgesonedge( inputFilename, nEdges, maxEdges*2, edgesOnEdgeOld );
	netcdf_mpas_read_weightsonedge( inputFilename, nEdges, maxEdges*2, weightsOnEdgeOld );

	for(int iEdge = 0; iEdge < nEdges; iEdge++){
		int edgeCount = 0;
		if(edgeMap.at(iEdge) != -1){
			int cell1, cell2;

			cell1 = cellsOnEdgeOld[iEdge * 2] - 1;
			cell2 = cellsOnEdgeOld[iEdge * 2 + 1] - 1;

			if(cell1 != -1){
				cell1 == cellMap.at(cell1);
			}

			if(cell2 != -1){
				cell2 == cellMap.at(cell1);
			}

			if(cell1 != -1 && cell2 != -1){
				for(int j = 0; j < nEdgesOnEdgeOld[iEdge]; j++){
					int eoe = edgesOnEdgeOld[iEdge*maxEdges*2 + j] - 1;

					if(eoe != -1 && eoe < edgeMap.size()){
						edgesOnEdgeNew[edgeMap.at(iEdge)*maxEdges2New + j] = edgeMap.at(eoe) + 1;
						weightsOnEdgeNew[edgeMap.at(iEdge)*maxEdges2New + j] = weightsOnEdgeOld[iEdge*maxEdges*2 + j];
						edgeCount++;
					} else {
						edgesOnEdgeNew[edgeMap.at(iEdge)*maxEdges2New + j] = 0;
						weightsOnEdgeNew[edgeMap.at(iEdge)*maxEdges2New + j] = 0;
					}
				}
			} else if ( cell1 == -1  || cell2 == -1){
				for(int j = 0; j < nEdgesOnEdgeOld[iEdge]; j++){
					edgesOnEdgeNew[edgeMap.at(iEdge)*maxEdges2New + j] = 0;
					weightsOnEdgeNew[edgeMap.at(iEdge)*maxEdges2New + j] = 0;
				}
			}

			for(int j = edgeCount; j < maxEdges2New; j++){
				edgesOnEdgeNew[edgeMap.at(iEdge)*maxEdges2New + j] = 0;
				weightsOnEdgeNew[edgeMap.at(iEdge)*maxEdges2New + j] = 0;
			}
			nEdgesOnEdgeNew[edgeMap.at(iEdge)] = edgeCount;
		}
	}

	if (!(nEoeVar = grid.add_var("nEdgesOnEdge", ncInt, nEdgesDim))) return NC_ERR;
	if (!nEoeVar->put(nEdgesOnEdgeNew,nEdgesNew)) return NC_ERR;
	if (!(eoeVar = grid.add_var("edgesOnEdge", ncInt, nEdgesDim, maxEdges2Dim))) return NC_ERR;
	if (!eoeVar->put(edgesOnEdgeNew,nEdgesNew,maxEdges2New)) return NC_ERR;
	if (!(woeVar = grid.add_var("weightsOnEdge", ncDouble, nEdgesDim, maxEdges2Dim))) return NC_ERR;
	if (!woeVar->put(weightsOnEdgeNew,nEdgesNew,maxEdges2New)) return NC_ERR;

	delete[] nEdgesOnEdgeOld;
	delete[] cellsOnEdgeOld;
	delete[] edgesOnEdgeOld;
	delete[] edgesOnEdgeNew;
	delete[] weightsOnEdgeOld;
	delete[] weightsOnEdgeNew;

	// Map dvEdge, dcEdge, and angleEdge
	dvEdgeOld = new double[nEdges];
	dcEdgeOld = new double[nEdges];
	angleEdgeOld = new double[nEdges];
	dvEdgeNew = new double[nEdgesNew];
	dcEdgeNew = new double[nEdgesNew];
	angleEdgeNew = new double[nEdgesNew];

	netcdf_mpas_read_dvedge ( inputFilename, nEdges, dvEdgeOld );
	netcdf_mpas_read_dcedge ( inputFilename, nEdges, dcEdgeOld );
	netcdf_mpas_read_angleedge ( inputFilename, nEdges, angleEdgeOld );

	for(int iEdge = 0; iEdge < nEdges; iEdge++){
		if(edgeMap.at(iEdge) != -1){
			dvEdgeNew[edgeMap.at(iEdge)] = dvEdgeOld[iEdge];
			dcEdgeNew[edgeMap.at(iEdge)] = dcEdgeOld[iEdge];
			angleEdgeNew[edgeMap.at(iEdge)] = angleEdgeOld[iEdge];
		}
	}

	if (!(dcEdgeVar = grid.add_var("dcEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!dcEdgeVar->put(dcEdgeNew,nEdgesNew)) return NC_ERR;
	if (!(dvEdgeVar = grid.add_var("dvEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!dvEdgeVar->put(dvEdgeNew,nEdgesNew)) return NC_ERR;
	if (!(angleVar = grid.add_var("angleEdge", ncDouble, nEdgesDim))) return NC_ERR;
	if (!angleVar->put(angleEdgeNew,nEdgesNew)) return NC_ERR;

	delete[] dvEdgeOld;
	delete[] dvEdgeNew;
	delete[] dcEdgeOld;
	delete[] dcEdgeNew;
	delete[] angleEdgeOld;
	delete[] angleEdgeNew;

	return 0;
}/*}}}*/
int mapAndOutputVertexFields( const string inputFilename, const string outputFilename) {/*{{{*/
	/*****************************************************************
	 *
	 * This function maps and writes all of the vertex related fields. Including
	 * cellsOnVertex
	 * edgesOnVertex
	 * areaTriangle
	 * kiteAreasOnVertex
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
	NcVar *covVar, *eovVar, *kaovVar, *atVar;

	int nVerticesNew = nVerticesDim->size();

	int *cellsOnVertexOld, *cellsOnVertexNew;
	int *edgesOnVertexOld, *edgesOnVertexNew;
	double *kiteAreasOnVertexOld, *kiteAreasOnVertexNew;
	double *areaTriangleOld, *areaTriangleNew;

	cellsOnVertexOld = new int[nVertices * vertexDegree];
	cellsOnVertexNew = new int[nVerticesNew * vertexDegree];
	edgesOnVertexOld = new int[nVertices * vertexDegree];
	edgesOnVertexNew = new int[nVerticesNew * vertexDegree];
	kiteAreasOnVertexOld = new double[nVertices * vertexDegree];
	kiteAreasOnVertexNew = new double[nVerticesNew * vertexDegree];
	areaTriangleOld = new double[nVertices];
	areaTriangleNew = new double[nVerticesNew];

	netcdf_mpas_read_cellsonvertex ( inputFilename, nVertices, vertexDegree, cellsOnVertexOld );
	netcdf_mpas_read_edgesonvertex ( inputFilename, nVertices, vertexDegree, edgesOnVertexOld );
	netcdf_mpas_read_kiteareasonvertex ( inputFilename, nVertices, vertexDegree, kiteAreasOnVertexOld );
	netcdf_mpas_read_areatriangle ( inputFilename, nVertices, areaTriangleOld );

	for(int iVertex = 0; iVertex < nVertices; iVertex++){
		double area = 0.0;
		if(vertexMap.at(iVertex) != -1){
			for(int j = 0; j < vertexDegree; j++){
				int iCell, iEdge;

				iCell = cellsOnVertexOld[iVertex*vertexDegree + j] - 1;
				iEdge = edgesOnVertexOld[iVertex*vertexDegree + j] - 1;

				if(iCell != -1){
					cellsOnVertexNew[ vertexMap.at(iVertex) * vertexDegree + j] = cellMap.at(iCell) + 1;
					if(cellMap.at(iCell) == -1){
						kiteAreasOnVertexNew[ vertexMap.at(iVertex) * vertexDegree + j] = 0.0;
					} else {
						kiteAreasOnVertexNew[ vertexMap.at(iVertex) * vertexDegree + j] = kiteAreasOnVertexOld[iVertex*vertexDegree + j];
					}
					area += kiteAreasOnVertexNew[ vertexMap.at(iVertex) * vertexDegree + j];
				} else {
					cellsOnVertexNew[ vertexMap.at(iVertex) * vertexDegree + j] = 0;
					kiteAreasOnVertexNew[ vertexMap.at(iVertex) * vertexDegree + j] = 0.0;
				}

				if(iEdge != -1){
					edgesOnVertexNew[ vertexMap.at(iVertex) * vertexDegree + j] = edgeMap.at(iEdge) + 1;
				} else {
					edgesOnVertexNew[ vertexMap.at(iVertex) * vertexDegree + j] = 0;
				}
			}

			areaTriangleNew[vertexMap.at(iVertex)] = area;
		}
	}

#ifdef _DEBUG
	cout << "   Writing edgesOnVertex" << endl;
#endif
	if (!(eovVar = grid.add_var("edgesOnVertex", ncInt, nVerticesDim, vertexDegreeDim))) return NC_ERR;
	if (!eovVar->put(edgesOnVertexNew,nVerticesNew,vertexDegree)) return NC_ERR;
#ifdef _DEBUG
	cout << "   Writing cellsOnVertex" << endl;
#endif
	if (!(covVar = grid.add_var("cellsOnVertex", ncInt, nVerticesDim, vertexDegreeDim))) return NC_ERR;
	if (!covVar->put(cellsOnVertexNew,nVerticesNew,vertexDegree)) return NC_ERR;
#ifdef _DEBUG
	cout << "   Writing areaTriangle" << endl;
#endif
	if (!(atVar = grid.add_var("areaTriangle", ncDouble, nVerticesDim))) return NC_ERR;
	if (!atVar->put(areaTriangleNew, nVerticesNew)) return NC_ERR;
#ifdef _DEBUG
	cout << "   Writing kiteAreasOnVertex" << endl;
#endif
	if (!(kaovVar = grid.add_var("kiteAreasOnVertex", ncDouble, nVerticesDim, vertexDegreeDim))) return NC_ERR;
	if (!kaovVar->put(kiteAreasOnVertexNew, nVerticesNew, vertexDegree)) return NC_ERR;

	delete[] cellsOnVertexOld;
	delete[] cellsOnVertexNew;
	delete[] edgesOnVertexOld;
	delete[] edgesOnVertexNew;
	delete[] kiteAreasOnVertexOld;
	delete[] kiteAreasOnVertexNew;
	delete[] areaTriangleOld;
	delete[] areaTriangleNew;

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

