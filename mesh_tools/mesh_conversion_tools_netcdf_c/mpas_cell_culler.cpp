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
#include "string_utils.h"

#define ID_LEN 10

using namespace std;

enum { mergeOp, invertOp, preserveOp };

size_t nCells, nVertices, nEdges, vertexDegree, maxEdges;
bool spherical, periodic;
bool cullMasks = false;
double sphere_radius, xPeriod, yPeriod;
string in_history = "";
string in_file_id = "";
string in_parent_id = "";
string in_mesh_spec = "1.0";
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
int mapAndOutputCellFields(const string inputFilename, const string outputPath,
                           const string outputFilename);
int mapAndOutputEdgeFields(const string inputFilename, const string outputFilename);
int mapAndOutputVertexFields(const string inputFilename, const string outputFilename);
int outputCellMap(const string outputPath);
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
    string out_path = "";
    string out_file = "";
    string out_fext = "";
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
    else
    {
        cullMasks = true;
        in_name = argv[1];
        out_name = argv[2];
        bool foundOperation;

        int i = 3;
        while ( i < argc ) {
            foundOperation = false;
            if (strcmp(argv[i], "-m") == 0 ) {
                mask_ops.push_back(static_cast<int>(mergeOp));
                foundOperation = true;
            } else if ( strcmp(argv[i], "-i") == 0 ){
                mask_ops.push_back(static_cast<int>(invertOp));
                foundOperation = true;
            } else if ( strcmp(argv[i], "-p") == 0 ){
                mask_ops.push_back(static_cast<int>(preserveOp));
                foundOperation = true;
            } else if ( strcmp(argv[i], "-c") == 0 ){
                outputMap = true;
            } else {
                cout << " ERROR: Invalid option passed on the command line " << argv[i] << ". Exiting..." << endl;
                print_usage();
                exit(1);
            }
            i++;

            if (foundOperation) {
                mask_names.push_back( argv[i] );
                i++;
            }
        }
    }

    if(out_name == in_name){
        cout << "   ERROR: Input and Output names are the same." << endl;
        return 1;
    }

    file_part(out_name, out_path, out_file, out_fext);

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
    if(error = mapAndOutputCellFields(in_name, out_path, out_name)){
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
        if(error = outputCellMap(out_path)){
            cout << "Error - " << error << endl;
            exit(error);
        }
    }

    return 0;
}

int outputCellMap(const string outputPath){/*{{{*/

    int iCell;
    ofstream outputfileForward, outputfileBackward;

    // forwards mapping
    outputfileForward.open(path_join(outputPath, "cellMapForward.txt"));

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

    outputfileBackward.open(path_join(outputPath, "cellMapBackward.txt"));

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
    string on_a_sphere, is_periodic;

#ifdef _DEBUG
    cout << endl << endl << "Begin function: readGridInput" << endl << endl;
#endif

    ncutil::get_dim(inputFilename, "nCells", nCells);
    ncutil::get_dim(inputFilename, "nEdges", nEdges);
    ncutil::get_dim(inputFilename, "nVertices", nVertices);
    ncutil::get_dim(inputFilename, "vertexDegree", vertexDegree);
    ncutil::get_dim(inputFilename, "maxEdges", maxEdges);

    try {
#ifdef _DEBUG
    cout << "   Reading on_a_sphere" << endl;
#endif
    ncutil::get_str(inputFilename, "on_a_sphere", on_a_sphere);
    spherical = (on_a_sphere.find("YES") != string::npos);
    } catch (...) {
    // allow errors for optional attr. not found
    }
    try {
#ifdef _DEBUG
    cout << "   Reading sphere_radius" << endl;
#endif
    ncutil::get_att(inputFilename, "sphere_radius", &sphere_radius);
    } catch (...) {
    // allow errors for optional attr. not found
    }
    try {
#ifdef _DEBUG
    cout << "   Reading history" << endl;
#endif
    ncutil::get_str(inputFilename, "history", in_history);
    } catch (...) {
    // allow errors for optional attr. not found
    }
    try {
#ifdef _DEBUG
    cout << "   Reading file_id" << endl;
#endif
    ncutil::get_str(inputFilename, "file_id", in_file_id);
    } catch (...) {
    // allow errors for optional attr. not found
    }
    try {
#ifdef _DEBUG
    cout << "   Reading parent_id" << endl;
#endif
    ncutil::get_str(inputFilename, "parent_id", in_parent_id);
    } catch (...) {
    // allow errors for optional attr. not found
    }
    try {
#ifdef _DEBUG
    cout << "   Reading parent_id" << endl;
#endif
    ncutil::get_str(inputFilename, "mesh_spec", in_mesh_spec);
    } catch (...) {
    // allow errors for optional attr. not found
    }
    try {
#ifdef _DEBUG
    cout << "   Reading is_periodic" << endl;
#endif
    ncutil::get_str(inputFilename, "is_periodic", is_periodic);
    periodic = (is_periodic.find("YES") != string::npos);
    } catch (...) {
    // allow errors for optional attr. not found
    }
    try {
#ifdef _DEBUG
    cout << "   Reading x_period" << endl;
#endif
    ncutil::get_att(inputFilename, "x_period", &xPeriod);
    } catch (...) {
    // allow errors for optional attr. not found
    }
    try {
#ifdef _DEBUG
    cout << "   Reading y_period" << endl;
#endif
    ncutil::get_att(inputFilename, "y_period", &yPeriod);

    } catch (...) {
    // allow errors for optional attr. not found
    }

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
    ncutil::get_var(inputFilename, "areaCell", &areaCell[0]);

#ifdef _DEBUG
    cout << " Read cullCell" << endl;
#endif
    cullCell.clear();
    cullCell.resize(nCells);
    try {
        ncutil::get_var(inputFilename, "cullCell", &cullCell[0]);
    } catch (...) {
    // allow errors for optional vars. not found
    }

    // Build cellsOnVertex information
    cellsonvertex_list = new int[nVertices * vertexDegree];

#ifdef _DEBUG
    cout << " Read cellsOnVertex" << endl;
#endif
    ncutil::get_var(inputFilename, "cellsOnVertex", cellsonvertex_list);
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
    ncutil::get_var(inputFilename, "verticesOnEdge", verticesonedge_list);
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
    ncutil::get_var(inputFilename, "cellsOnEdge", cellsonedge_list);
    for(int i = 0; i < nEdges; i++){
        // Subtract 1 to convert into base 0 (c index space).
        cellsOnEdge.at(i).push_back(cellsonedge_list[i*2] - 1);
        cellsOnEdge.at(i).push_back(cellsonedge_list[i*2+1] - 1);
    }

    delete[] cellsonedge_list;

    return 0;
}/*}}}*/
int mergeCellMasks(const string masksFilename, const int maskOp){/*{{{*/
    size_t nRegions = 0, nTransects = 0;
    int *regionCellMasks, *transectCellMasks, *cellSeedMask, *flattenedMask;
    int i, j;

    try {
        ncutil::get_dim(masksFilename, "nRegions", nRegions);
    } catch (...) {
    // allow errors for optional attr. not found
    }
    try {
        ncutil::get_dim(masksFilename, "nTransects", nTransects);
    } catch (...) {
    // allow errors for optional attr. not found
    }

    cout << "mask filename: " << masksFilename << endl;
    cout << "  nRegions: " << nRegions << endl;
    cout << "  nTransects: " << nTransects << endl;
    cout << "  maskOp: " << maskOp << endl;
    cout << "  nCells: " << nCells << endl;

    regionCellMasks = new int[nCells*nRegions];
    transectCellMasks = new int[nCells*nTransects];
    cellSeedMask = new int[nCells];
    flattenedMask = new int[nCells];

    for ( i = 0; i < nCells; i++){
        cellSeedMask[i] = 0;
    }
    try {
        ncutil::get_var(masksFilename, "cellSeedMask", cellSeedMask);
    } catch (...) {
    // allow errors for optional vars. not found
    }

    if (nRegions > 0) {
        cout << "  Reading regionCellMasks" << endl;
        ncutil::get_var(masksFilename, "regionCellMasks", regionCellMasks);
    }
    if (nTransects > 0) {
        cout << "  Reading transectCellMasks" << endl;
        ncutil::get_var(masksFilename, "transectCellMasks", transectCellMasks);
    }

    cout << "  Flattening seed, region and/or transect masks" << endl;
    for ( i = 0; i < nCells; i++){
        flattenedMask[i] = cellSeedMask[i];
        for ( j = 0; j < nRegions; j++){
            flattenedMask[i] = max(flattenedMask[i], regionCellMasks[i * nRegions + j]);
        }

        for ( j = 0; j < nTransects; j++ ) {
            flattenedMask[i] = max(flattenedMask[i], transectCellMasks[i * nTransects + j]);
        }
    }

    cout << "  Applying flattened mask to cullCell" << endl;

    if ( maskOp == invertOp || maskOp == mergeOp ) {
        if ( maskOp == invertOp ) {
            cout << "  Inverting flattened mask" << endl;
            for (i = 0; i < nCells; i++){
                flattenedMask[i] = (flattenedMask[i] + 1) % 2;
            }
        }

        cout << "  Masking cullCell with the flattened mask" << endl;
        for ( i = 0; i < nCells; i++ ){
            cullCell[i] = max(cullCell[i], flattenedMask[i]);
        }
    } else if ( maskOp == preserveOp ) {
        cout << "  Preserving cells in cullCell with the flattened mask" << endl;
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

    int grid, retv;
    if ((retv = nc_create(outputFilename.c_str(), NC_CLOBBER|NC_NETCDF4, &grid)))
    {
        std::cout << "Can't create file: " << outputFilename << std::endl;
        return retv ;
    }

    size_t nCellsNew, nEdgesNew, nVerticesNew;

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

    ncutil::def_dim(outputFilename, "nCells", nCellsNew);
    ncutil::def_dim(outputFilename, "nEdges", nEdgesNew);
    ncutil::def_dim(outputFilename, "nVertices", nVerticesNew);
    ncutil::def_dim(outputFilename, "TWO", 2);
    ncutil::def_dim(outputFilename, "vertexDegree", vertexDegree);
    ncutil::def_dim(outputFilename, "Time", NC_UNLIMITED);

    return 0;
}/*}}}*/
int outputGridAttributes( const string inputFilename, const string outputFilename ){/*{{{*/
    /************************************************************************
     *
     * This function writes the grid dimensions to the netcdf file named
     * outputFilename
     *
     * **********************************************************************/

    string history_str = "";
    string id_str = "";
    string parent_str = "";

    // write attributes
    if(!spherical){
        ncutil::put_str(outputFilename, "on_a_sphere", "NO");
        ncutil::put_att(outputFilename, "sphere_radius", NC_DOUBLE, 0.);
    } else {
        ncutil::put_str(outputFilename, "on_a_sphere", "YES");
        ncutil::put_att(outputFilename,
            "sphere_radius", NC_DOUBLE, sphere_radius);
    }

    if(!periodic){
        ncutil::put_str(outputFilename, "is_periodic", "NO");
    } else {
        ncutil::put_str(outputFilename, "is_periodic", "YES");
        ncutil::put_att(outputFilename, "x_period", NC_DOUBLE, xPeriod);
        ncutil::put_att(outputFilename, "y_period", NC_DOUBLE, yPeriod);
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
        ncutil::put_str(outputFilename, "parent_id", parent_str);
    }

    id_str = gen_random(ID_LEN);

    ncutil::put_str(outputFilename, "history", history_str);
    ncutil::put_str(outputFilename, "mesh_spec", in_mesh_spec);
    ncutil::put_str(outputFilename, "Conventions", "MPAS");
    ncutil::put_str(outputFilename, "source", "MpasCellCuller.x");
    ncutil::put_str(outputFilename, "file_id", id_str);

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

    size_t nCellsNew, nEdgesNew, nVerticesNew;
    ncutil::get_dim(outputFilename, "nCells", nCellsNew);
    ncutil::get_dim(outputFilename, "nEdges", nEdgesNew);
    ncutil::get_dim(outputFilename, "nVertices", nVerticesNew);

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

    ncutil::get_var(inputFilename, "xCell", xOld);
    ncutil::get_var(inputFilename, "yCell", yOld);
    ncutil::get_var(inputFilename, "zCell", zOld);
    ncutil::get_var(inputFilename, "latCell", latOld);
    ncutil::get_var(inputFilename, "lonCell", lonOld);

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

    ncutil::def_var(outputFilename, "latCell",
        NC_DOUBLE, "latitudes of cell centres", {"nCells"});
    ncutil::def_var(outputFilename, "lonCell",
        NC_DOUBLE, "longitudes of cell centres", {"nCells"});

    ncutil::put_var(outputFilename, "latCell", &latNew[0]);
    ncutil::put_var(outputFilename, "lonCell", &lonNew[0]);

    ncutil::def_var(outputFilename, "xCell",
        NC_DOUBLE, "x-coordinates of cell centres", {"nCells"});
    ncutil::def_var(outputFilename, "yCell",
        NC_DOUBLE, "y-coordinates of cell centres", {"nCells"});
    ncutil::def_var(outputFilename, "zCell",
        NC_DOUBLE, "z-coordinates of cell centres", {"nCells"});

    ncutil::put_var(outputFilename, "xCell", &xNew[0]);
    ncutil::put_var(outputFilename, "yCell", &yNew[0]);
    ncutil::put_var(outputFilename, "zCell", &zNew[0]);

    ncutil::def_var(outputFilename, "indexToCellID",
        NC_INT, "index to cell ID mapping", {"nCells"});

    ncutil::put_var(outputFilename, "indexToCellID", &idxToNew[0]);

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

    ncutil::get_var(inputFilename, "xEdge", xOld);
    ncutil::get_var(inputFilename, "yEdge", yOld);
    ncutil::get_var(inputFilename, "zEdge", zOld);
    ncutil::get_var(inputFilename, "latEdge", latOld);
    ncutil::get_var(inputFilename, "lonEdge", lonOld);

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

    ncutil::def_var(outputFilename, "latEdge",
        NC_DOUBLE, "latitudes of edge centres", {"nEdges"});
    ncutil::def_var(outputFilename, "lonEdge",
        NC_DOUBLE, "longitudes of edge centres", {"nEdges"});

    ncutil::put_var(outputFilename, "latEdge", &latNew[0]);
    ncutil::put_var(outputFilename, "lonEdge", &lonNew[0]);

    ncutil::def_var(outputFilename, "xEdge",
        NC_DOUBLE, "x-coordinates of edge centres", {"nEdges"});
    ncutil::def_var(outputFilename, "yEdge",
        NC_DOUBLE, "y-coordinates of edge centres", {"nEdges"});
    ncutil::def_var(outputFilename, "zEdge",
        NC_DOUBLE, "z-coordinates of edge centres", {"nEdges"});

    ncutil::put_var(outputFilename, "xEdge", &xNew[0]);
    ncutil::put_var(outputFilename, "yEdge", &yNew[0]);
    ncutil::put_var(outputFilename, "zEdge", &zNew[0]);

    ncutil::def_var(outputFilename, "indexToEdgeID",
        NC_INT, "index to edge ID mapping", {"nEdges"});

    ncutil::put_var(outputFilename, "indexToEdgeID", &idxToNew[0]);

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

    ncutil::get_var(inputFilename, "xVertex", xOld);
    ncutil::get_var(inputFilename, "yVertex", yOld);
    ncutil::get_var(inputFilename, "zVertex", zOld);
    ncutil::get_var(inputFilename, "latVertex", latOld);
    ncutil::get_var(inputFilename, "lonVertex", lonOld);

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

    ncutil::def_var(outputFilename, "latVertex",
        NC_DOUBLE, "latitudes of vertices", {"nVertices"});
    ncutil::def_var(outputFilename, "lonVertex",
        NC_DOUBLE, "longitudes of vertices", {"nVertices"});

    ncutil::put_var(outputFilename, "latVertex", &latNew[0]);
    ncutil::put_var(outputFilename, "lonVertex", &lonNew[0]);

    ncutil::def_var(outputFilename, "xVertex",
        NC_DOUBLE, "x-coordinates of vertices", {"nVertices"});
    ncutil::def_var(outputFilename, "yVertex",
        NC_DOUBLE, "y-coordinates of vertices", {"nVertices"});
    ncutil::def_var(outputFilename, "zVertex",
        NC_DOUBLE, "z-coordinates of vertices", {"nVertices"});

    ncutil::put_var(outputFilename, "xVertex", &xNew[0]);
    ncutil::put_var(outputFilename, "yVertex", &yNew[0]);
    ncutil::put_var(outputFilename, "zVertex", &zNew[0]);

    ncutil::def_var(outputFilename, "indexToVertexID",
        NC_INT, "index to vertex ID mapping", {"nVertices"});

    ncutil::put_var(outputFilename, "indexToVertexID", &idxToNew[0]);

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

    return 0;
}/*}}}*/
int mapAndOutputCellFields( const string inputFilename, const string outputPath,
                            const string outputFilename) {/*{{{*/
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

    size_t nCellsNew, nEdgesNew, maxEdgesNew, edgeCount;
    ncutil::get_dim(outputFilename, "nCells", nCellsNew);
    ncutil::get_dim(outputFilename, "nEdges", nEdgesNew);

    double *meshDensityOld, *meshDensityNew;
    double *areaCellNew;
    int *tmp_arr_old, *nEdgesOnCellOld, *nEdgesOnCellNew;
    int *tmp_arr_new;

    tmp_arr_old = new int[nCells*maxEdges];
    nEdgesOnCellOld = new int[nCells];
    nEdgesOnCellNew = new int[nCellsNew];

    ncutil::get_var(inputFilename, "edgesOnCell", tmp_arr_old);
    ncutil::get_var(inputFilename, "nEdgesOnCell", nEdgesOnCellOld);

    // Need to map nEdgesOnCell to get maxEdges
    maxEdgesNew = 0;
    for(int iCell = 0; iCell < nCells; iCell++){
        if(cellMap.at(iCell) != -1){
            nEdgesOnCellNew[cellMap.at(iCell)] = nEdgesOnCellOld[iCell];
            maxEdgesNew = max(maxEdgesNew, (size_t)nEdgesOnCellNew[cellMap.at(iCell)]);
        }
    }
    tmp_arr_new = new int[nCells * maxEdgesNew];

    // Write maxEdges and maxEdges2 to output file
    ncutil::def_dim(outputFilename, "maxEdges", maxEdgesNew);
    ncutil::def_dim(outputFilename, "maxEdges2", maxEdgesNew * 2);

    // Write nEdgesOncell to output file
    ncutil::def_var(outputFilename, "nEdgesOnCell",
        NC_INT, "number of edges on each cell", {"nCells"});

    ncutil::put_var(outputFilename, "nEdgesOnCell", &nEdgesOnCellNew[0]);

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
                cout << "    Mapping edge: " << iEdge << " to " <<
                    tmp_arr_new[cellMap.at(iCell)*maxEdgesNew + j] << " dbg info: " <<
                        tmp_arr_old[iCell*maxEdges + j] << " " << j << endl;
#endif
            }
        }
    }

    ncutil::def_var(outputFilename, "edgesOnCell",
        NC_INT, "edges on each cell", {"nCells", "maxEdges"});

    ncutil::put_var(outputFilename, "edgesOnCell", &tmp_arr_new[0]);

    ncutil::get_var(inputFilename, "cellsOnCell", tmp_arr_old);

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
    ofstream graph(path_join(outputPath, "culled_graph.info"));
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

    ncutil::def_var(outputFilename, "cellsOnCell",
        NC_INT, "cells adj. to each cell", {"nCells", "maxEdges"});

    ncutil::put_var(outputFilename, "cellsOnCell", &tmp_arr_new[0]);

    delete[] nEdgesOnCellNew;
    delete[] nEdgesOnCellOld;

    ncutil::get_var(inputFilename, "verticesOnCell", tmp_arr_old);

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

    ncutil::def_var(outputFilename, "verticesOnCell",
        NC_INT, "vertices on each cell", {"nCells", "maxEdges"});

    ncutil::put_var(outputFilename, "verticesOnCell", &tmp_arr_new[0]);

    delete[] tmp_arr_old;
    delete[] tmp_arr_new;

    // Map areaCell
    areaCellNew = new double[nCellsNew];

    for(int iCell = 0; iCell < nCells; iCell++){
        if(cellMap.at(iCell) != -1){
            areaCellNew[cellMap.at(iCell)] = areaCell.at(iCell);
        }
    }

    ncutil::def_var(outputFilename, "areaCell",
        NC_DOUBLE, "surface area of each cell", {"nCells"});

    ncutil::put_var(outputFilename, "areaCell", &areaCellNew[0]);

    delete[] areaCellNew;

    // Map meshDensity
    meshDensityOld = new double[nCells];
    meshDensityNew = new double[nCellsNew];

    ncutil::get_var(inputFilename, "meshDensity", meshDensityOld);

    for(int iCell = 0; iCell < nCells; iCell++){
        if(cellMap.at(iCell) != -1){
            meshDensityNew[cellMap.at(iCell)] = meshDensityOld[iCell];
        }
    }

    ncutil::def_var(outputFilename, "meshDensity",
        NC_DOUBLE, "mesh density distribution", {"nCells"});

    ncutil::put_var(outputFilename, "meshDensity", &meshDensityNew[0]);

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

    size_t nEdgesNew, maxEdges2New, two = 2;
    ncutil::get_dim(outputFilename, "nEdges", nEdgesNew);
    ncutil::get_dim(outputFilename, "maxEdges2", maxEdges2New);

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

    ncutil::get_var(inputFilename, "cellsOnEdge", cellsOnEdgeOld);
    ncutil::get_var(inputFilename, "verticesOnEdge", verticesOnEdgeOld);

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

    ncutil::def_var(outputFilename, "verticesOnEdge",
        NC_INT, "vertices on each edge", {"nEdges", "TWO"});
    ncutil::def_var(outputFilename, "cellsOnEdge",
        NC_INT, "cells adj. to each edge", {"nEdges", "TWO"});

    ncutil::put_var(outputFilename, "verticesOnEdge", &verticesOnEdgeNew[0]);
    ncutil::put_var(outputFilename, "cellsOnEdge", &cellsOnEdgeNew[0]);

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

    ncutil::get_var(inputFilename, "nEdgesOnEdge", nEdgesOnEdgeOld);
    ncutil::get_var(inputFilename, "edgesOnEdge", edgesOnEdgeOld);
    ncutil::get_var(inputFilename, "weightsOnEdge", weightsOnEdgeOld);

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
                        weightsOnEdgeNew[edgeMap.at(iEdge)*maxEdges2New + j] =
                            weightsOnEdgeOld[iEdge*maxEdges*2 + j];
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

    ncutil::def_var(outputFilename, "nEdgesOnEdge",
        NC_INT, "number of edges adj. to each edge", {"nEdges"});
    ncutil::def_var(outputFilename, "edgesOnEdge",
        NC_INT, "edges adj. to each edge", {"nEdges", "maxEdges2"});

    ncutil::put_var(outputFilename, "nEdgesOnEdge", &nEdgesOnEdgeNew[0]);
    ncutil::put_var(outputFilename, "edgesOnEdge", &edgesOnEdgeNew[0]);

    ncutil::def_var(outputFilename, "weightsOnEdge",
        NC_DOUBLE, "tangential flux reconstruction weights", {"nEdges", "maxEdges2"});

    ncutil::put_var(outputFilename, "weightsOnEdge", &weightsOnEdgeNew[0]);

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

    ncutil::get_var(inputFilename, "dvEdge", dvEdgeOld);
    ncutil::get_var(inputFilename, "dcEdge", dcEdgeOld);
    ncutil::get_var(inputFilename, "angleEdge", angleEdgeOld);

    for(int iEdge = 0; iEdge < nEdges; iEdge++){
        if(edgeMap.at(iEdge) != -1){
            dvEdgeNew[edgeMap.at(iEdge)] = dvEdgeOld[iEdge];
            dcEdgeNew[edgeMap.at(iEdge)] = dcEdgeOld[iEdge];
            angleEdgeNew[edgeMap.at(iEdge)] = angleEdgeOld[iEdge];
        }
    }

    ncutil::def_var(outputFilename, "dvEdge",
        NC_DOUBLE, "length of arc between centres", {"nEdges"});
    ncutil::def_var(outputFilename, "dcEdge",
        NC_DOUBLE, "length of arc between centres", {"nEdges"});
    ncutil::def_var(outputFilename, "angleEdge",
        NC_DOUBLE, "angle to edges", {"nEdges"}) ;

    ncutil::put_var(outputFilename, "dvEdge", &dvEdgeNew[0]);
    ncutil::put_var(outputFilename, "dcEdge", &dcEdgeNew[0]);
    ncutil::put_var(outputFilename, "angleEdge", &angleEdgeNew[0]);

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

    size_t nVerticesNew;
    ncutil::get_dim(outputFilename, "nVertices", nVerticesNew);

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

    ncutil::get_var(inputFilename, "cellsOnVertex", cellsOnVertexOld);
    ncutil::get_var(inputFilename, "edgesOnVertex", edgesOnVertexOld);
    ncutil::get_var(inputFilename, "kiteAreasOnVertex", kiteAreasOnVertexOld);
    ncutil::get_var(inputFilename, "areaTriangle", areaTriangleOld);

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
                        kiteAreasOnVertexNew[ vertexMap.at(iVertex) * vertexDegree + j] =
                            kiteAreasOnVertexOld[iVertex*vertexDegree + j];
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

    ncutil::def_var(outputFilename, "edgesOnVertex",
        NC_INT, "edges adj. to each vertex", {"nVertices", "vertexDegree"});
    ncutil::def_var(outputFilename, "cellsOnVertex",
        NC_INT, "cells adj. to each vertex", {"nVertices", "vertexDegree"});

    ncutil::put_var(outputFilename, "edgesOnVertex", &edgesOnVertexNew [0]);
    ncutil::put_var(outputFilename, "cellsOnVertex", &cellsOnVertexNew [0]);

    ncutil::def_var(outputFilename, "areaTriangle",
        NC_DOUBLE, "surface area of dual cells", {"nVertices"});
    ncutil::def_var(outputFilename, "kiteAreasOnVertex",
        NC_DOUBLE,
    "surface areas of overlap between cells and dual cells", {"nVertices", "vertexDegree"});

    ncutil::put_var(outputFilename, "areaTriangle", &areaTriangleNew [0]);
    ncutil::put_var(outputFilename, "kiteAreasOnVertex", &kiteAreasOnVertexNew [0]);

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
//        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";

    string rand_str = "";

    for (int i = 0; i < len; ++i) {
        rand_str += alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    return rand_str;
}/*}}}*/

