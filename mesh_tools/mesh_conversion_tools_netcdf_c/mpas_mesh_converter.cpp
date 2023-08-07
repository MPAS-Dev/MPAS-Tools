#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <math.h>
#include <assert.h>
#include <unordered_set>
#include <time.h>
#include <float.h>

#include "netcdf_utils.h"

#include "pnt.h"
#include "edge.h"

#include "string_utils.h"

#define MESH_SPEC 1.0
#define ID_LEN 10

using namespace std;

int nCells, nVertices, vertex_degree;
int maxEdges;
int obtuseTriangles = 0;
bool spherical, periodic;
double sphereRadius, xPeriod, yPeriod;
double xCellRange[2];
double yCellRange[2];
double zCellRange[2];
double xVertexRange[2];
double yVertexRange[2];
double zVertexRange[2];
double xCellDistance, yCellDistance, zCellDistance;
double xVertexDistance, yVertexDistance, zVertexDistance;
double xPeriodicFix, yPeriodicFix;
string in_history = "";
string in_file_id = "";
string in_parent_id = "";

// Connectivity and location information {{{

vector<pnt> cells;
vector<pnt> edges;
vector<pnt> vertices;
vector<int> completeCellMask;
vector<int> nEdgesOnCell;
vector< vector<int> > cellsOnEdge;
vector< vector<int> > verticesOnEdge;
vector< vector<int> > edgesOnVertex;
vector< vector<int> > cellsOnVertex;
vector< vector<int> > cellsOnCell;
vector< vector<int> > edgesOnCell;
vector< vector<int> > verticesOnCell;
vector< vector<int> > edgesOnEdge;
vector< vector<double> > weightsOnEdge;
vector< vector<double> > kiteAreasOnVertex;
vector<double> dvEdge;
vector<double> dcEdge;
vector<double> areaCell;
vector<double> areaTriangle;
vector<double> angleEdge;
vector<double> meshDensity;
vector<double> cellQuality;
vector<double> gridSpacing;
vector<double> triangleQuality;
vector<double> triangleAngleQuality;
vector<int> obtuseTriangle;

// }}}

// Iterators {{{
vector<pnt>::iterator pnt_itr;
vector<int>::iterator int_itr;
vector< vector<int> >::iterator vec_int_itr;
vector< vector<double> >::iterator vec_dbl_itr;
vector<double>::iterator dbl_itr;
// }}}
struct int_hasher {/*{{{*/
    size_t operator()(const int i) const {
        return (size_t)i;
    }
};/*}}}*/

/* Building/Ordering functions {{{ */
int readGridInput(const string inputFilename);
int buildUnorderedCellConnectivity();
int firstOrderingVerticesOnCell();
int buildCompleteCellMask();
int buildEdges();
int orderVertexArrays();
int orderCellArrays();
int buildAreas();
int buildEdgesOnEdgeArrays();
int buildAngleEdge();
int buildMeshQualities();
/*}}}*/

/* Output functions {{{*/
int outputGridDimensions(const string outputFilename);
int outputGridAttributes(const string outputFilename, const string inputFilename);
int outputGridCoordinates(const string outputFilename);
int outputCellConnectivity(const string outputFilename);
int outputEdgeConnectivity(const string outputFilename);
int outputVertexConnectivity(const string outputFilename);
int outputCellParameters(const string outputFilename);
int outputVertexParameters(const string outputFilename);
int outputEdgeParameters(const string outputFilename);
int outputMeshDensity(const string outputFilename);
int outputMeshQualities(const string outputFilename);
int writeGraphFile(const string outputFilename);
/*}}}*/

string gen_random(const int len);

int main ( int argc, char *argv[] ) {
    int error;
    string out_name = "mesh.nc";
    string in_name = "grid.nc";
    string out_path = "";
    string out_file = "";
    string out_fext = "";

    cout << endl << endl;
    cout << "************************************************************" << endl;
    cout << "MPAS_MESH_CONVERTER:\n";
    cout << "  C++ version\n";
    cout << "  Convert a NetCDF file describing Cell Locations, \n";
    cout << "  Vertex Location, and Connectivity into a valid MPAS mesh.\n";
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
        cout << "MPAS_MESH_CONVERTER:\n";
        cout << "  Please enter the NetCDF input filename.\n";

        cin >> in_name;

        cout << "\n";
        cout << "MPAS_MESH_CONVERTER:\n";
        cout << "  Please enter the output NetCDF MPAS Mesh filename.\n";

        cin >> out_name;
    }
    else if (argc == 2)
    {
        in_name = argv[1];

        cout << "\n";
        cout << "MPAS_MESH_CONVERTER:\n";
        cout << "  Output name not specified. Using default of mesh.nc\n";
    }
    else if (argc == 3)
    {
        in_name = argv[1];
        out_name = argv[2];
    }

    if(in_name == out_name){
        cout << "   ERROR: Input and output names are the same." << endl;
        return 1;
    }

    file_part(out_name, out_path, out_file, out_fext);

    srand(time(NULL));

    cout << "Reading input grid." << endl;
    error = readGridInput(in_name);
    if(error) return 1;

    cout << "Build prelimiary cell connectivity." << endl;
    error = buildUnorderedCellConnectivity();
    if(error) return 1;

    cout << "Order vertices on cell." << endl;
    error = firstOrderingVerticesOnCell();
    if(error) return 1;

    cout << "Build complete cell mask." << endl;
    error = buildCompleteCellMask();
    if(error) return 1;

    cout << "Build and order edges, dvEdge, and dcEdge." << endl;
    error = buildEdges();
    if(error) return 1;

    cout << "Build and order vertex arrays," << endl;
    error = orderVertexArrays();
    if(error) return 1;

    cout << "Build and order cell arrays," << endl;
    error = orderCellArrays();
    if(error) return 1;

    cout << "Build areaCell, areaTriangle, and kiteAreasOnVertex." << endl;
    error = buildAreas();
    if(error) return 1;

    cout << "Build edgesOnEdge and weightsOnEdge." << endl;
    error = buildEdgesOnEdgeArrays();
    if(error) return 1;

    cout << "Build angleEdge." << endl;
    error = buildAngleEdge();
    if(error) return 1;

    cout << "Building mesh qualities." << endl;
    error = buildMeshQualities();
    if(error) return 1;

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
    cout << "Writing cell connectivity" << endl;
    if(error = outputCellConnectivity(out_name)){
        cout << "Error - " << error << endl;
        exit(error);
    }
    cout << "Writing edge connectivity" << endl;
    if(error = outputEdgeConnectivity(out_name)){
        cout << "Error - " << error << endl;
        exit(error);
    }
    cout << "Writing vertex connectivity" << endl;
    if(error = outputVertexConnectivity(out_name)){
        cout << "Error - " << error << endl;
        exit(error);
    }
    cout << "Writing cell parameters" << endl;
    if(error = outputCellParameters(out_name)){
        cout << "Error - " << error << endl;
        exit(error);
    }
    cout << "Writing edge parameters" << endl;
    if(error = outputEdgeParameters(out_name)){
        cout << "Error - " << error << endl;
        exit(error);
    }
    cout << "Writing vertex parameters" << endl;
    if(error = outputVertexParameters(out_name)){
        cout << "Error - " << error << endl;
        exit(error);
    }

    cout << "Writing mesh qualities" << endl;
    if(error = outputMeshQualities(out_name)){
        cout << "Error - " << error << endl;
        exit(error);
    }

    cout << "Reading and writing meshDensity" << endl;
    if(error = outputMeshDensity(out_name)){
        cout << "Error - " << error << endl;
        exit(error);
    }

    cout << "Write graph.info file" << endl;
    if(error = writeGraphFile(path_join(out_path, "graph.info"))){
        cout << "Error - " << error << endl;
        exit(error);
    }
}

/* Building/Ordering functions {{{ */
int readGridInput(const string inputFilename){/*{{{*/
    int ncid, retv;
    size_t nCells, nVertices, vertexDegree;
    string on_a_sphere, is_periodic;
    double *xcell, *ycell,*zcell;
    double *xvertex, *yvertex,*zvertex;
    int *cellsonvertex_list;
    pnt new_location;

#ifdef _DEBUG
    cout << endl << endl << "Begin function: readGridInput" << endl << endl;
#endif

    // Initialize range mins with huge values.
    xCellRange[0] = 1E10;
    yCellRange[0] = 1E10;
    zCellRange[0] = 1E10;
    xVertexRange[0] = 1E10;
    yVertexRange[0] = 1E10;
    zVertexRange[0] = 1E10;

    // Initialize range maxes with small values.
    xCellRange[1] = -1E10;
    yCellRange[1] = -1E10;
    zCellRange[1] = -1E10;
    xVertexRange[1] = -1E10;
    yVertexRange[1] = -1E10;
    zVertexRange[1] = -1E10;

    ncutil::get_dim(inputFilename, "nCells", nCells);
    ncutil::get_dim(inputFilename, "nVertices", nVertices);
    ncutil::get_dim(inputFilename, "vertexDegree", vertexDegree);

    vertex_degree = vertexDegree;

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
    ncutil::get_att(inputFilename, "sphere_radius", &sphereRadius);
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
    cout << "    vertexDegree = " << vertexDegree << endl;
    cout << "    Spherical? = " << spherical << endl;
    cout << "    Periodic? = " << periodic << endl;
    if ( periodic ) {
        cout << "    x_period = " << xPeriod << endl;
        cout << "    y_period = " << yPeriod << endl;
    }

    // Build cell center location information
    xcell = new double[nCells];
    ycell = new double[nCells];
    zcell = new double[nCells];

    ncutil::get_var(inputFilename, "xCell", xcell);
    ncutil::get_var(inputFilename, "yCell", ycell);
    ncutil::get_var(inputFilename, "zCell", zcell);

    cells.clear();
    for(int i = 0; i < nCells; i++){
        xCellRange[0] = min(xCellRange[0], xcell[i]);
        xCellRange[1] = max(xCellRange[1], xcell[i]);

        yCellRange[0] = min(yCellRange[0], ycell[i]);
        yCellRange[1] = max(yCellRange[1], ycell[i]);

        zCellRange[0] = min(zCellRange[0], zcell[i]);
        zCellRange[1] = max(zCellRange[1], zcell[i]);

        new_location = pnt(xcell[i], ycell[i], zcell[i], i);

        if(spherical) new_location.normalize();
        cells.push_back(new_location);
    }

    cout << "Built " << cells.size() << " cells." << endl;
    delete[] xcell;
    delete[] ycell;
    delete[] zcell;

    // Build vertex location information
    xvertex = new double[nVertices];
    yvertex = new double[nVertices];
    zvertex = new double[nVertices];

    ncutil::get_var(inputFilename, "xVertex", xvertex);
    ncutil::get_var(inputFilename, "yVertex", yvertex);
    ncutil::get_var(inputFilename, "zVertex", zvertex);

    vertices.clear();
    for(int i = 0; i < nVertices; i++){
        xVertexRange[0] = min(xVertexRange[0], xvertex[i]);
        xVertexRange[1] = max(xVertexRange[1], xvertex[i]);

        yVertexRange[0] = min(yVertexRange[0], yvertex[i]);
        yVertexRange[1] = max(yVertexRange[1], yvertex[i]);

        zVertexRange[0] = min(zVertexRange[0], zvertex[i]);
        zVertexRange[1] = max(zVertexRange[1], zvertex[i]);

        new_location = pnt(xvertex[i], yvertex[i], zvertex[i], i);
        if(spherical) new_location.normalize();
        vertices.push_back(new_location);
    }

    cout << "Built " << vertices.size() << " vertices." << endl;
    delete[] xvertex;
    delete[] yvertex;
    delete[] zvertex;

    // Build unordered cellsOnVertex information
    cellsonvertex_list = new int[nVertices * vertexDegree];

    ncutil::get_var(inputFilename, "cellsOnVertex", cellsonvertex_list);
    cellsOnVertex.resize(nVertices);

    for(int i = 0; i < nVertices; i++){
        for(int j = 0; j < vertexDegree; j++){
            // Subtract 1 to convert into base 0 (c index space).
            cellsOnVertex.at(i).push_back(cellsonvertex_list[i*vertexDegree + j] - 1);
        }
    }

    delete[] cellsonvertex_list;

    meshDensity.clear();
    meshDensity.resize(cells.size());
    ncutil::get_var(inputFilename, "meshDensity", &meshDensity[0]);

    xCellDistance = fabs(xCellRange[1] - xCellRange[0]);
    yCellDistance = fabs(yCellRange[1] - yCellRange[0]);
    zCellDistance = fabs(zCellRange[1] - zCellRange[0]);

    xVertexDistance = fabs(xVertexRange[1] - xVertexRange[0]);
    yVertexDistance = fabs(yVertexRange[1] - yVertexRange[0]);
    zVertexDistance = fabs(zVertexRange[1] - zVertexRange[0]);

    if(periodic){
        // Quads are not staggered in the y direction, so need to period by
        // max distance + min distance
        if(vertexDegree == 4){
            if(xPeriod < 0.0){
                xPeriodicFix = xCellRange[0] + xCellRange[1];
            } else {
                xPeriodicFix = xPeriod;
            }

            if(yPeriod < 0.0){
                yPeriodicFix = yCellRange[0] + yCellRange[1];
            } else {
                yPeriodicFix = yPeriod;
            }
        // Triangles can be staggered, so only period my max distance
        } else {
            if(xPeriod < 0.0){
                xPeriodicFix = xCellRange[1];
            } else {
                xPeriodicFix = xPeriod;
            }

            if(yPeriod < 0.0){
                yPeriodicFix = yCellRange[1];
            } else {
                yPeriodicFix = yPeriod;
            }
        }
    } else {
        xPeriodicFix = 1.0e5 * abs(xCellRange[1] + xCellRange[0]);
        yPeriodicFix = 1.0e5 * abs(yCellRange[1] + yCellRange[0]);
    }

#ifdef _DEBUG
    xPeriod = xPeriodicFix;
    yPeriod = yPeriodicFix;
    cout << "cell Mins: " << xCellRange[0] << " " << yCellRange[0] << " " << zCellRange[0] << endl;
    cout << "vertex Mins: " << xVertexRange[0] << " " << yVertexRange[0] << " " << zVertexRange[0] << endl;
    cout << "cell Maxes: " << xCellRange[1] << " " << yCellRange[1] << " " << zCellRange[1] << endl;
    cout << "vertex Maxes: " << xVertexRange[1] << " " << yVertexRange[1] << " " << zVertexRange[1] << endl;
    cout << "cell Distances: " << xCellDistance << " " << yCellDistance << " " << zCellDistance << endl;
    cout << "vertex Distances: " << xVertexDistance << " " << yVertexDistance << " " << zVertexDistance << endl;
    cout << "xPeriodicFix: " << xPeriodicFix << endl;
    cout << "yPeriodicFix: " << yPeriodicFix << endl;
    cout << "xPeriod: " << xPeriod << endl;
    cout << "yPeriod: " << yPeriod << endl;
#endif

    if(!spherical && (zCellDistance > 0.0 || zVertexDistance > 0.0)){
        cout << "ERROR:" << endl;
        cout << "  This point set is defined in the plane, but has non-zero Z coordinates." << endl;
        cout << endl;
        return 1;
    }

    return 0;
}/*}}}*/
int buildUnorderedCellConnectivity(){/*{{{*/
    // buildUnorderedCellConnectivity should assume that cellsOnVertex hasn't been ordered properly yet.
    //    It should compute the inverse of cellsOnVertex (verticesOnCell) unordered.
    //    It will also compute an unordered cellsOnCell that can be considered invalid for actual use (e.g. quad grids).
    //    This ordering should happen regardless of it the mesh is planar or spherical.

    int iVertex, iCell, newCell, j, k, l, m, matches;
    bool add;

#ifdef _DEBUG
    cout << endl << endl << "Begin function: buildUnorderedCellConnectivity" << endl << endl;
#endif

    verticesOnCell.clear();
    verticesOnCell.resize(cells.size());

    for(iVertex = 0; iVertex < vertices.size(); iVertex++){
        for(j = 0; j < cellsOnVertex.at(iVertex).size(); j++){
            iCell = cellsOnVertex.at(iVertex).at(j);
            if(iCell != -1){
                add = true;
                for(int k = 0; k < verticesOnCell.at(iCell).size(); k++){
                    if(verticesOnCell.at(iCell).at(k) == iVertex){
                        add = false;
                    }
                }

                if(add) {
                    verticesOnCell.at(iCell).push_back(iVertex);
                }
            }
        }
    }


    cellsOnCell.clear();
    cellsOnCell.resize(cells.size());
    for(iCell = 0; iCell < cells.size(); iCell++){
        for(j = 0; j < verticesOnCell.at(iCell).size(); j++){
            iVertex = verticesOnCell.at(iCell).at(j);
            for(k = 0; k < cellsOnVertex.at(iVertex).size(); k++){
                newCell = cellsOnVertex.at(iVertex).at(k);

                if(newCell != iCell){
                    add = true;
                    for(l = 0; l < cellsOnCell.at(iCell).size(); l++){
                        if(cellsOnCell.at(iCell).at(l) == newCell){
                            add = false;
                        }
                    }

                    if(add) {
                        if(newCell != -1){
                            matches = 0;
                            for(l = 0; l < verticesOnCell.at(iCell).size(); l++){
                                for(m = 0; m < verticesOnCell.at(newCell).size(); m++){
                                    if(verticesOnCell.at(iCell).at(l) == verticesOnCell.at(newCell).at(m)){
                                        matches++;
                                    }
                                }
                            }

                            if(matches == 2){
                                cellsOnCell.at(iCell).push_back(newCell);
#ifdef _DEBUG
                                cout << "   Found two shared vertices for cell edge. Adding cell." << endl;
#endif
                            } else {
#ifdef _DEBUG
                                cout << "   Only found " << matches <<
                                    " shared vertices for cell edge. Not adding cell." << endl;
#endif

                            }
                        } else {
                            cellsOnCell.at(iCell).push_back(newCell);
                        }
                    }
                }
            }
        }
    }

    return 0;
}/*}}}*/
int firstOrderingVerticesOnCell(){/*{{{*/
    /*
     * firstOrderingVerticesOnCell should order the vertices around a cell such that they are connected.
     *        i.e. verticesOnCell.at(iCell).at(i) should be the tail of a vector pointing to
     *             verticesOnCell.at(iCell).at(i+1)
     *
     *        This is done by computing angles between vectors pointing from the cell center
     *        to each individual vertex. The next vertex in the list, has the smallest positive angle.
     *
     */

    pnt vec1, vec2, cross, normal;
    pnt vertex1, vertex2;
    int iCell, iVertex1, iVertex2, j, k, l, swp_idx, swp;
    double dot, mag1, mag2, angle, min_angle;

#ifdef _DEBUG
    cout << endl << endl << "Begin function: firstOrderingVerticesOnCell" << endl << endl;
#endif

    if(!spherical){
        normal = pnt(0.0, 0.0, 1.0);
    }

    for(iCell = 0; iCell < cells.size(); iCell++){
#ifdef _DEBUG
        cout << "new cell: " << iCell << endl;
        cout << "    " << cells.at(iCell) << endl;
#endif
        if(spherical){
            normal = cells.at(iCell);
        }

#ifdef _DEBUG
        cout << "  Unsorted verticesOnCell: ";
        for(j = 0; j < verticesOnCell.at(iCell).size(); j++){
            cout << verticesOnCell.at(iCell).at(j) << " ";
        }
        cout << endl;
        for(j = 0; j < verticesOnCell.at(iCell).size(); j++){
            cout << "  cellsOnVertex " << verticesOnCell.at(iCell).at(j) << ": ";
            for(k = 0; k < cellsOnVertex.at(verticesOnCell.at(iCell).at(j)).size(); k++){
                cout <<  cellsOnVertex.at(verticesOnCell.at(iCell).at(j)).at(k) << " ";
            }
            cout << endl;
        }
#endif

        // Loop over all vertices except the last two as a starting vertex
        for(j = 0; j < verticesOnCell.at(iCell).size()-1; j++){
            iVertex1 = verticesOnCell.at(iCell).at(j);

            vertex1 = vertices.at(iVertex1);

            //Check for, and fix up periodicity
            if(!spherical){
                vertex1.fixPeriodicity(cells.at(iCell), xPeriodicFix, yPeriodicFix);
            }
            vec1 = vertex1 - cells.at(iCell);

            mag1 = vec1.magnitude();
            min_angle = 2.0*M_PI;
            angle = 0.0;
            swp_idx = -1;

            // Don't sort any vertices that have already been sorted.
            for(k = j+1; k < verticesOnCell.at(iCell).size(); k++){
#ifdef _DEBUG
                cout << "    Comparing " << j << " " << k << endl;
#endif
                iVertex2 = verticesOnCell.at(iCell).at(k);

                vertex2 = vertices.at(iVertex2);

                //Check for, and fix up periodicity
                if(!spherical){
                    vertex2.fixPeriodicity(cells.at(iCell), xPeriodicFix, yPeriodicFix);
                }
                vec2 = vertex2 - cells.at(iCell);

                mag2 = vec2.magnitude();

                cross = vec1.cross(vec2);
                dot = cross.dot(normal) / (cross.magnitude() * normal.magnitude());
#ifdef _DEBUG
                cout << "    dot = " << dot << endl;
#endif

                // Only look at vertices that are CCW from the current vertex.
                if(dot > 0){
                    angle = acos(vec1.dot(vec2) / (mag1 * mag2));
#ifdef _DEBUG
                    cout << "    angle = " << angle << endl;
#endif
                    if(angle < min_angle){
                        min_angle = angle;
                        swp_idx = k;
                    }
                }
            }

            if(swp_idx != -1 && swp_idx != j+1){
                swp = verticesOnCell.at(iCell).at(j+1);
                verticesOnCell.at(iCell).at(j+1) = verticesOnCell.at(iCell).at(swp_idx);
                verticesOnCell.at(iCell).at(swp_idx) = swp;
            }

#ifdef _DEBUG
            if(swp_idx != -1 && swp_idx != j+1){
                cout << "      swapped " << j << " " << swp_idx << endl;
            } else {
                cout << "      no swap" << endl;
            }
            cout << "      cellsOnVertex(" << iVertex1 << "): ";
            for(k = 0; k < cellsOnVertex.at(iVertex1).size(); k++){
                cout << cellsOnVertex.at(iVertex1).at(k) << " ";
            }
            cout << endl;
            iVertex2 = verticesOnCell.at(iCell).at(j+1);
            cout << "      cellsOnVertex(" << iVertex2 << "): ";
            for(k = 0; k < cellsOnVertex.at(iVertex2).size(); k++){
                cout << cellsOnVertex.at(iVertex2).at(k) << " ";
            }
            cout << endl;
#endif
        }

#ifdef _DEBUG
        cout << "  Sorted verticesOnCell: ";
        for(j = 0; j < verticesOnCell.at(iCell).size(); j++){
            cout << verticesOnCell.at(iCell).at(j) << " ";
        }
        cout << endl;
        for(j = 0; j < verticesOnCell.at(iCell).size(); j++){
            cout << "  cellsOnVertex " << verticesOnCell.at(iCell).at(j) << ": ";
            for(k = 0; k < cellsOnVertex.at(verticesOnCell.at(iCell).at(j)).size(); k++){
                cout <<  cellsOnVertex.at(verticesOnCell.at(iCell).at(j)).at(k) << " ";
            }
            cout << endl;
        }
#endif
    }

    return 0;
}/*}}}*/
int buildCompleteCellMask(){/*{{{*/
    /*
     * The buildCompleteCellMask function parses an ordered
     * verticesOnCell field to determine if a cell is complete. It takes each
     * vertex-vertex pair (from verticesOnCell) and determines the angle
     * between them. If only one vertex is shared, the "edge" is skipped.
     *
     * These angles are summed up around each cell. If the total
     * angle is close to 2.0*Pi then the cell is marked as complete. Otherwise
     * it is marked as incomplete.
     *
     * The complete check is used when building edges. If an edge only has one
     * vertex, but is connected to two complete cells, the edge is not added to
     * the set. This resolves an issue related to quad grids, where a possible
     * edge connects two cells across a vertex (where the edge mid point would
     * be the vertex location). Using this check removes these edges from the
     * possible set of edges.
     *
     */
    int iCell, iCell2, vertex1, vertex2;
    int j, k, l;
    pnt vert_loc1, vert_loc2;
    pnt vec1, vec2;
    pnt cross, normal;
    double angle, angle_sum, dot;
    bool complete;

    completeCellMask.clear();
    completeCellMask.resize(cells.size());

    if(!spherical) {
        normal = pnt(0.0, 0.0, 1.0);
    }

    // Iterate over all cells
    for(iCell = 0; iCell < cells.size(); iCell++){
        complete = false;
        angle_sum = 0.0;

        if(spherical){
            normal = cells.at(iCell);
        }

#ifdef _DEBUG
        cout << "  Checking " << iCell << " for completeness" << endl;
        cout << "           " << cells.at(iCell) << endl;
#endif

        // Since verticesOnEdge is ordered already, compute angles between
        // neighboring vertex/vertex pairs if they are both valid vertices. Sum
        // the angles, and if the angles are close to 2.0 * Pi then the cell is
        // "complete"
        if(verticesOnCell.at(iCell).size() >= cellsOnCell.at(iCell).size()){
            for(j = 0; j < verticesOnCell.at(iCell).size()-1; j++){
                vertex1 = verticesOnCell.at(iCell).at(j);
                vertex2 = verticesOnCell.at(iCell).at(j+1);

                if(vertex1 != -1 && vertex2 != -1){
                    vert_loc1 = vertices.at(vertex1);
                    vert_loc2 = vertices.at(vertex2);

                    if(!spherical){
                        vert_loc1.fixPeriodicity(cells.at(iCell), xPeriodicFix, yPeriodicFix);
                        vert_loc2.fixPeriodicity(cells.at(iCell), xPeriodicFix, yPeriodicFix);
                    }

                    vec1 = vert_loc1 - cells.at(iCell);
                    vec2 = vert_loc2 - cells.at(iCell);

                    cross = vec1.cross(vec2);
                    dot = cross.dot(normal);

#ifdef _DEBUG
                    cout << "         vec1 : " << vec1 << endl;
                    cout << "         vec2 : " << vec2 << endl;
                    cout << "        cross : " << cross << endl;
                    cout << "          dot : " << dot << endl;
#endif

                    if(dot > 0){

                        angle = acos( vec2.dot(vec1) / (vec1.magnitude() * vec2.magnitude()));
                        angle_sum += angle;

#ifdef _DEBUG
                        cout << "        adding angle (rad) : " << angle << endl;
                        cout << "        adding angle (deg) : " << angle * 180.0 / M_PI << endl;
                        cout << "        new sum : " << angle_sum << endl;
#endif
                    } else {
                        cout << "     complete check has non CCW edge..." << endl;
                    }
                }
            }

            vertex1 = verticesOnCell.at(iCell).at(verticesOnCell.at(iCell).size() - 1);
            vertex2 = verticesOnCell.at(iCell).at(0);

            if(vertex1 != -1 && vertex2 != -1){
                vert_loc1 = vertices.at(vertex1);
                vert_loc2 = vertices.at(vertex2);

                if(!spherical){
                    vert_loc1.fixPeriodicity(cells.at(iCell), xPeriodicFix, yPeriodicFix);
                    vert_loc2.fixPeriodicity(cells.at(iCell), xPeriodicFix, yPeriodicFix);
                }

                vec1 = vert_loc1 - cells.at(iCell);
                vec2 = vert_loc2 - cells.at(iCell);

                cross = vec1.cross(vec2);
                dot = cross.dot(normal);

#ifdef _DEBUG
                cout << "         vec1 : " << vec1 << endl;
                cout << "         vec2 : " << vec2 << endl;
                cout << "        cross : " << cross << endl;
                cout << "          dot : " << dot << endl;
#endif

                if(dot > 0) {

                    angle = acos( vec2.dot(vec1) / (vec1.magnitude() * vec2.magnitude()));
                    angle_sum += angle;

#ifdef _DEBUG
                    cout << "        adding angle (rad) : " << angle << endl;
                    cout << "        adding angle (deg) : " << angle * 180.0 / M_PI << endl;
                    cout << "        new sum : " << angle_sum << endl;
#endif
                } else {
                    cout << "     complete check has non CCW edge..." << endl;
                }
            }
        }

        if(angle_sum > 2.0 * M_PI * 0.98){
            complete = true;
        }

#ifdef _DEBUG
        cout << "   Vertices on cell: " << verticesOnCell.at(iCell).size();
        cout << "   Cells on cell: " << cellsOnCell.at(iCell).size();
        if(complete){
            cout << "   Is complete!" << endl;
        } else {
            cout << "   Is not complete!" << endl;
        }
#endif

        if(complete){
            completeCellMask.at(iCell) = 1;
        } else {
            completeCellMask.at(iCell) = 0;
        }
    }

    return 0;
}/*}}}*/
int buildEdges(){/*{{{*/
    /*
     * buildEdges is intended to build a hash table that contains the
     *    cellsOnEdge and verticesonEdge pairs for each edge.  The actual edge
     *    location is not constructed at this point in time, as we need to
     *    determine the four constituent locations for each edge before we can
     *    compute it.
     *
     *    This function assumes the cellsOnVertex array is ordered such that
     *    two consecutive cell indices make up an edge.  This isn't an issue
     *    for triangular dual grids, but quadrilateral dual grids have an issue
     *    where every cell is not connected with every other cell. So, the tool
     *    that generates input data for this program needs to ensure that
     *    ordering.
     *
     *    After the hash table is created, it is iterated over to create each
     *    individual edge and ensure ordering is properly right handed.
     *
     */
    unordered_set<edge, edge::edge_hasher> edge_idx_hash;
    unordered_set<edge, edge::edge_hasher>::iterator edge_idx_itr;
    unordered_set<edge, edge::edge_hasher>::const_iterator edge_idx_search;

    int iCell, iVertex, iEdge, j, k, l;
    int cell1, cell2;
    int vertex1, vertex2, swp;
    int land;
    pair<unordered_set<edge,edge::edge_hasher>::iterator,bool> out_pair;
    edge new_edge;
    pnt edge_loc, normal;
    pnt cell_loc1, cell_loc2, vert_loc1, vert_loc2, dist_vec;
    pnt u_vec, v_vec, cross;
    double vert1_x_movement, vert1_y_movement;
    double vert2_x_movement, vert2_y_movement;
    double temp, dot;
    bool fixed_edge, add_edge;

#ifdef _DEBUG
    cout << endl << endl << "Begin function: buildEdges" << endl << endl;
#endif

    edge_idx_hash.clear();

    land = 0;

    // Build all edges
    for(iCell = 0; iCell < cells.size(); iCell++){
        if(completeCellMask.at(iCell) == 1) {
            // Build edges from every vertex/vertex pair around a cell if the cell is complete
            for(l = 0; l < verticesOnCell.at(iCell).size()-1; l++){
                vertex1 = verticesOnCell.at(iCell).at(l);
                vertex2 = verticesOnCell.at(iCell).at(l+1);

                // Find cell shaerd by vertices, that's not iCell
                cell1 = iCell;
                cell2 = -1;
                for(j = 0; j < cellsOnVertex.at(vertex1).size(); j++){
                    if(cellsOnVertex.at(vertex1).at(j) != -1 && cellsOnVertex.at(vertex1).at(j) != iCell){
                        for(k = 0; k < cellsOnVertex.at(vertex2).size(); k++){
                            if(cellsOnVertex.at(vertex1).at(j) == cellsOnVertex.at(vertex2).at(k)){
                                cell2 = cellsOnVertex.at(vertex1).at(j);
                            }
                        }
                    }
                }

                add_edge = true;

#ifdef _DEBUG
                cout << "1 Starting edge: " << endl;
                cout << "      vertex 1 : " << vertex1 << endl;
                cout << "      vertex 2 : " << vertex2 << endl;
                cout << "      cell 1 : " << cell1 << endl;
                cout << "      cell 2 : " << cell2 << endl;
#endif

                if(vertex1 != -1 && vertex2 != -1){
                    new_edge.vertex1 = min(vertex1, vertex2);
                    new_edge.vertex2 = max(vertex1, vertex2);

                    if(cell2 != -1){
                        new_edge.cell1 = min(cell1, cell2);
                        new_edge.cell2 = max(cell1, cell2);
                    } else {
                        new_edge.cell1 = cell1;
                        new_edge.cell2 = cell2;
                    }

                } else {
                    add_edge = false;
                }

                if(add_edge){
#ifdef _DEBUG
                    cout << " Adding edge" << endl;
#endif
                    out_pair = edge_idx_hash.insert(new_edge);

#ifdef _DEBUG
                    if(out_pair.second == false){
                        cout << " Failed to add edge." << endl;
                    }
#endif
                } else {
#ifdef _DEBUG
                    cout << " Not adding edge" << endl;
#endif
                }

            }

            vertex1 = verticesOnCell.at(iCell).at(verticesOnCell.at(iCell).size() - 1);
            vertex2 = verticesOnCell.at(iCell).at(0);

            // Find cell shaerd by vertices, that's not iCell
            cell1 = iCell;
            cell2 = -1;
            for(j = 0; j < cellsOnVertex.at(vertex1).size(); j++){
                if(cellsOnVertex.at(vertex1).at(j) != -1 && cellsOnVertex.at(vertex1).at(j) != iCell){
                    for(k = 0; k < cellsOnVertex.at(vertex2).size(); k++){
                        if(cellsOnVertex.at(vertex1).at(j) == cellsOnVertex.at(vertex2).at(k)){
                            cell2 = cellsOnVertex.at(vertex1).at(j);
                        }
                    }
                }
            }

            add_edge = true;

#ifdef _DEBUG
            cout << "  Starting edge: " << endl;
            cout << "      vertex 1 : " << vertex1 << endl;
            cout << "      vertex 2 : " << vertex2 << endl;
            cout << "      cell 1 : " << cell1 << endl;
            cout << "      cell 2 : " << cell2 << endl;
#endif

            if(vertex1 != -1 && vertex2 != -1){
                new_edge.vertex1 = min(vertex1, vertex2);
                new_edge.vertex2 = max(vertex1, vertex2);

                if(cell2 != -1){
                    new_edge.cell1 = min(cell1, cell2);
                    new_edge.cell2 = max(cell1, cell2);
                } else {
                    new_edge.cell1 = cell1;
                    new_edge.cell2 = cell2;
                }

            } else {
                add_edge = false;
            }

            if(add_edge){
#ifdef _DEBUG
                cout << " Adding edge" << endl;
#endif
                out_pair = edge_idx_hash.insert(new_edge);
#ifdef _DEBUG
                if(out_pair.second == false){
                    cout << " Failed to add edge." << endl;
                }
#endif
            } else {
#ifdef _DEBUG
                cout << " Not adding edge" << endl;
#endif
            }

        } else {
            // Build edges from every cell/cell pair only if cell is not complete
            for(l = 0; l < cellsOnCell.at(iCell).size(); l++){
                cell1 = iCell;
                cell2 = cellsOnCell.at(iCell).at(l);

                // Find vertex pair for cell pair
                vertex1 = -1;
                vertex2 = -1;

                for(j = 0; j < verticesOnCell.at(cell1).size(); j++){
                    if(cell2 != -1){
                        for(k = 0; k < verticesOnCell.at(cell2).size(); k++){
                            if(verticesOnCell.at(cell1).at(j) == verticesOnCell.at(cell2).at(k)) {
                                if(vertex1 == -1){
                                    vertex1 = verticesOnCell.at(cell1).at(j);
                                } else if(vertex2 == -1) {
                                    vertex2 = verticesOnCell.at(cell1).at(j);
                                } else {
                                    cout << " Found more than 2 vertices for edge? " << endl;
                                }
                            }
                        }
                    }
                }

                if(vertex2 == -1) {
                    new_edge.vertex1 = vertex1;
                    new_edge.vertex2 = vertex2;
                } else {
                    new_edge.vertex1 = min(vertex1, vertex2);
                    new_edge.vertex2 = max(vertex1, vertex2);
                }

#ifdef _DEBUG
                cout << "  Starting edge: " << endl;
                cout << "      vertex 1 : " << vertex1 << endl;
                cout << "      vertex 2 : " << vertex2 << endl;
                cout << "      cell 1 : " << cell1 << endl;
                cout << "      cell 2 : " << cell2 << endl;

                if(new_edge.vertex1 == -1){
                    cout << "  Edge is missing vertex 1" << endl;
                } else if(new_edge.vertex2 == -1){
                    cout << "  Edge is missing vertex 2" << endl;
                }
#endif

                if(cell2 == -1){
                    new_edge.cell1 = cell1;
                    new_edge.cell2 = -1;
                } else {
                    new_edge.cell1 = min(cell1, cell2);
                    new_edge.cell2 = max(cell1, cell2);
                }
                add_edge = true;

                if(new_edge.vertex1 == -1 || new_edge.vertex2 == -1){
                    if(new_edge.cell1 == -1 || new_edge.cell2 == -1){
                        add_edge = false;
                    }

#ifdef _DEBUG
                    cout << "   Cell 1 complete: " << completeCellMask.at(new_edge.cell1) << endl;
                    if(new_edge.cell2 != -1) {
                        cout << "   Cell 2 complete: " << completeCellMask.at(new_edge.cell2) << endl;
                    }
#endif

                    if(completeCellMask.at(new_edge.cell1) != 0 ||
                        (new_edge.cell2 != -1 && completeCellMask.at(new_edge.cell2) != 0) ){
                        add_edge = false;
                    }
                }


                if(add_edge){
#ifdef _DEBUG
                    cout << " Adding edge" << endl;
#endif
                    out_pair = edge_idx_hash.insert(new_edge);
#ifdef _DEBUG
                    if(out_pair.second == false){
                        cout << " Failed to add edge." << endl;
                    }
#endif
                } else {
#ifdef _DEBUG
                    cout << " Not adding edge" << endl;
#endif
                }
            }
        }
    }

    cout << "Built " << edge_idx_hash.size() << " edge indices..." << endl;

    cellsOnEdge.clear();
    verticesOnEdge.clear();
    dvEdge.clear();
    dcEdge.clear();
    cellsOnEdge.resize(edge_idx_hash.size());
    verticesOnEdge.resize(edge_idx_hash.size());
    dvEdge.resize(edge_idx_hash.size());
    dcEdge.resize(edge_idx_hash.size());

    if(!spherical){
        normal = pnt(0.0, 0.0, 1.0);
    }

    iEdge = 0;
    for(edge_idx_itr = edge_idx_hash.begin(); edge_idx_itr != edge_idx_hash.end(); ++edge_idx_itr){
#ifdef _DEBUG
        cout << "new edge: " << endl;
#endif
        fixed_edge = false;
        cell1 = (*edge_idx_itr).cell1;
        cell2 = (*edge_idx_itr).cell2;
        vertex1 = (*edge_idx_itr).vertex1;
        vertex2 = (*edge_idx_itr).vertex2;

        cell_loc1 = cells.at(cell1);
        vert_loc1 = vertices.at(vertex1);
        if(vertex2 != -1) {
            vert_loc2 = vertices.at(vertex2);
        } else {
            vert_loc2 = vert_loc1;
        }

        if(spherical){
            normal = cells.at(cell1);
        } else {
            // Clean up periodic edges. See mesh specification document for how periodicity is defined.
            // These edges are special in the sense that they must be across a "uniform" edge.
            // So, they are exactly the mid point between vertices.
            // Since the edge will be "owned" by whatever processor owns cell1, make sure edge is close to cell1.
            // So, fix periodicity relative to cell1.
            vert_loc1.fixPeriodicity(cell_loc1, xPeriodicFix, yPeriodicFix);
        }

        if(cell2 != -1){
            cell_loc2 = cells.at(cell2);
            if(vertex2 != -1){
                vert_loc2 = vertices.at(vertex2);

                if(!spherical) {
                    cell_loc2.fixPeriodicity(cell_loc1, xPeriodicFix, yPeriodicFix);
                    vert_loc2.fixPeriodicity(cell_loc1, xPeriodicFix, yPeriodicFix);
                    edge_loc = planarIntersect(cell_loc1, cell_loc2, vert_loc1, vert_loc2);
                    dvEdge.at(iEdge) = (vert_loc2 - vert_loc1).magnitude();
                    dcEdge.at(iEdge) = (cell_loc2 - cell_loc1).magnitude();
                } else { // If spherical
                    edge_loc = gcIntersect(cell_loc1, cell_loc2, vert_loc1, vert_loc2);
                    dvEdge.at(iEdge) = vert_loc2.sphereDistance(vert_loc1);
                    dcEdge.at(iEdge) = cell_loc2.sphereDistance(cell_loc1);
                }
            } else {
                edge_loc = (cell_loc1 + cell_loc2) * 0.5;
                vert_loc2 = edge_loc;
                if(!spherical){
                    dcEdge.at(iEdge) = (cell_loc2 - cell_loc1).magnitude();
                    dvEdge.at(iEdge) = dcEdge.at(iEdge) / sqrt(3.0);
                } else {
                    dcEdge.at(iEdge) = cell_loc2.sphereDistance(cell_loc1);
                    dvEdge.at(iEdge) = dcEdge.at(iEdge) / sqrt(3.0);
                }
            }
        } else {
            if(vertex2 != -1){
                vert_loc2 = vertices.at(vertex2);
                if(!spherical){
                    // Set dv and dc Edge values on plane
                    vert_loc2.fixPeriodicity(cell_loc1, xPeriodicFix, yPeriodicFix);
                    edge_loc = (vert_loc1 + vert_loc2) * 0.5;
                    dvEdge.at(iEdge) = (vert_loc2 - vert_loc1).magnitude();
                    dcEdge.at(iEdge) = sqrt(3.0) * dvEdge.at(iEdge);
                } else {
                    // Set dv and dc Edge values on Sphere
                    edge_loc = (vert_loc1 + vert_loc2) * 0.5;
                    edge_loc.normalize();
                    dvEdge.at(iEdge) = vert_loc2.sphereDistance(vert_loc1);
                    dcEdge.at(iEdge) = sqrt(3.0) * dvEdge.at(iEdge);
                }
                cell_loc2 = edge_loc;
            } else {
                cout << " ERROR: Edge found with only 1 cell and 1 vertex...." << endl;
                return 1;
            }
        }

        if(spherical) edge_loc.normalize();

        edge_loc.idx = iEdge;

#ifdef _DEBUG
        cout << "New Edge At: " << edge_loc << endl;
        cout << "         c1: " << cells.at(cell1) << endl;
        if(cell2 > -1) {
            cout << "         c2: " << cells.at(cell2) << endl;
            cout << "    mod? c2: " << cell_loc2 << endl;
        } else {
            cout << "         c2: land" << endl;
            cout << "    mod? c2: " << cell_loc2 << endl;
        }
        cout << "         v1: " << vertices.at(vertex1) << endl;
        cout << "    mod? v1: " << vert_loc1 << endl;
        if(vertex2 != -1) {
            cout << "         v2: " << vertices.at(vertex2) << endl;
            cout << "    mod? v2: " << vert_loc2 << endl;
        } else {
            cout << "         v2: land" << endl;
            cout << "    mod? v2: " << vert_loc2 << endl;
        }
#endif

        u_vec = cell_loc2 - cell_loc1;
        v_vec = vert_loc2 - vert_loc1;

        cross = u_vec.cross(v_vec);
        dot = cross.dot(normal);

        if(dot < 0){
            if(vertex2 != -1){
#ifdef _DEBUG
                cout << "   swapping vertex " << vertex1 << " and " << vertex2 << endl;
#endif
                swp = vertex2;
                vertex2 = vertex1;
                vertex1 = swp;
            } else {
#ifdef _DEBUG
                cout << "   swapping cell " << cell1 << " and " << cell2 << endl;
#endif
                swp = cell2;
                cell2 = cell1;
                cell1 = swp;
            }
        }

        edges.push_back(edge_loc);
        cellsOnEdge.at(iEdge).clear();
        cellsOnEdge.at(iEdge).push_back(cell1);
        cellsOnEdge.at(iEdge).push_back(cell2);
        verticesOnEdge.at(iEdge).clear();
        verticesOnEdge.at(iEdge).push_back(vertex1);
        verticesOnEdge.at(iEdge).push_back(vertex2);

        iEdge++;
    }

    edge_idx_hash.clear();

    return 0;
}/*}}}*/
int orderVertexArrays(){/*{{{*/
    /*
     * orderVertexArrays builds and orders the connectivity arrays for vertices.
     *        This includes edgesOnVertex and cellsOnvertex.
     *      First, edgeSOnVertex is built and ordered correctly, then using
     *      that ordering cellsOnVertex is built and ordered correctly as well.
     */
    int iEdge, iVertex, vertex1, vertex2, cell1, cell2, edge1, edge2;
    int i, j, k, l, swp_idx, swp;
    pnt normal;
    pnt cross;
    pnt vec1, vec2;
    pnt cell_loc, edge_loc1, edge_loc2;
    double mag1, mag2, min_angle, angle;
    double dvEdge, dot, area;
    bool fixed_area;

#ifdef _DEBUG
    cout << endl << endl << "Begin function: orderVertexArrays" << endl << endl;
#endif

    edgesOnVertex.clear();
    edgesOnVertex.resize(vertices.size());
    cellsOnVertex.clear();
    cellsOnVertex.resize(vertices.size());

    // Get lists of edges for each vertex
    for(iEdge = 0; iEdge < edges.size(); iEdge++){/*{{{*/
        vertex1 = verticesOnEdge.at(iEdge).at(0);
        vertex2 = verticesOnEdge.at(iEdge).at(1);

        if(vertex1 != -1) {
            edgesOnVertex.at(vertex1).push_back(iEdge);
        }

        if(vertex2 != -1){
            edgesOnVertex.at(vertex2).push_back(iEdge);
        }
    }/*}}}*/

    if(!spherical){
        normal = pnt(0.0, 0.0, 1.0);
    }

    // Order edges counter-clockwise
    for(iVertex = 0; iVertex < vertices.size(); iVertex++){/*{{{*/
        if(spherical){
            normal = vertices.at(iVertex);
        }

        //Loop over all edges except the last two as a starting edge.
        for(j = 0; j < edgesOnVertex.at(iVertex).size(); j++){
            edge1 = edgesOnVertex.at(iVertex).at(j);

            edge_loc1 = edges.at(edge1);

            // Fix edge1 periodicity
            if(!spherical){
                edge_loc1.fixPeriodicity(vertices.at(iVertex), xPeriodicFix, yPeriodicFix);
            }
            vec1 = edge_loc1 - vertices.at(iVertex);

            mag1 = vec1.magnitude();
            min_angle = 2.0*M_PI;
            angle = 0.0;
            swp_idx = -1;

            //Don't sort any edges that have already been sorted.
            for(k = j+1; k < edgesOnVertex.at(iVertex).size(); k++){
                edge2 = edgesOnVertex.at(iVertex).at(k);

                edge_loc2 = edges.at(edge2);

                // Fix edge2 periodicity
                if(!spherical){
                    edge_loc2.fixPeriodicity(vertices.at(iVertex), xPeriodicFix, yPeriodicFix);
                }
                vec2 = edge_loc2 - vertices.at(iVertex);
                mag2 = vec2.magnitude();

                cross = vec1.cross(vec2);
                dot = cross.dot(normal) / (cross.magnitude() * normal.magnitude());
                angle = acos(vec1.dot(vec2) / (mag1 * mag2));

                // Only look at vertices that are CCW from the current edge.
                if(dot > 0){
                    angle = acos(vec1.dot(vec2) / (mag1 * mag2));
                    if(angle < min_angle){
                        min_angle = angle;
                        swp_idx = k;
                    }
                }
            }

            if(swp_idx != -1 && swp_idx != j+1){
                swp = edgesOnVertex.at(iVertex).at(j+1);
                edgesOnVertex.at(iVertex).at(j+1) = edgesOnVertex.at(iVertex).at(swp_idx);
                edgesOnVertex.at(iVertex).at(swp_idx) = swp;
            }
        }

#ifdef _DEBUG
        cout << "edgesOnVertex("<< iVertex <<"): ";
        for(j = 0; j < edgesOnVertex.at(iVertex).size(); j++){
            cout << edgesOnVertex.at(iVertex).at(j) << " ";
        }
        cout << endl;
#endif

        cellsOnVertex.at(iVertex).clear();
#ifdef _DEBUG
        cout << "CellsOnVertex("<< iVertex << "): ";
#endif

        // Using the ordered edges. Buld cellsOnVertex in the correct order.
        for(j = 0; j < edgesOnVertex.at(iVertex).size(); j++){
            edge1 = edgesOnVertex.at(iVertex).at(j);
            fixed_area = false;

            // Get cell id and add it to list of cells
            if(iVertex == verticesOnEdge.at(edge1).at(0)){
                cell1 = cellsOnEdge.at(edge1).at(0);
                cellsOnVertex.at(iVertex).push_back( cellsOnEdge.at(edge1).at(0) );
#ifdef _DEBUG
                cout << cellsOnEdge.at(edge1).at(0) << " ";
#endif
            } else {
                cell1 = cellsOnEdge.at(edge1).at(1);
                cellsOnVertex.at(iVertex).push_back( cellsOnEdge.at(edge1).at(1) );
#ifdef _DEBUG
                cout << cellsOnEdge.at(edge1).at(1) << " ";
#endif
            }
        }
#ifdef _DEBUG
        cout << endl << endl;
#endif
    }/*}}}*/

    return 0;
}/*}}}*/
int orderCellArrays(){/*{{{*/
    /*
     * orderCellArrays assumes verticesOnCell are ordered CCW already.
     *
     */
    int iCell, iVertex, iEdge, iEdge2, add_cell;
    int cell1, cell2, vertex1, vertex2, prev_vertex;
    int edge_idx, loc_edge_idx;
    int i, j, k, swp_idx;
    pnt normal, cross;
    pnt vec1, vec2;
    pnt vert_loc, edge_loc, cell_loc;
    pnt edge_loc1, edge_loc2;
    pnt next_edge_loc;
    double swp;
    double dot, mag1, mag2, angle, min_angle;
    bool found, bad_vertices;

#ifdef _DEBUG
    cout << endl << endl << "Begin function: orderCellArrays" << endl << endl;
#endif

    maxEdges = 0;

    verticesOnCell.clear();
    edgesOnCell.clear();
    cellsOnCell.clear();

    verticesOnCell.resize(cells.size());
    edgesOnCell.resize(cells.size());
    cellsOnCell.resize(cells.size());

    if(!spherical){
        normal = pnt(0.0, 0.0, 1.0);
    }

    // First, build full list of edges on cell.
    for(iEdge = 0; iEdge < edges.size(); iEdge++){
        cell1 = cellsOnEdge.at(iEdge).at(0);
        cell2 = cellsOnEdge.at(iEdge).at(1);

        edgesOnCell.at(cell1).push_back(iEdge);
        if(cell2 != -1){
            edgesOnCell.at(cell2).push_back(iEdge);
        }
    }

    // Loop over all cells.
    for(iCell = 0; iCell < cells.size(); iCell++){
#ifdef _DEBUG
        cout << "New Cell " << cells.at(iCell) << endl;
#endif
        if(spherical){
            normal = cells.at(iCell);
        }

        loc_edge_idx = 0; // if no edges on cell, apparently occurs for quad meshes
        if ( edgesOnCell.at(iCell).size() != 0 ) {
#ifdef _DEBUG
        cout << endl;
            cout << "   Starting edgesOnCell on cell: " << iCell << endl;
            cout << "   Starting edgesOnCell size: " << edgesOnCell.at(iCell).size() << endl;
            cout << "   Starting edgesOnCell: ";
            for(i = 0; i < edgesOnCell.at(iCell).size(); ++i){
                cout << edgesOnCell.at(iCell).at(i) << " ";
            }
            cout << endl;
#endif

            // /*
            // Determine starting edge. It should either be the first edge in the set,
            // or the only edge such that all other edges are CCW from it.
            edge_idx = edgesOnCell.at(iCell).at(0);
    #ifdef _DEBUG
                    cout << "Finding starting edge for cell " << iCell << endl;
    #endif
            for(j = 0; j < edgesOnCell.at(iCell).size(); j++){
                iEdge = edgesOnCell.at(iCell).at(j);
                vertex1 = verticesOnEdge.at(iEdge).at(0);
                vertex2 = verticesOnEdge.at(iEdge).at(1);

                if(vertex2 == -1){
    #ifdef _DEBUG
                    cout << "   Start edge: " << iEdge << endl;
                    cout << "               " << edges.at(iEdge) << endl;
                    cout << "           v1: " << vertex1 << endl;
                    cout << "           v2: " << vertex2 << endl;
    #endif
                    edge_loc1 = edges.at(iEdge);
                    edge_loc1 = vertices.at(vertex1);
                    if(!spherical){
                        edge_loc1.fixPeriodicity(cells.at(iCell), xPeriodicFix, yPeriodicFix);
                    }

                    vec1 = edge_loc1 - cells.at(iCell);
                    mag1 = vec1.magnitude();

                    // If edge only has one vertex. Need to find edge that shares the vertex.
                    // This edge is kept if the neighboring edge is CCW from it.
                    for(k = 0; k < edgesOnCell.at(iCell).size(); k++){
                        if(j != k){
                            iEdge2 = edgesOnCell.at(iCell).at(k);
    #ifdef _DEBUG
                            cout << "    Test edge: " << iEdge2 << endl;
                            cout << "               " << edges.at(iEdge2) << endl;
                            cout << "           v1: " << verticesOnEdge.at(iEdge2).at(0) << endl;
                            cout << "           v2: " << verticesOnEdge.at(iEdge2).at(1) << endl;
    #endif

                            if(vertex1 == verticesOnEdge.at(iEdge2).at(0) || vertex1 == verticesOnEdge.at(iEdge2).at(1)){
                                // This edge is a neighboring edge. Check for CCW ordering.
                                if(vertex1 == verticesOnEdge.at(iEdge2).at(0)) {
                                    vertex2 = verticesOnEdge.at(iEdge2).at(1);
                                } else {
                                    vertex2 = verticesOnEdge.at(iEdge2).at(0);
                                }
                                edge_loc2 = edges.at(iEdge2);
                                edge_loc2 = vertices.at(vertex2);

                                if(!spherical){
                                    edge_loc2.fixPeriodicity(cells.at(iCell), xPeriodicFix, yPeriodicFix);
                                }

                                vec2 = edge_loc2 - cells.at(iCell);
                                mag2 = vec2.magnitude();

                                cross = vec1.cross(vec2);
                                dot = cross.dot(normal) / (cross.magnitude() * normal.magnitude());

    #ifdef _DEBUG
                                cout << "     Vec1: " << vec1 << endl;
                                cout << "     Vec2: " << vec2 << endl;
                                cout << "    Cross: " << cross << endl;
                                cout << "      dot: " << dot << endl;
    #endif
                                if(dot > 0){
    #ifdef _DEBUG
                                    cout << "       Found edge: " << iEdge << endl;
    #endif
                                    edge_idx = iEdge;
                                    loc_edge_idx = j;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Swap edge_idx with first edge.
        if(loc_edge_idx != 0){
#ifdef _DEBUG
            cout << "     Swapping edges: " <<
                edgesOnCell.at(iCell).at(loc_edge_idx) <<
                    " and " << edgesOnCell.at(iCell).at(0) << endl;
#endif
            edgesOnCell.at(iCell).at(loc_edge_idx) = edgesOnCell.at(iCell).at(0);
            edgesOnCell.at(iCell).at(0) = edge_idx;
        }
        // */

        // Order all edges in CCW relative to the first edge.
        // Loop over all vertices except the last two as a starting vertex.
        for(j = 0; j < edgesOnCell.at(iCell).size(); j++){
            iEdge = edgesOnCell.at(iCell).at(j);
            edge_loc = edges.at(iEdge);

            // Add cell and vertex from first edge.
            // Add cell across edge to cellsOnCell
            // Also, add vertex that is CCW relative to current edge location.
            if(cellsOnEdge.at(iEdge).at(0) == iCell){
                cellsOnCell.at(iCell).push_back(cellsOnEdge.at(iEdge).at(1));

                if(verticesOnEdge.at(iEdge).at(1) != -1){
                    verticesOnCell.at(iCell).push_back(verticesOnEdge.at(iEdge).at(1));
                }
            } else {
                cellsOnCell.at(iCell).push_back(cellsOnEdge.at(iEdge).at(0));
                verticesOnCell.at(iCell).push_back(verticesOnEdge.at(iEdge).at(0));
            }


            if(!spherical){
                edge_loc.fixPeriodicity(cells.at(iCell), xPeriodicFix, yPeriodicFix);
            }

            vec1 = edge_loc - cells.at(iCell);
            mag1 = vec1.magnitude();

            min_angle = 2.0 * M_PI;
            angle = 0.0;
            swp_idx = -1;

            for(k = j+1; k < edgesOnCell.at(iCell).size(); k++){
                iEdge2 = edgesOnCell.at(iCell).at(k);

                next_edge_loc = edges.at(iEdge2);

                if(!spherical){
                    next_edge_loc.fixPeriodicity(cells.at(iCell), xPeriodicFix, yPeriodicFix);
                }

                vec2 = next_edge_loc - cells.at(iCell);
                mag2 = vec2.magnitude();

                cross = vec1.cross(vec2);
                dot = cross.dot(normal) / (cross.magnitude() * normal.magnitude());

                // Only look at edges that are CCW from the current edge
                if(dot > 0){
                    angle = acos(vec1.dot(vec2) / (mag1 * mag2));
                    if(angle < min_angle){
                        min_angle = angle;
                        swp_idx = k;
                    }
                }
            }

            if(swp_idx != -1 && swp_idx != j+1){
                swp = edgesOnCell.at(iCell).at(j+1);
                edgesOnCell.at(iCell).at(j+1) = edgesOnCell.at(iCell).at(swp_idx);
                edgesOnCell.at(iCell).at(swp_idx) = swp;
            }
        }


#ifdef _DEBUG
        cout << "   cellsOnCell: ";
        for(i = 0; i < cellsOnCell.at(iCell).size(); i++){
            cout << cellsOnCell.at(iCell).at(i) << " ";
        }
        cout << endl;
        cout << "   verticesOnCell: ";
        for(i = 0; i < verticesOnCell.at(iCell).size(); i++){
            cout << verticesOnCell.at(iCell).at(i) << " ";
        }
        cout << endl;
        cout << "   edgesOnCell: ";
        for(i = 0; i < edgesOnCell.at(iCell).size(); i++){
            cout << edgesOnCell.at(iCell).at(i) << " ";
        }
        cout << endl;
#endif

        maxEdges = max(maxEdges, (int)edgesOnCell.at(iCell).size());
    }

    return 0;
}/*}}}*/
int buildAreas(){/*{{{*/
    /*
     * buildAreas constructs the area arrays.
     *    This includes areaCell, areaTriangle, and kiteAreasOnVertex
     *
     * Before constructing areaCell, the cell is checked for completeness.
     *    This is accomplished by looping over each edge. For each edge,
     *    the angle between it's vertices is computed and added to a running total.
     *    If the total is significantly less than 2*Pi, the cell is not complete.
     *
     * Non-complete cells are given a negative area, to ease removal at a later stage.
     */
    int iVertex, iCell, iEdge, i, j;
    int vertex1, vertex2;
    int edge1, edge2;
    int incomplete_cells;
    pnt vert_loc1, vert_loc2;
    pnt edge_loc1, edge_loc2;
    pnt cell_loc;
    pnt vec1, vec2;

#ifdef _DEBUG
    cout << endl << endl << "Begin function: buildAreas" << endl << endl;
#endif

    areaCell.clear();
    areaTriangle.clear();
    kiteAreasOnVertex.clear();

    areaCell.resize(cells.size());
    areaTriangle.resize(vertices.size());
    kiteAreasOnVertex.resize(vertices.size());

    incomplete_cells = 0;

    for(iCell = 0; iCell < cells.size(); iCell++){
        areaCell.at(iCell) = 0.0;

        if(completeCellMask.at(iCell) == 1){
            for(j = 0; j < edgesOnCell.at(iCell).size(); j++){
                iEdge = edgesOnCell.at(iCell).at(j);

                if(cellsOnEdge.at(iEdge).at(0) == iCell){
                    vertex1 = verticesOnEdge.at(iEdge).at(0);
                    vertex2 = verticesOnEdge.at(iEdge).at(1);
                } else {
                    vertex1 = verticesOnEdge.at(iEdge).at(1);
                    vertex2 = verticesOnEdge.at(iEdge).at(0);
                }

                if(vertex1 != -1 && vertex2 != -1){
                    vert_loc1 = vertices.at(vertex1);
                    vert_loc2 = vertices.at(vertex2);

                    if(!spherical){
                        vert_loc1.fixPeriodicity(cells.at(iCell), xPeriodicFix, yPeriodicFix);
                        vert_loc2.fixPeriodicity(cells.at(iCell), xPeriodicFix, yPeriodicFix);
                        areaCell.at(iCell) += planarTriangleArea(cells.at(iCell), vert_loc1, vert_loc2);
                    } else {
                        areaCell.at(iCell) += sphericalTriangleArea(cells.at(iCell), vert_loc1, vert_loc2);
                    }
                }
            }
        } else {
            incomplete_cells++;
#ifdef _DEBUG
            cout << "   Non complete cell found at " << iCell << endl;
#endif
            areaCell.at(iCell) = -1;
        }
    }

    for(iVertex = 0; iVertex < vertices.size(); iVertex++){
        areaTriangle.at(iVertex) = 0.0;
        kiteAreasOnVertex.at(iVertex).resize(cellsOnVertex.at(iVertex).size());
        for(j = 0; j < cellsOnVertex.at(iVertex).size(); j++){
            kiteAreasOnVertex.at(iVertex).at(j) = 0.0;
            iCell = cellsOnVertex.at(iVertex).at(j);

            if(iCell != -1){
                edge1 = edgesOnVertex.at(iVertex).at(j);

                if(j == cellsOnVertex.at(iVertex).size()-1){
                    edge2 = edgesOnVertex.at(iVertex).at(0);
                } else {
                    edge2 = edgesOnVertex.at(iVertex).at(j+1);
                }

                cell_loc = cells.at(iCell);
                edge_loc1 = edges.at(edge1);
                edge_loc2 = edges.at(edge2);

                if(!spherical){
                    cell_loc.fixPeriodicity(vertices.at(iVertex), xPeriodicFix, yPeriodicFix);
                    edge_loc1.fixPeriodicity(vertices.at(iVertex), xPeriodicFix, yPeriodicFix);
                    edge_loc2.fixPeriodicity(vertices.at(iVertex), xPeriodicFix, yPeriodicFix);
                    kiteAreasOnVertex.at(iVertex).at(j) +=
                        planarTriangleArea(vertices.at(iVertex), edge_loc1, cell_loc);
                    kiteAreasOnVertex.at(iVertex).at(j) +=
                        planarTriangleArea(vertices.at(iVertex), cell_loc, edge_loc2);
                } else {
                    kiteAreasOnVertex.at(iVertex).at(j) +=
                        sphericalTriangleArea(vertices.at(iVertex), edge_loc1, cell_loc);
                    kiteAreasOnVertex.at(iVertex).at(j) +=
                        sphericalTriangleArea(vertices.at(iVertex), cell_loc, edge_loc2);
                }

                areaTriangle.at(iVertex) += kiteAreasOnVertex.at(iVertex).at(j);
            }
        }
    }

    cout << "     Found " << incomplete_cells << " incomplete cells. Each is marked with an area of -1." << endl;

    return 0;
}/*}}}*/
int buildEdgesOnEdgeArrays(){/*{{{*/
    /*
     * buildEdgesOnEdgeArrays builds both edgesOnEdge and weightsOnEdge
     *
     * The weights correspond to edgesOnEdge and allow MPAS to reconstruct the
     * edge perpendicular (previously defined as v) velocity using the edge
     * neighbors
     *
     * Weight formulation is defined in J. Thurburn, et al. JCP 2009
     * Numerical representation of geostrophic modes on arbitrarily
     * structured C-grids
     *
     * Before using this function, dcEdge, dvEdge, areaCell, and kiteAreasOnVertex
     * all need to be computed correctly.
     */
    int iEdge, iCell, cell1, cell2;
    int i, j, k;
    int shared_vertex;
    int last_edge, cur_edge;
    double area_sum;
    bool found;

#ifdef _DEBUG
    cout << endl << endl << "Begin function: buildEdgesOnEdgeArrays" << endl << endl;
#endif

    edgesOnEdge.clear();
    edgesOnEdge.resize(edges.size());

    weightsOnEdge.clear();
    weightsOnEdge.resize(edges.size());

    for(iEdge = 0; iEdge < edges.size(); iEdge++){
#ifdef _DEBUG
        cout << "New edge: " << edges.at(iEdge) << endl;
#endif
        cell1 = cellsOnEdge.at(iEdge).at(0);
        cell2 = cellsOnEdge.at(iEdge).at(1);
        found = false;

        // Loop over cell 1. Starting from the edge after the current edge, add
        // all edges counter clockwise around the cell.  Don't add the current
        // edge to the list.
#ifdef _DEBUG
        cout << "   On cell1: " << cell1 << endl;
#endif
        last_edge = iEdge;
        area_sum = 0;
        for(i = 0; i < edgesOnCell.at(cell1).size(); i++){
#ifdef _DEBUG
            cout << "      checking edge: " << edgesOnCell.at(cell1).at(i) << endl;
#endif
            if(edgesOnCell.at(cell1).at(i) == iEdge){
                found = true;
#ifdef _DEBUG
                cout << "        -- found -- " << endl;
#endif
            }

            if(found && edgesOnCell.at(cell1).at(i) != iEdge){
                cur_edge = edgesOnCell.at(cell1).at(i);
                edgesOnEdge.at(iEdge).push_back(cur_edge);

                // Find shared vertex between newly added edge
                // and last_edge
                if(verticesOnEdge.at(last_edge).at(0) == verticesOnEdge.at(cur_edge).at(0) ||
                        verticesOnEdge.at(last_edge).at(0) == verticesOnEdge.at(cur_edge).at(1)){

                    shared_vertex = verticesOnEdge.at(last_edge).at(0);

                } else if(verticesOnEdge.at(last_edge).at(1) == verticesOnEdge.at(cur_edge).at(0) ||
                        verticesOnEdge.at(last_edge).at(1) == verticesOnEdge.at(cur_edge).at(1)){

                    shared_vertex = verticesOnEdge.at(last_edge).at(1);
                }

                if(shared_vertex != -1) {
                    // Find cell 1 on shared vertex (to get kite area)
                    for(j = 0; j < cellsOnVertex.at(shared_vertex).size(); j++){
                        iCell = cellsOnVertex.at(shared_vertex).at(j);

                        if(iCell == cell1){
                            area_sum += kiteAreasOnVertex.at(shared_vertex).at(j) / areaCell.at(cell1);
                        }
                    }
                }

                if(cell1 == cellsOnEdge.at(cur_edge).at(0)){
                    weightsOnEdge.at(iEdge).push_back(
                        +1.0 * (0.5 - area_sum) * dvEdge.at(cur_edge) / dcEdge.at(iEdge) );
                } else {
                    weightsOnEdge.at(iEdge).push_back(
                        -1.0 * (0.5 - area_sum) * dvEdge.at(cur_edge) / dcEdge.at(iEdge) );
                }

                last_edge = edgesOnCell.at(cell1).at(i);
#ifdef _DEBUG
                cout << "        added " << edgesOnCell.at(cell1).at(i) << endl;
#endif
            }
        }
        if(!found) {
            cout << "Error finding edge " << iEdge << " on cell " << cell1 << " 1" << endl;
            return 1;
        }

        for(i = 0; i < edgesOnCell.at(cell1).size() && found; i++){
#ifdef _DEBUG
            cout << "      checking edge: " << edgesOnCell.at(cell1).at(i) << endl;
#endif
            if(edgesOnCell.at(cell1).at(i) == iEdge){
#ifdef _DEBUG
                cout << "        -- found -- " << endl;
#endif
                found = false;
            }

            if(found && edgesOnCell.at(cell1).at(i) != iEdge){
                cur_edge = edgesOnCell.at(cell1).at(i);
                edgesOnEdge.at(iEdge).push_back(cur_edge);

                // Find shared vertex between newly added edge
                // and last_edge
                if(verticesOnEdge.at(last_edge).at(0) == verticesOnEdge.at(cur_edge).at(0) ||
                        verticesOnEdge.at(last_edge).at(0) == verticesOnEdge.at(cur_edge).at(1)){

                    shared_vertex = verticesOnEdge.at(last_edge).at(0);

                } else if(verticesOnEdge.at(last_edge).at(1) == verticesOnEdge.at(cur_edge).at(0) ||
                        verticesOnEdge.at(last_edge).at(1) == verticesOnEdge.at(cur_edge).at(1)){

                    shared_vertex = verticesOnEdge.at(last_edge).at(1);
                }

                // Find cell 1 on shared vertex (to get kite area)
                if(shared_vertex != -1){
                    for(j = 0; j < cellsOnVertex.at(shared_vertex).size(); j++){
                        iCell = cellsOnVertex.at(shared_vertex).at(j);

                        if(iCell == cell1){
                            area_sum += kiteAreasOnVertex.at(shared_vertex).at(j) / areaCell.at(cell1);
                        }
                    }
                }

                if(cell1 == cellsOnEdge.at(cur_edge).at(0)){
                    weightsOnEdge.at(iEdge).push_back(
                        +1.0 * (0.5 - area_sum) * dvEdge.at(cur_edge) / dcEdge.at(iEdge) );
                } else {
                    weightsOnEdge.at(iEdge).push_back(
                        -1.0 * (0.5 - area_sum) * dvEdge.at(cur_edge) / dcEdge.at(iEdge) );
                }

                last_edge = edgesOnCell.at(cell1).at(i);
#ifdef _DEBUG
                cout << "        added " << edgesOnCell.at(cell1).at(i) << endl;
#endif
            }
        }
        if(found) {
            cout << "Error finding edge " << iEdge << " on cell " << cell1 << " 1" << endl;
            return 1;
        }

        // Check if cell1 is a real cell or not.
        if(cell2 > -1){
            last_edge = iEdge;
            area_sum = 0.0;
#ifdef _DEBUG
            cout << "   On cell2: " << cell2 << endl;
#endif
            found = false;
            // Loop over cell 2. Starting from the edge after the current edge,
            // add all edges counter clockwise around the cell.  Don't add the
            // current edge to the list.
            for(i = 0; i < edgesOnCell.at(cell2).size(); i++){
#ifdef _DEBUG
                cout << "      checking edge: " << edgesOnCell.at(cell2).at(i) << endl;
#endif
                if(edgesOnCell.at(cell2).at(i) == iEdge){
                    found = true;
#ifdef _DEBUG
                    cout << "        -- found -- " << endl;
#endif
                }

                if(found && edgesOnCell.at(cell2).at(i) != iEdge){
                    cur_edge = edgesOnCell.at(cell2).at(i);
                    edgesOnEdge.at(iEdge).push_back(cur_edge);

                    // Find shared vertex between newly added edge
                    // and last_edge
                    if(verticesOnEdge.at(last_edge).at(0) == verticesOnEdge.at(cur_edge).at(0) ||
                            verticesOnEdge.at(last_edge).at(0) == verticesOnEdge.at(cur_edge).at(1)){

                        shared_vertex = verticesOnEdge.at(last_edge).at(0);

                    } else if(verticesOnEdge.at(last_edge).at(1) == verticesOnEdge.at(cur_edge).at(0) ||
                            verticesOnEdge.at(last_edge).at(1) == verticesOnEdge.at(cur_edge).at(1)){

                        shared_vertex = verticesOnEdge.at(last_edge).at(1);
                    }

                    // Find cell 1 on shared vertex (to get kite area)
                    if(shared_vertex != -1){
                        for(j = 0; j < cellsOnVertex.at(shared_vertex).size(); j++){
                            iCell = cellsOnVertex.at(shared_vertex).at(j);

                            if(iCell == cell2){
                                area_sum += kiteAreasOnVertex.at(shared_vertex).at(j) / areaCell.at(cell2);
                            }
                        }
                    }

                    if(cell2 == cellsOnEdge.at(cur_edge).at(0)){
                        weightsOnEdge.at(iEdge).push_back(
                            -1.0 * (0.5 - area_sum) * dvEdge.at(cur_edge) / dcEdge.at(iEdge) );
                    } else {
                        weightsOnEdge.at(iEdge).push_back(
                            +1.0 * (0.5 - area_sum) * dvEdge.at(cur_edge) / dcEdge.at(iEdge) );
                    }

                    last_edge = edgesOnCell.at(cell2).at(i);
#ifdef _DEBUG
                    cout << "        added " << edgesOnCell.at(cell2).at(i) << endl;
#endif
                }
            }
            if(!found) {
                cout << "Error finding edge " << iEdge << " on cell " << cell2 << " 2" << endl;
                return 1;
            }

            for(i = 0; i < edgesOnCell.at(cell2).size(); i++){
#ifdef _DEBUG
                cout << "      checking edge: " << edgesOnCell.at(cell2).at(i) << endl;
#endif
                if(edgesOnCell.at(cell2).at(i) == iEdge){
                    found = false;
#ifdef _DEBUG
                    cout << "        -- found -- " << endl;
#endif
                }

                if(found && edgesOnCell.at(cell2).at(i) != iEdge){
                    cur_edge = edgesOnCell.at(cell2).at(i);
                    edgesOnEdge.at(iEdge).push_back(cur_edge);

                    // Find shared vertex between newly added edge
                    // and last_edge
                    if(verticesOnEdge.at(last_edge).at(0) == verticesOnEdge.at(cur_edge).at(0) ||
                            verticesOnEdge.at(last_edge).at(0) == verticesOnEdge.at(cur_edge).at(1)){

                        shared_vertex = verticesOnEdge.at(last_edge).at(0);

                    } else if(verticesOnEdge.at(last_edge).at(1) == verticesOnEdge.at(cur_edge).at(0) ||
                            verticesOnEdge.at(last_edge).at(1) == verticesOnEdge.at(cur_edge).at(1)){

                        shared_vertex = verticesOnEdge.at(last_edge).at(1);
                    }

                    // Find cell 1 on shared vertex (to get kite area)
                    if(shared_vertex != -1){
                        for(j = 0; j < cellsOnVertex.at(shared_vertex).size(); j++){
                            iCell = cellsOnVertex.at(shared_vertex).at(j);

                            if(iCell == cell2){
                                area_sum += kiteAreasOnVertex.at(shared_vertex).at(j) / areaCell.at(cell2);
                            }
                        }
                    }

                    if(cell2 == cellsOnEdge.at(cur_edge).at(0)){
                        weightsOnEdge.at(iEdge).push_back(
                            -1.0 * (0.5 - area_sum) * dvEdge.at(cur_edge) / dcEdge.at(iEdge) );
                    } else {
                        weightsOnEdge.at(iEdge).push_back(
                            +1.0 * (0.5 - area_sum) * dvEdge.at(cur_edge) / dcEdge.at(iEdge) );
                    }

                    last_edge = edgesOnCell.at(cell2).at(i);
#ifdef _DEBUG
                    cout << "        added " << edgesOnCell.at(cell2).at(i) << endl;
#endif
                }
            }
            if(found) {
                cout << "Error finding edge " << iEdge << " on cell " << cell2 << " 2" << endl;
                return 1;
            }
        }
    }

    return 0;
}/*}}}*/
int buildAngleEdge(){/*{{{*/
    /*
     * buildAngleEdge constructs angle edge for each edge.
     *
     * angleEdge is either:
     * 1. The angle the positive tangential direction (v)
     *    makes with the local northward direction.
     *    or
     * 2. The angles the positive normal direction (u)
     *       makes with the local eastward direction.
     *
     *       In a plane, local eastward is defined as the x axis, and
     *       nortward is defined as the y axis.
     */
    int iEdge;
    int cell1, cell2;
    int vertex1, vertex2;
    pnt np, x_axis, normal;
    pnt cell_loc1, cell_loc2;
    pnt vertex_loc1, vertex_loc2;
    double angle, sign;

#ifdef _DEBUG
    cout << endl << endl << "Begin function: buildAngleEdge" << endl << endl;
#endif

    angleEdge.clear();
    angleEdge.resize(edges.size());

    x_axis = pnt(1.0, 0.0, 0.0);

    for(iEdge = 0; iEdge < edges.size(); iEdge++){

#ifdef _DEBUG
        cout << "New edge: " << edges.at(iEdge) << endl;
#endif
        cell1 = cellsOnEdge.at(iEdge).at(0);
        cell2 = cellsOnEdge.at(iEdge).at(1);

        vertex1 = verticesOnEdge.at(iEdge).at(0);
        vertex2 = verticesOnEdge.at(iEdge).at(1);

        cell_loc1 = cells.at(cell1);
        if(cell2 != -1){
            cell_loc2 = cells.at(cell2);
        } else {
            cell_loc2 = edges.at(iEdge);
        }

        vertex_loc1 = vertices.at(vertex1);
        if(vertex2 != -1){
            vertex_loc2 = vertices.at(vertex2);
        } else {
            vertex_loc2 = edges.at(iEdge);
        }

        if(!spherical){
            cell_loc2.fixPeriodicity(cell_loc1, xPeriodicFix, yPeriodicFix);

            normal = cell_loc2 - cell_loc1;
            angleEdge.at(iEdge) = acos( x_axis.dot(normal) / (x_axis.magnitude() * normal.magnitude()));
            if (cell_loc2.y < cell_loc1.y) angleEdge.at(iEdge) = 2.0 * M_PI - angleEdge.at(iEdge);
        } else {

            np = pntFromLatLon(edges.at(iEdge).getLat()+0.05, edges.at(iEdge).getLon());
            np.normalize();

#ifdef _DEBUG
            cout << "     NP: " << np << endl;
#endif

            angle = (vertex_loc2.getLat() - vertex_loc1.getLat()) / dvEdge.at(iEdge);
            angle = max( min(angle, 1.0), -1.0);
            angle = acos(angle);

#ifdef _DEBUG
            cout << "  angle: " << angle << endl;
#endif

            sign = planeAngle(edges.at(iEdge), np, vertex_loc2, edges.at(iEdge));
            if(sign != 0.0){
                sign = sign / fabs(sign);
            } else {
                sign = 1.0;
            }


#ifdef _DEBUG
            cout << "  sign : " << sign << endl;
            cout << "  a*s  : " << angle * sign << endl;
#endif

            angle = angle * sign;
            if(angle > M_PI) angle = angle - 2.0 * M_PI;
            if(angle < -M_PI) angle = angle + 2.0 * M_PI;

#ifdef _DEBUG
            cout << " fangle: " << angle << endl;
#endif

            angleEdge.at(iEdge) = angle;
        }
    }

    return 0;
}/*}}}*/
int buildMeshQualities(){/*{{{*/
    /*
     * buildMeshQualities constructs fields describing the quality of the mesh, including:
     *        - cellQuality: a double between 0 and 1 describing how uniform the cell is
     *        - gridSpacing: an estimate of the grid spacing of each cell
     *        - triangleQuality: a double between 0 and 1 describing how uniform the cell is based on edge lenghts
     *        - triangleAngleQuality: a double between 0 and 1 describing how uniform the cell is based on the angles of the triangle
     *        - obtuseTriangle: an integer either 0 or 1 that flags triangles with an obtuse angle (set to 1)
     */

    int iCell, iVertex, iEdge;
    int i, j;
    double maxEdge, minEdge, spacing;
    double angle, maxAngle, minAngle;
    int obtuse;

    cellQuality.clear();
    gridSpacing.clear();
    triangleQuality.clear();
    triangleAngleQuality.clear();
    obtuseTriangle.clear();

    cellQuality.resize(cells.size());
    gridSpacing.resize(cells.size());

    triangleQuality.resize(vertices.size());
    triangleAngleQuality.resize(vertices.size());
    obtuseTriangle.resize(vertices.size());

#ifdef _DEBUG
    cout << "Starting mesh quality calcs." <<endl;
#endif

    for ( iCell = 0; iCell < cells.size(); iCell++ ) {
        maxEdge = 0.0;
        minEdge = DBL_MAX;
        spacing = 0.0;

        for ( j = 0; j < edgesOnCell.at(iCell).size(); j++ ) {
            iEdge = edgesOnCell.at(iCell).at(j);

            maxEdge = max(maxEdge, dvEdge.at(iEdge));
            minEdge = min(minEdge, dvEdge.at(iEdge));

            spacing += dcEdge.at(iEdge);
        }
        if (maxEdge > 0.0) {
            cellQuality.at(iCell) = minEdge / maxEdge;
        } else {
            cellQuality.at(iCell) = 0.0;
        }
        gridSpacing.at(iCell) = spacing / edgesOnCell.at(iCell).size();
    }

#ifdef _DEBUG
    cout << "Starting loop over vertices." <<endl;
#endif

    for ( iVertex = 0; iVertex < vertices.size(); iVertex++ ) {
        maxEdge = 0.0;
        minEdge = DBL_MAX;
        maxAngle = 0.0;
        minAngle = DBL_MAX;
        obtuse = 0;

#ifdef _DEBUG
    cout << "vertex:" << iVertex <<endl;
#endif

        for ( j = 0; j < edgesOnVertex.at(iVertex).size(); j++ ) {
#ifdef _DEBUG
            cout << "  edge:" << j <<endl;
#endif
            iEdge = edgesOnVertex.at(iVertex).at(j);

            maxEdge = max(maxEdge, dcEdge.at(iEdge));
            minEdge = min(minEdge, dcEdge.at(iEdge));
            triangleQuality.at(iVertex) = minEdge / maxEdge;

#ifdef _DEBUG
            cout << "  triangleQuality:" << minEdge / maxEdge <<endl;
#endif

            if ( vertex_degree == 3 ) {
                double a_len, b_len, c_len;
                double angle1, angle2, angle3;

                if (edgesOnVertex.at(iVertex).size() == 3) {
                    a_len = dcEdge.at(edgesOnVertex.at(iVertex).at(0));
#ifdef _DEBUG
                    cout << "  length 1:" << a_len <<endl;
#endif

                    b_len = dcEdge.at(edgesOnVertex.at(iVertex).at(1));
#ifdef _DEBUG
                    cout << "  length 2:" << b_len<< endl;
#endif

                    c_len = dcEdge.at(edgesOnVertex.at(iVertex).at(2));
#ifdef _DEBUG
                    cout << "  length 3:" << c_len <<endl;
#endif

                    angle1 = acos( max(-1.0, min(1.0,
                        (b_len * b_len + c_len * c_len - a_len * a_len) / (2 * b_len * c_len))));
                    angle2 = acos( max(-1.0, min(1.0,
                        (a_len * a_len + c_len * c_len - b_len * c_len) / (2 * a_len * c_len))));
                    angle3 = acos( max(-1.0, min(1.0,
                        (a_len * a_len + b_len * b_len - c_len * c_len) / (2 * a_len * b_len))));

                    minAngle = min(angle1, min(angle2, angle3));
                    maxAngle = max(angle1, max(angle2, angle3));

                    if ( maxAngle > M_PI_2 ) {
                        obtuse = 1;
                        obtuseTriangles++;
                    }

                    triangleAngleQuality.at(iVertex) = minAngle / maxAngle;
                    obtuseTriangle.at(iVertex) = obtuse;
                } else {
                    triangleAngleQuality.at(iVertex) = 1.0;
                    obtuseTriangle.at(iVertex) = 0;
                }
            }
        }
    }

    cout << "\tMesh contains: " << obtuseTriangles << " obtuse triangles." << endl;

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

    int grid, retv;
    if ((retv = nc_create(outputFilename.c_str(), NC_CLOBBER|NC_NETCDF4, &grid)))
    {
        std::cout << "Can't create file: " << outputFilename << std::endl;
        return retv ;
    }

    nCells = cells.size();

    /*
    for(vec_int_itr = edgesOnCell.begin(); vec_int_itr != edgesOnCell.end(); ++vec_int_itr){
        maxEdges = std::max(maxEdges, (int)(*vec_int_itr).size());
    }*/

    ncutil::def_dim(outputFilename, "nCells", cells.size());
    ncutil::def_dim(outputFilename, "nEdges", edges.size());
    ncutil::def_dim(outputFilename, "nVertices", vertices.size());
    ncutil::def_dim(outputFilename, "maxEdges", maxEdges);
    ncutil::def_dim(outputFilename, "maxEdges2", maxEdges * 2);
    ncutil::def_dim(outputFilename, "TWO", 2);
    ncutil::def_dim(outputFilename, "vertexDegree", vertex_degree);
    ncutil::def_dim(outputFilename, "Time", NC_UNLIMITED);

    return 0;
}/*}}}*/
int outputGridAttributes( const string outputFilename, const string inputFilename ){/*{{{*/
    /************************************************************************
     *
     * This function writes the grid dimensions to the netcdf file named
     * outputFilename
     *
     * **********************************************************************/

    char mesh_spec_str[1024];

    string history_str = "";
    string id_str = "";
    string parent_str ="";

    // write attributes
    if(!spherical){
        ncutil::put_str(outputFilename, "on_a_sphere", "NO");
        ncutil::put_att(outputFilename, "sphere_radius", NC_DOUBLE, 0.);
    } else {
        ncutil::put_str(outputFilename, "on_a_sphere", "YES");
        ncutil::put_att(outputFilename,
            "sphere_radius", NC_DOUBLE, sphereRadius);
    }

    if(!periodic){
        ncutil::put_str(outputFilename, "is_periodic", "NO");
    } else {
        ncutil::put_str(outputFilename, "is_periodic", "YES");
        ncutil::put_att(outputFilename, "x_period", NC_DOUBLE, xPeriod);
        ncutil::put_att(outputFilename, "y_period", NC_DOUBLE, yPeriod);
    }

    history_str += "MpasMeshConverter.x ";
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

    sprintf(mesh_spec_str, "%2.1lf", (double)MESH_SPEC);

    ncutil::put_str(outputFilename, "history", history_str);
    ncutil::put_str(outputFilename, "mesh_spec", mesh_spec_str);
    ncutil::put_str(outputFilename, "Conventions", "MPAS");
    ncutil::put_str(outputFilename, "source", "MpasMeshConverter.x");
    ncutil::put_str(outputFilename, "file_id", id_str);

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

    int nCells = cells.size();
    int nEdges = edges.size();
    int nVertices = vertices.size();

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
    for(pnt_itr = cells.begin(); pnt_itr != cells.end(); ++pnt_itr){
        if(!spherical){
            x[i] = (*pnt_itr).x;
            y[i] = (*pnt_itr).y;
            z[i] = (*pnt_itr).z;
            lat[i] = 0.0;
            lon[i] = 0.0;
        } else {
            x[i] = (*pnt_itr).x * sphereRadius;
            y[i] = (*pnt_itr).y * sphereRadius;
            z[i] = (*pnt_itr).z * sphereRadius;
            lat[i] = (*pnt_itr).getLat();
            lon[i] = (*pnt_itr).getLon();
        }
        idxTo[i] = (*pnt_itr).idx+1;

        i++;
    }

    ncutil::def_var(outputFilename, "latCell",
        NC_DOUBLE, "latitudes of cell centres", {"nCells"});
    ncutil::def_var(outputFilename, "lonCell",
        NC_DOUBLE, "longitudes of cell centres", {"nCells"});

    ncutil::put_var(outputFilename, "latCell", &lat[0]);
    ncutil::put_var(outputFilename, "lonCell", &lon[0]);

    ncutil::def_var(outputFilename, "xCell",
        NC_DOUBLE, "x-coordinates of cell centres", {"nCells"});
    ncutil::def_var(outputFilename, "yCell",
        NC_DOUBLE, "y-coordinates of cell centres", {"nCells"});
    ncutil::def_var(outputFilename, "zCell",
        NC_DOUBLE, "z-coordinates of cell centres", {"nCells"});

    ncutil::put_var(outputFilename, "xCell", &x[0]);
    ncutil::put_var(outputFilename, "yCell", &y[0]);
    ncutil::put_var(outputFilename, "zCell", &z[0]);

    ncutil::def_var(outputFilename, "indexToCellID",
        NC_INT, "index to cell ID mapping", {"nCells"});

    ncutil::put_var(outputFilename, "indexToCellID", &idxTo[0]);

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] lat;
    delete[] lon;
    delete[] idxTo;

    //Build and write edge coordinate arrays
    x = new double[nEdges];
    y = new double[nEdges];
    z = new double[nEdges];
    lat = new double[nEdges];
    lon = new double[nEdges];
    idxTo = new int[nEdges];

    i = 0;
    for(pnt_itr = edges.begin(); pnt_itr != edges.end(); ++pnt_itr){
        if(!spherical){
            x[i] = (*pnt_itr).x;
            y[i] = (*pnt_itr).y;
            z[i] = (*pnt_itr).z;
            lat[i] = 0.0;
            lon[i] = 0.0;
        } else {
            x[i] = (*pnt_itr).x * sphereRadius;
            y[i] = (*pnt_itr).y * sphereRadius;
            z[i] = (*pnt_itr).z * sphereRadius;
            lat[i] = (*pnt_itr).getLat();
            lon[i] = (*pnt_itr).getLon();
        }
        idxTo[i] = (*pnt_itr).idx+1;

        i++;
    }

    ncutil::def_var(outputFilename, "latEdge",
        NC_DOUBLE, "latitudes of edge centres", {"nEdges"});
    ncutil::def_var(outputFilename, "lonEdge",
        NC_DOUBLE, "longitudes of edge centres", {"nEdges"});

    ncutil::put_var(outputFilename, "latEdge", &lat[0]);
    ncutil::put_var(outputFilename, "lonEdge", &lon[0]);

    ncutil::def_var(outputFilename, "xEdge",
        NC_DOUBLE, "x-coordinates of edge centres", {"nEdges"});
    ncutil::def_var(outputFilename, "yEdge",
        NC_DOUBLE, "y-coordinates of edge centres", {"nEdges"});
    ncutil::def_var(outputFilename, "zEdge",
        NC_DOUBLE, "z-coordinates of edge centres", {"nEdges"});

    ncutil::put_var(outputFilename, "xEdge", &x[0]);
    ncutil::put_var(outputFilename, "yEdge", &y[0]);
    ncutil::put_var(outputFilename, "zEdge", &z[0]);

    ncutil::def_var(outputFilename, "indexToEdgeID",
        NC_INT, "index to edge ID mapping", {"nEdges"});

    ncutil::put_var(outputFilename, "indexToEdgeID", &idxTo[0]);

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] lat;
    delete[] lon;
    delete[] idxTo;

    //Build and write vertex coordinate arrays
    x = new double[nVertices];
    y = new double[nVertices];
    z = new double[nVertices];
    lat = new double[nVertices];
    lon = new double[nVertices];
    idxTo = new int[nVertices];

    i = 0;
    for(pnt_itr = vertices.begin(); pnt_itr != vertices.end(); ++pnt_itr){
        if(!spherical){
            x[i] = (*pnt_itr).x;
            y[i] = (*pnt_itr).y;
            z[i] = (*pnt_itr).z;
            lat[i] = 0.0;
            lon[i] = 0.0;
        } else {
            x[i] = (*pnt_itr).x * sphereRadius;
            y[i] = (*pnt_itr).y * sphereRadius;
            z[i] = (*pnt_itr).z * sphereRadius;
            lat[i] = (*pnt_itr).getLat();
            lon[i] = (*pnt_itr).getLon();
        }
        idxTo[i] = (*pnt_itr).idx+1;

        i++;
    }

    ncutil::def_var(outputFilename, "latVertex",
        NC_DOUBLE, "latitudes of vertices", {"nVertices"});
    ncutil::def_var(outputFilename, "lonVertex",
        NC_DOUBLE, "longitudes of vertices", {"nVertices"});

    ncutil::put_var(outputFilename, "latVertex", &lat[0]);
    ncutil::put_var(outputFilename, "lonVertex", &lon[0]);

    ncutil::def_var(outputFilename, "xVertex",
        NC_DOUBLE, "x-coordinates of vertices", {"nVertices"});
    ncutil::def_var(outputFilename, "yVertex",
        NC_DOUBLE, "y-coordinates of vertices", {"nVertices"});
    ncutil::def_var(outputFilename, "zVertex",
        NC_DOUBLE, "z-coordinates of vertices", {"nVertices"});

    ncutil::put_var(outputFilename, "xVertex", &x[0]);
    ncutil::put_var(outputFilename, "yVertex", &y[0]);
    ncutil::put_var(outputFilename, "zVertex", &z[0]);

    ncutil::def_var(outputFilename, "indexToVertexID",
        NC_INT, "index to vertex ID mapping", {"nVertices"});

    ncutil::put_var(outputFilename, "indexToVertexID", &idxTo[0]);

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] lat;
    delete[] lon;
    delete[] idxTo;

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

    int nCells = cells.size();
    int i, j;

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

    ncutil::def_var(outputFilename, "cellsOnCell",
        NC_INT, "cells adj. to each cell", {"nCells", "maxEdges"});

    ncutil::put_var(outputFilename, "cellsOnCell", &tmp_arr[0]);

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

    ncutil::def_var(outputFilename, "edgesOnCell",
        NC_INT, "edges on each cell", {"nCells", "maxEdges"});

    ncutil::put_var(outputFilename, "edgesOnCell", &tmp_arr[0]);

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

    ncutil::def_var(outputFilename, "verticesOnCell",
        NC_INT, "vertices on each cell", {"nCells", "maxEdges"});

    ncutil::put_var(outputFilename, "verticesOnCell", &tmp_arr[0]);

    delete[] tmp_arr;

    //Build and write nEOC array
    tmp_arr = new int[nCells];

    i = 0;
    for(vec_int_itr = edgesOnCell.begin(); vec_int_itr != edgesOnCell.end(); ++vec_int_itr){
        tmp_arr[i] = (*vec_int_itr).size();
        i++;
    }

    ncutil::def_var(outputFilename, "nEdgesOnCell",
        NC_INT, "number of edges on each cell", {"nCells"});

    ncutil::put_var(outputFilename, "nEdgesOnCell", &tmp_arr[0]);

    verticesOnCell.clear();
    edgesOnCell.clear();

    delete[] tmp_arr;

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

    int nEdges = edges.size();
    int maxEdges2 = maxEdges * 2;
    int vertexDegree = vertex_degree;
    int two = 2;
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

    ncutil::def_var(outputFilename, "edgesOnEdge",
        NC_INT, "edges adj. to each edge", {"nEdges", "maxEdges2"});

    ncutil::put_var(outputFilename, "edgesOnEdge", &tmp_arr[0]);

    delete[] tmp_arr;

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

    ncutil::def_var(outputFilename, "cellsOnEdge",
        NC_INT, "cells adj. to each edge", {"nEdges", "TWO"});

    ncutil::put_var(outputFilename, "cellsOnEdge", &tmp_arr[0]);

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

    ncutil::def_var(outputFilename, "verticesOnEdge",
        NC_INT, "vertices on each edge", {"nEdges", "TWO"});

    ncutil::put_var(outputFilename, "verticesOnEdge", &tmp_arr[0]);

    delete[] tmp_arr;

    // Build and write nEoe array
    tmp_arr = new int[nEdges];
    i = 0;
    for(vec_int_itr = edgesOnEdge.begin(); vec_int_itr != edgesOnEdge.end(); ++vec_int_itr){
        tmp_arr[i] = (*vec_int_itr).size();
        i++;
    }

    ncutil::def_var(outputFilename, "nEdgesOnEdge",
        NC_INT, "number of edges on each edge", {"nEdges"});

    ncutil::put_var(outputFilename, "nEdgesOnEdge", &tmp_arr[0]);

    delete[] tmp_arr;

    cellsOnEdge.clear();
//    verticesOnEdge.clear(); // Needed for Initial conditions.
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

    int nVertices = vertices.size();
    int vertexDegree = vertex_degree;
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

    ncutil::def_var(outputFilename, "cellsOnVertex",
        NC_INT, "vertices adj. to each vertex", {"nVertices", "vertexDegree"});

    ncutil::put_var(outputFilename, "cellsOnVertex", &tmp_arr[0]);

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

    ncutil::def_var(outputFilename, "edgesOnVertex",
        NC_INT, "edges adj. to each vertex", {"nVertices", "vertexDegree"});

    ncutil::put_var(outputFilename, "edgesOnVertex", &tmp_arr[0]);

    delete[] tmp_arr;

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

    ncutil::def_var(outputFilename, "boundaryVertex",
        NC_INT, "non-zero for each vertex on mesh boundary", {"nVertices"});

    ncutil::put_var(outputFilename, "boundaryVertex", &tmp_arr[0]);

    delete[] tmp_arr;

    cellsOnVertex.clear();
    edgesOnVertex.clear();

    return 0;
}/*}}}*/
int outputCellParameters( const string outputFilename) {/*{{{*/
    /*********************************************************
     *
     * This function writes all cell parameters, including
     *     areaCell
     *
     * *******************************************************/

    int nCells = cells.size();
    int i, j;

    if(spherical){
        for(i = 0; i < nCells; i ++){
            areaCell[i] = areaCell[i] * sphereRadius * sphereRadius;
        }
    }

    ncutil::def_var(outputFilename, "areaCell",
        NC_DOUBLE, "surface areas of cells", {"nCells"});

    ncutil::put_var(outputFilename, "areaCell", &areaCell[0]);

    areaCell.clear();

    return 0;
}/*}}}*/
int outputVertexParameters( const string outputFilename) {/*{{{*/
    /*********************************************************
     *
     * This function writes all vertex parameters, including
     *     areaTriangle
     *     kiteAreasOnVertex
     *
     * *******************************************************/

    int nVertices = vertices.size();
    int vertexDegree = vertex_degree;
    int i, j;
    int count;

    double *tmp_arr;

    if(spherical){
        for(i = 0; i < nVertices; i++){
            areaTriangle[i] = areaTriangle[i] * sphereRadius * sphereRadius;
        }
    }

    ncutil::def_var(outputFilename, "areaTriangle",
        NC_DOUBLE, "surface areas of dual cells", {"nVertices"});

    ncutil::put_var(outputFilename, "areaTriangle", &areaTriangle[0]);

    // Build and write kiteAreasOnVertex
    // TODO: Fix kite area for quads?
    tmp_arr = new double[nVertices*vertexDegree];
    for(i = 0; i < nVertices*vertexDegree; i++){
        tmp_arr[i] = 0.0;
    }
    i = 0;
    count = 0;
    for(vec_dbl_itr = kiteAreasOnVertex.begin(); vec_dbl_itr != kiteAreasOnVertex.end(); ++vec_dbl_itr){
        j = 0;
        for(dbl_itr = (*vec_dbl_itr).begin(); dbl_itr != (*vec_dbl_itr).end(); ++dbl_itr){
            if(!spherical){
                tmp_arr[i*vertexDegree + j] = (*dbl_itr);
            } else {
                tmp_arr[i*vertexDegree + j] = (*dbl_itr) * sphereRadius * sphereRadius;
            }
            j++;
            count++;
        }
        i++;
    }

    ncutil::def_var(outputFilename, "kiteAreasOnVertex",
        NC_DOUBLE,
    "surface areas of overlap between cells and dual cells", {"nVertices", "vertexDegree"});

    ncutil::put_var(outputFilename, "kiteAreasOnVertex", &tmp_arr[0]);

    delete[] tmp_arr;

    kiteAreasOnVertex.clear();

    return 0;
}/*}}}*/
int outputEdgeParameters( const string outputFilename) {/*{{{*/
    /*********************************************************
     *
     * This function writes all grid parameters, including
     *    angleEdge
     *    dcEdge
     *    dvEdge
     *    weightsOnEdge
     *
     * *******************************************************/

    int nEdges = edges.size();
    int maxEdges2 = maxEdges * 2;
    int i, j;

    double *tmp_arr;

    if(spherical){
        for(i = 0; i < nEdges; i++){
            dcEdge[i] = dcEdge[i] * sphereRadius;
            dvEdge[i] = dvEdge[i] * sphereRadius;
        }
    }

    ncutil::def_var(outputFilename, "angleEdge",
        NC_DOUBLE, "angle to edges", {"nEdges"}) ;

    ncutil::put_var(outputFilename, "angleEdge", &angleEdge[0]);

    ncutil::def_var(outputFilename, "dcEdge",
        NC_DOUBLE, "length of arc between centres", {"nEdges"});
    ncutil::def_var(outputFilename, "dvEdge",
        NC_DOUBLE, "length of arc between vertices", {"nEdges"});

    ncutil::put_var(outputFilename, "dcEdge", &dcEdge[0]) ;
    ncutil::put_var(outputFilename, "dvEdge", &dvEdge[0]) ;

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

    ncutil::def_var(outputFilename, "weightsOnEdge",
        NC_DOUBLE, "tangential flux reconstruction weights", {"nEdges", "maxEdges2"});

    ncutil::put_var(outputFilename, "weightsOnEdge", &tmp_arr[0]);

    delete[] tmp_arr;

    angleEdge.clear();
    dcEdge.clear();
    dvEdge.clear();
    weightsOnEdge.clear();

    return 0;
}/*}}}*/
int outputMeshDensity( const string outputFilename) {/*{{{*/
    /***************************************************************************
     *
     * This function writes the meshDensity variable. Read in from the file SaveDensity
     *
     * *************************************************************************/

    ncutil::def_var(outputFilename, "meshDensity",
        NC_DOUBLE, "mesh density distribution", {"nCells"});

    ncutil::put_var(outputFilename, "meshDensity", &meshDensity[0]);

    return 0;
}/*}}}*/
int outputMeshQualities( const string outputFilename) {/*{{{*/
    /***************************************************************************
     *
     * This function writes the mesh quality variables.
     *        - cellQuality
     *        - gridSpacing
     *        - triangleQuality
     *        - triangleAnalgeQuality
     *        - obtuseTriangle
     *
     * *************************************************************************/

    ncutil::def_var(outputFilename, "cellQuality",
        NC_DOUBLE, "quality of mesh cells", {"nCells"});

    ncutil::put_var(outputFilename, "cellQuality", &cellQuality[0]);

    ncutil::def_var(outputFilename, "gridSpacing",
        NC_DOUBLE, "grid spacing distribution", {"nCells"});

    ncutil::put_var(outputFilename, "gridSpacing", &cellQuality[0]);

    ncutil::def_var(outputFilename, "triangleQuality",
        NC_DOUBLE, "quality of mesh dual cells", {"nVertices"});

    ncutil::put_var(outputFilename, "triangleQuality", &triangleQuality[0]);

    ncutil::def_var(outputFilename, "triangleAngleQuality",
        NC_DOUBLE, "quality of mesh dual cells", {"nVertices"});

    ncutil::put_var(outputFilename,
        "triangleAngleQuality", &triangleAngleQuality [0]);

    ncutil::def_var(outputFilename, "obtuseTriangle",
        NC_INT,
    "non-zero for any dual cell containing obtuse angles", {"nVertices"});

    ncutil::put_var(outputFilename, "obtuseTriangle", &obtuseTriangle[0]);

    cellQuality.clear();
    gridSpacing.clear();
    triangleQuality.clear();
    triangleAngleQuality.clear();
    obtuseTriangle.clear();

    return 0;
}/*}}}*/
int writeGraphFile(const string outputFilename){/*{{{*/
    ofstream graph(outputFilename.c_str());

    int edgeCount = 0;

    for (int iCell = 0; iCell < nCells; iCell++){
        for ( int j = 0; j < cellsOnCell.at(iCell).size(); j++){
            int coc = cellsOnCell.at(iCell).at(j);

            if ( coc >= 0 && coc < nCells ) {
                edgeCount++;
            }
        }
    }

    edgeCount = edgeCount / 2;
    graph << cells.size() << " " << edgeCount << endl;

    for(int i = 0; i < cellsOnCell.size(); i++){
        for(int j = 0; j < cellsOnCell.at(i).size(); j++){
            if ( cellsOnCell.at(i).at(j) >= 0 ) {
                graph << cellsOnCell.at(i).at(j)+1 << " ";
            }
        }
        graph << endl;
    }

    graph.close();

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

