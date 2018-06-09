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
#include <float.h>
#include <json/json.h>
#include <omp.h>
#include <string.h>
#include <tuple>
#include <sstream>
#include <sys/stat.h>

#include "netcdf_utils.h"
#include "pnt.h"
#include "edge.h"

#define MESH_SPEC 1.0
#define ID_LEN 10

//#define _WRITE_POINT_MASK

using namespace std;

bool lonRangePositive = false;
int nCells, nVertices, nEdges;
int maxEdges, vertexDegree;
int nRegions = 0;
int nTransects = 0;
int nPoints = 0;
int maxRegPolygons = 0;
int maxRegVertices = 0;
int maxRegionsInGroup = 0;
int maxPointsInGroup = 0;
int maxTrnLines = 0;
int maxTrnVertices = 0;
int maxTransectsInGroup = 0;
int maxCellsInTransect = 0;
int maxVerticesInTransect = 0;
int maxEdgesInTransect = 0;
double sphereRadius;
bool spherical, periodic;
string in_history = "";
string in_file_id = "";
string in_parent_id = "";

enum types { num_int, num_double, text };

// Mask and location information {{{

// Seed information
vector<pnt> seedPoints;
vector<int> cellSeedMask;

// Region information:
int *cellMasks, *vertexMasks;
vector< vector< vector<pnt> > > regionPolygons;
vector<string> regionNames;
vector< vector< tuple<string, string, int>> > regionProperties;
vector< vector<int> > regionsInGroup;
vector<string> regionGroupNames;
vector< vector< vector<double> > > polygonConstants;
vector< vector< vector<double> > > polygonMultiples;

// Transect information
vector< vector< vector<pnt> > > transectPoints;
vector<string> transectNames;
vector< vector< tuple<string, string, int>> > transectProperties;
vector< vector<int> > transectsInGroup;
vector<string> transectGroupNames;
vector< vector< vector<int> > > closestIndexToTransectPoint;
vector< vector<int> > cellPaths, edgePaths, vertexPaths;
vector< vector<int> > edgePathSigns;

// Point information:
int *pointCellIndices, *pointVertexIndices;
int *cellPointMasks;
vector<pnt> pointLocations;
vector<string> pointNames;
vector< vector< tuple<string, string, int> > > pointProperties;
vector< vector<int> > pointsInGroup;
vector<string> pointGroupNames;

// Mesh locations:
vector<pnt> cells;
vector<pnt> vertices;
vector<pnt> edges;
vector< vector< pair<int, double> > > connectivityGraph;

// }}}

// Iterators {{{
vector<pnt>::iterator pnt_itr, vert_itr;
vector<string>::iterator str_itr;
vector< vector<int> >::iterator mask_itr, vec_int_itr;
vector< vector< vector<pnt> > >::iterator reg_itr;
vector< vector<pnt> >::iterator poly_itr;
vector<int>::iterator int_itr;
// }}}

/* Utility functions {{{*/
int featureIndex( const string featureName, const vector<string> featureNames );
int stringType ( const string testStr );
vector< tuple<string, string, int> > extractFeatureProperties( Json::Value properties );
vector<int> walkGraph( const vector< vector< pair<int, double> > > graph, const int start, const int end );
void buildMaskFromFloodFill( const vector<pnt> seedLocations, const vector<pnt> staticLocations, const vector< vector< pair<int, double> > > graph, vector<int> *mask );
/*}}}*/

/* Building/Reading functions {{{ */
int readGridInfo(const string inputFilename);
int readCells(const string inputFilename);
int readVertices(const string inputFilename);
int readEdges(const string inputFilename);
vector< vector< pair<int, double> > > readCellGraph(const string inputFilename);
vector< vector< pair<int, double> > > readVertexGraph(const string inputFilename);
vector< vector< pair<int, double> > > readEdgeGraph(const string inputFilename);
int resetFeatureInfo();
int getFeatureInfo(const string featureFilename);
int getSeedInfo(const string seedFilename);
vector< vector< vector<int> > > buildClosestValuesForTransects(const vector< vector< vector<pnt> > > lineSegments, const vector<pnt> staticLocations);
vector< vector<int> > buildLinePaths( const vector< vector< vector<int> > > closestIndices, const vector< vector< pair<int, double> > > graph );
vector< vector<int> > buildEdgePathSigns( const string inputFilename, const vector< vector<int> > edgePaths );
int buildPolygonValues();
int buildMasks(vector<pnt> locations, int *masks);
int buildPointIndices(vector<pnt> testLocations, vector<pnt> staticLocations, int *indices);
int buildAllFeatureGroups();
/*}}}*/

/* Output functions {{{*/
int outputMaskDimensions( const string outputFilename );
int outputMaskAttributes( const string outputFilename, const string inputFilename );
int outputMaskFields( const string outputFilename);
/*}}}*/

void print_usage() {/*{{{*/
	cout << endl << endl;
	cout << " USAGE:" << endl;
	cout << "\tMpasMaskCreator.x in_file out_file [ [-f/-s] file.geojson ] [--positive_lon]" << endl;
	cout << "\t\tin_file: This argument defines the input file that masks will be created for." << endl;
	cout << "\t\tout_file: This argument defines the file that masks will be written to." << endl;
	cout << "\t\t-s file.geojson: This argument pair defines a set of points (from the geojson point definition)" << endl;
	cout << "\t\t\tthat will be used as seed points in a flood fill algorithim. This is useful when trying to remove isolated cells from a mesh." << endl;
	cout << "\t\t-f file.geojson: This argument pair defines a set of geojson features (regions, transects, or points)" << endl;
	cout << "\t\t\tthat will be converted into masks / lists." << endl;
	cout << "\t\t--positive_lon: This argument causes the logitude range to be 0-360 degrees with the prime meridian at 0 degrees." << endl;
	cout << "\t\t\tIf this flag is not set, the logitude range is -180-180 with 0 degrees being the prime meridian." << endl;
	cout << "\t\t\tWhether this flag is passed in or not, any longitudes written are in the 0-360 range." << endl;
}/*}}}*/

string gen_random(const int len);

int main ( int argc, char *argv[] ) {
	int error;
	vector<string> mask_files;
	vector<string> seed_files;
	string out_name = "masks.nc";
	string in_name = "grid.nc";

	cout << endl << endl;
	cout << "************************************************************" << endl;
	cout << "MPAS_MASK_CREATOR:\n";
	cout << "  C++ version\n";
	cout << "  Creates a set of masks for a given MPAS mesh and a set of feature files. " << endl;
	cout << endl << endl;
	cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
	cout << "************************************************************" << endl;
	cout << endl << endl;

	if ( argc < 5 ) {
		cout << " ERROR: Incorrect usage. See usage statement." << endl;
		print_usage();
		exit(1);
	} else {
		string str_flag;
		string str_file;
		in_name = argv[1];
		out_name = argv[2];

		int i = 3;
		while ( i < argc ) {
			str_flag = argv[i];

			if ( str_flag == "-s" ) {
				str_file = argv[i+1];
				seed_files.push_back(str_file);
				i += 2;
			} else if ( str_flag == "-f" ) {
				str_file = argv[i+1];
				mask_files.push_back(str_file);
				i += 2;
			} else if ( str_flag == "--positive_lon" ) {
				lonRangePositive = true;
				i++;
			} else {
				cout << " ERROR: Invalid flag " << str_flag << " passed in. See usage statement." << endl;
				print_usage();
				exit(1);
			}
		}
	}

	if(in_name == out_name){
		cout << "   ERROR: Input and output names are the same." << endl;
		return 1;
	}

	// Check if input files exist
	{
		bool error;
		struct stat buffer;

		error = false;

		if ( stat( in_name.c_str(), &buffer ) != 0 ){
			cout << "ERROR: Input file '" << in_name << "' does not exist." << endl;
			error = true;
		}

		for ( str_itr = mask_files.begin(); str_itr != mask_files.end(); str_itr++ ) {
			if ( stat( (*str_itr).c_str(), &buffer ) != 0 ) {
				cout << "ERROR: Input file '" << (*str_itr) << "' does not exist." << endl;
				error = true;
			}
		}

		for ( str_itr = seed_files.begin(); str_itr != seed_files.end(); str_itr++ ) {
			if ( stat( (*str_itr).c_str(), &buffer ) != 0 ) {
				cout << "ERROR: Input file '" << (*str_itr) << "' does not exist." << endl;
				error = true;
			}
		}

		if ( error ) {
			cout << "ERROR: One or more missing input files. Exiting..." << endl;
			exit(1);
		}


	}

	// Write out the longitude range that is being used
	cout << endl;
	if ( lonRangePositive ) {
		cout << "Using 0 to 360 degrees for longitude range." << endl;
	} else {
		cout << "Using -180 to 180 degrees for longitude range." << endl;
	}
	cout << endl;
	cout << endl;

	srand(time(NULL));

	cout << "Reading input grid." << endl;
	error = readGridInfo(in_name);
	if(error) exit(1);

	error = resetFeatureInfo();

	if ( mask_files.size() > 0 ) {
		cout << "Building feature information." << endl;
		for ( vector<string>::iterator mask_itr = mask_files.begin(); mask_itr != mask_files.end(); mask_itr++ ) {
			error = getFeatureInfo( (*mask_itr) );
			if(error) exit(1);
		}
	}

	if ( seed_files.size() > 0 ) {
		cout << "Building seed locations." << endl;
		for ( vector<string>::iterator seed_itr = seed_files.begin(); seed_itr != seed_files.end(); seed_itr++ ) {
			error = getSeedInfo( (*seed_itr) );
			if(error) exit(1);
		}
	}

	cout << "Building 'all' feature groups." << endl;
	if ( error = buildAllFeatureGroups()){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Building polygon edge lines" << endl;
	if(error = buildPolygonValues()){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Reading cell center locations" << endl;
	if(error = readCells(in_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Marking cells based on region definitions" << endl;
	cellMasks = new int[cells.size() * regionPolygons.size()];
	if(error = buildMasks(cells, &cellMasks[0])){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Marking points based on cell centers" << endl;
	pointCellIndices = new int[pointLocations.size()];
	if(error = buildPointIndices(pointLocations, cells, &pointCellIndices[0])){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << " Building closest point location for line segments" << endl;
	closestIndexToTransectPoint = buildClosestValuesForTransects(transectPoints, cells);

	cout << " Reading cell graph" << endl;
	connectivityGraph = readCellGraph(in_name);

	cout << " Building cell transects" << endl;
	cellPaths.clear();
	cellPaths = buildLinePaths(closestIndexToTransectPoint, connectivityGraph);

	cellSeedMask.clear();
	buildMaskFromFloodFill( seedPoints, cells, connectivityGraph, &cellSeedMask );

	cout << "Deleting cell center information" << endl;
	cells.clear();
	connectivityGraph.clear();

	cout << "Reading vertex locations" << endl;
	if(error = readVertices(in_name)){
		cout << " Error - " << error << endl;
		exit(error);
	}

	cout << "Marking vertices based on region definitions" << endl;
	vertexMasks = new int[vertices.size() * regionPolygons.size()];
	if(error = buildMasks(vertices, &vertexMasks[0])){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Marking points based on vertices" << endl;
	pointVertexIndices = new int[pointLocations.size()];
	if(error = buildPointIndices(pointLocations, vertices, &pointVertexIndices[0])){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << " Building closest vertex location for line segments" << endl;
	closestIndexToTransectPoint = buildClosestValuesForTransects(transectPoints, vertices);

	cout << " Reading vertex graph" << endl;
	connectivityGraph = readVertexGraph(in_name);

	cout << " Building vertex transects" << endl;
	vertexPaths.clear();
	vertexPaths = buildLinePaths(closestIndexToTransectPoint, connectivityGraph);

	cout << "Deleting vertex information" << endl;
	vertices.clear();
	connectivityGraph.clear();

	cout << "Reading edge locations" << endl;
	if(error = readEdges(in_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << " Building closest edge location for line segments" << endl;
	closestIndexToTransectPoint = buildClosestValuesForTransects(transectPoints, edges);

	cout << " Reading edge graph" << endl;
	connectivityGraph = readEdgeGraph(in_name);

	cout << " Building edge transects" << endl;
	edgePaths.clear();
	edgePaths = buildLinePaths(closestIndexToTransectPoint, connectivityGraph);

	cout << " Building edge transect path signs" << endl;
	edgePathSigns.clear();
	edgePathSigns = buildEdgePathSigns( in_name, edgePaths );

	cout << "Deleting edge information" << endl;
	edges.clear();
	connectivityGraph.clear();

	cout << "Writing mask dimensions" << endl;
	if(error = outputMaskDimensions(out_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}

	cout << "Writing mask attributes" << endl;
	if(error = outputMaskAttributes(out_name, in_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}
	cout << "Writing mask fields" << endl;
	if(error = outputMaskFields(out_name)){
		cout << "Error - " << error << endl;
		exit(error);
	}
}

/* Utility functions {{{*/
int featureIndex( const string featureName, const vector<string> featureNames ) {/*{{{*/
	vector<string>::const_iterator name_itr;
	int idx;
	bool found;

	for ( name_itr = featureNames.begin(), idx = 0; name_itr != featureNames.end(); name_itr++, idx++ ){
		if ( (*name_itr) == featureName ) {
			return idx;
		}
	}

	return -1;
}/*}}}*/

int stringType ( const string testStr ){/*{{{*/
	if ( testStr.find_first_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ") != std::string::npos ) {
		return text;
	} else if ( testStr.find_first_of(".") != std::string::npos ) {
		return num_double;
	}
	return num_int;
}/*}}}*/

vector< tuple<string, string, int> > extractFeatureProperties( Json::Value properties ){/*{{{*/
	Json::Value::Members keys = properties.getMemberNames();
	size_t numKeys = keys.size();
	vector< tuple<string, string, int> > prop_vec;
	vector<string>::iterator prop_itr;
	tuple<string, string, int> prop_val;

	prop_vec.clear();

	for ( prop_itr = keys.begin(); prop_itr != keys.end(); prop_itr++ ) {
		// We ignore name, author, tags, component, and object properties, since those are already dealt with.
		if ( (*prop_itr) != "name" && (*prop_itr) != "author" && (*prop_itr) != "tags" && (*prop_itr) != "component" && (*prop_itr) != "object" ) {
			get<0>(prop_val) = (*prop_itr);
			get<1>(prop_val) = properties[ (*prop_itr) ].asString();
			get<2>(prop_val) = stringType( properties[(*prop_itr)].asString() );

			prop_vec.push_back( prop_val );
		}
	}

	return prop_vec;
}/*}}}*/

vector<int> walkGraph( const vector< vector< pair<int, double> > > graph, const int start, const int end ){/*{{{*/
	vector< pair<int, double> > toCheck; // Stack of nodes to check first is node index, second is node dist. Should be sorted so last element has the smallest dist.
	vector<int> visited;
	vector<int> arrivedFrom;
	vector<int> path;

	pair<int, double> weightedNode;

	int numChecked;
	int i, currNode, testNode;
	int sortIdx, swpIdx;
	bool foundNode;

	double min_dist, swpDist;
	double edgeWeight = 1.0;

	toCheck.clear();
	visited.clear();
	arrivedFrom.clear();

	toCheck.reserve( graph.size() / 2 );

	for ( i = 0; i < graph.size(); i++ ) {
		visited.push_back( 0 );
		arrivedFrom.push_back( -1 );
	}

	weightedNode.first = start;
	weightedNode.second = 0;

	toCheck.push_back( weightedNode );
	numChecked = 0;

	while ( toCheck.size() > 0 ) {
		currNode = toCheck.back().first;
		min_dist = toCheck.back().second;
		toCheck.pop_back();

		if ( visited[currNode] ) {
			cout << " I've already visited this node? Weird! " << currNode << endl;
			continue;
		}

		visited[currNode] = 1;

		for ( i = 0; i < graph[currNode].size(); i++ ) {
			testNode = graph[currNode][i].first;
			edgeWeight = graph[currNode][i].second;

			// Prevent cycles by only checking nodes that haven't been visited already.
			if ( !visited[testNode] ) {

				// See if node is in stack already
				foundNode = false;
				for ( int j = 0; j < toCheck.size() && !foundNode; j++ ){
					if ( toCheck[j].first == testNode ) {
						foundNode = true;

						// Check the distance, and re-sort
						// Only need to sort forward, as we would only make the distance smaller.
						if ( toCheck[j].second > min_dist + edgeWeight ) {
							toCheck[j].second = min_dist + edgeWeight;
							arrivedFrom[ testNode ] = currNode;

							if ( j + 1 < toCheck.size() ) {
								bool sorted;
								sorted = false;
								sortIdx = j;

								// Re-sort, until we find a neighboring index that has a smaller dist
								for ( int k = j + 1; k < toCheck.size() && !sorted; k++ ) {
									if ( toCheck[sortIdx].second < toCheck[k].second ) {
										//cout << " Swapping1: " << toCheck[sortIdx].first << " with " << toCheck[k].first << endl;
										swpIdx = toCheck[sortIdx].first;
										swpDist = toCheck[sortIdx].second;

										toCheck[sortIdx].first = toCheck[k].first;
										toCheck[sortIdx].second = toCheck[k].second;

										toCheck[k].first = swpIdx;
										toCheck[k].second = swpDist;

										sortIdx = k;
									} else {
										sorted = true;
									}
								}
							}
						}
					}
				}

				// If the node wasn't found, insert it at the end, and sort it backwards into place.
				if ( !foundNode ) {
					// Build and insert
					weightedNode.first = testNode;
					weightedNode.second = min_dist + edgeWeight;

					arrivedFrom[ testNode ] = currNode;

					toCheck.push_back( weightedNode );

					// Sort, since it's added at the end of the list, we need to swap it with neighbors that are smaller
					// Until all nodes to the left are larger than it.
					bool sorted = false;
					for ( int j = toCheck.size()-1; j > 0 && !sorted; j-- ) {
						if ( toCheck[j].second > toCheck[j-1].second ) {
							//cout << " Swapping2: " << toCheck[j].first << " with " << toCheck[j-1].first << endl;
							swpIdx = toCheck[j].first;
							swpDist = toCheck[j].second;

							toCheck[j].first = toCheck[j-1].first;
							toCheck[j].second = toCheck[j-1].second;

							toCheck[j-1].first = swpIdx;
							toCheck[j-1].second = swpDist;
						} else {
							// Assume the stack was already sorted, so if we find a neighbor that has a larger distance
							// we can stop storting.
							sorted = true;
						}
					}
				}
			}
		}

		if ( currNode == end ) {
			toCheck.clear();
		}
		numChecked++;
	}

	path.clear();

	if ( !visited[end] || arrivedFrom[end] == -1 ) {
		arrivedFrom.clear();
		visited.clear();
		return path;
	}

	vector<int> r_path;

	currNode = end;

	while ( currNode != -1 ) {
		r_path.push_back ( currNode );
		currNode = arrivedFrom[currNode];
	}

	vector<int>::reverse_iterator rint_itr;

	for ( rint_itr = r_path.rbegin(); rint_itr != r_path.rend(); rint_itr++ ) {
		path.push_back( (*rint_itr) );
	}

	r_path.clear();
	arrivedFrom.clear();
	visited.clear();

	return path;
}/*}}}*/

void buildMaskFromFloodFill( const vector<pnt> seedLocations, const vector<pnt> staticLocations, const vector< vector< pair<int, double> > > graph, vector<int> *mask ){/*{{{*/
	int i;
	int numLocs = staticLocations.size();
	double dist, min_dist;
	int min_loc, idx;

	vector<pnt>::const_iterator seed_itr, stat_itr;
	vector<int> toCheck;

	toCheck.clear();
	(*mask).clear();

	// Init mask with zeroes
	for ( i = 0; i < staticLocations.size(); i++ ){
		(*mask).push_back(0);
	}

	// Find seed locations, and add them to the toCheck vector
	for ( seed_itr = seedLocations.begin(); seed_itr != seedLocations.end(); seed_itr++ ) {
		min_dist = INFINITY;

		for ( stat_itr = staticLocations.begin(), idx = 0; stat_itr != staticLocations.end(); stat_itr++, idx++ ) {
			dist = (*seed_itr).dotForAngle( (*stat_itr) );

			if ( min_dist > dist ) {
				min_dist = dist;
				min_loc = idx;
			}
		}

		toCheck.push_back(min_loc);
	}

	// Mark all of the seed points as masked
	for ( vector<int>::iterator idx_itr = toCheck.begin(); idx_itr != toCheck.end(); idx_itr++ ) {
		(*mask)[ (*idx_itr) ] = 1;
	}

	// Start the flood fill. For a given location, add it's neighbors (and mask them) if they haven't been masked already.
	while ( toCheck.size() > 0 ) {
		int currIdx = toCheck.back();
		toCheck.pop_back();

		for ( i = 0; i < graph[currIdx].size(); i++ ) {
			idx = graph.at(currIdx).at(i).first;
			if ( !(*mask)[ idx ] ) {
				toCheck.push_back( idx );
				(*mask)[ idx ] = 1;
			}
		}
	}
}/*}}}*/
/*}}}*/

/* Building/Ordering functions {{{ */
int readGridInfo(const string inputFilename){/*{{{*/
#ifdef _DEBUG
	cout << endl << endl << "Begin function: readGridInput" << endl << endl;
#endif

	nCells = netcdf_mpas_read_dim(inputFilename, "nCells");
	nVertices = netcdf_mpas_read_dim(inputFilename, "nVertices");
	nEdges = netcdf_mpas_read_dim(inputFilename, "nEdges");
	maxEdges = netcdf_mpas_read_dim(inputFilename, "maxEdges");
	vertexDegree = netcdf_mpas_read_dim(inputFilename, "vertexDegree");
#ifdef _DEBUG
	cout << "   Reading on_a_sphere" << endl;
#endif
	spherical = netcdf_mpas_read_onsphere(inputFilename);
#ifdef _DEBUG
	cout << "   Reading sphere_radius" << endl;
#endif
	sphereRadius = netcdf_mpas_read_sphereradius(inputFilename);
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

	cout << "Read dimensions:" << endl;
	cout << "    nCells = " << nCells << endl;
	cout << "    nVertices = " << nVertices << endl;
	cout << "    nEdges = " << nEdges << endl;
	cout << "    maxEdges = " << maxEdges << endl;
	cout << "    vertexDegree = " << vertexDegree << endl;
	cout << "    Spherical? = " << spherical << endl;

	if ( ! spherical ) {
		cout << endl << " -- WARNING: This tool only works correctly with real latitude / longitude values for cell centers." << endl;
		cout << "             If these are incorrect, you'll likely see odd behavior." << endl;
	}
	return 0;
}/*}}}*/
int readCells(const string inputFilename){/*{{{*/
	double *latcell, *loncell;
	pnt new_location;

#ifdef _DEBUG
	cout << endl << endl << "Begin function: readCells" << endl << endl;
#endif

	// Build cell center location information
	latcell = new double[nCells];
	loncell = new double[nCells];

	netcdf_mpas_read_latloncell ( inputFilename, nCells, latcell, loncell );

	cells.clear();
	for(int i = 0; i < nCells; i++){
		new_location = pntFromLatLon(latcell[i], loncell[i]);
		new_location.setPositiveLonRange( lonRangePositive );
		new_location.idx = i;

		if(spherical) new_location.normalize();
		cells.push_back(new_location);
	}

	cout << "Built " << cells.size() << " cells." << endl;
	delete[] latcell;
	delete[] loncell;

	return 0;
}/*}}}*/
int readVertices(const string inputFilename){/*{{{*/
	double *latvertex, *lonvertex;
	pnt new_location;

#ifdef _DEBUG
	cout << endl << endl << "Begin function: readVertices" << endl << endl;
#endif

	// Build cell center location information
	latvertex = new double[nVertices];
	lonvertex = new double[nVertices];

	netcdf_mpas_read_latlonvertex ( inputFilename, nVertices, latvertex, lonvertex );

	vertices.clear();
	for(int i = 0; i < nVertices; i++){
		new_location = pntFromLatLon(latvertex[i], lonvertex[i]);
		new_location.idx = i;

		if(spherical) new_location.normalize();
		vertices.push_back(new_location);
	}

	cout << "Built " << vertices.size() << " vertices." << endl;
	delete[] latvertex;
	delete[] lonvertex;

	return 0;
}/*}}}*/
int readEdges(const string inputFilename){/*{{{*/
	double *latedge, *lonedge;
	pnt new_location;

#ifdef _DEBUG
	cout << endl << endl << "Begin function: readEdges" << endl << endl;
#endif

	// Build cell center location information
	latedge = new double[nEdges];
	lonedge = new double[nEdges];

	netcdf_mpas_read_latlonedge ( inputFilename, nEdges, latedge, lonedge );

	edges.clear();
	for(int i = 0; i < nEdges; i++){
		new_location = pntFromLatLon(latedge[i], lonedge[i]);
		new_location.idx = i;

		if(spherical) new_location.normalize();
		edges.push_back(new_location);
	}

	cout << "Built " << edges.size() << " edges." << endl;
	delete[] latedge;
	delete[] lonedge;

	return 0;
}/*}}}*/
vector< vector< pair<int, double> > > readCellGraph(const string inputFilename){/*{{{*/
	vector< vector< pair<int, double> > > graph;
	int *cellsoncell, *edgesoncell;
	int *nedgesoncell;
	double *dcedge;
	int cellid, edgeid;
	pair<int, double> connection;
	vector< pair<int, double> > cellList;

	graph.clear();
	nedgesoncell = new int[nCells];
	cellsoncell = new int[nCells * maxEdges];
	edgesoncell = new int[nCells * maxEdges];
	dcedge = new double[nEdges];
	netcdf_mpas_read_nedgesoncell ( inputFilename, nCells, nedgesoncell );
	netcdf_mpas_read_edgesoncell ( inputFilename, nCells, maxEdges, edgesoncell );
	netcdf_mpas_read_cellsoncell ( inputFilename, nCells, maxEdges, cellsoncell );
	netcdf_mpas_read_dcedge ( inputFilename, nEdges, dcedge );

	for ( int i = 0; i < nCells; i++ ) {
		cellList.clear();
		for ( int j = 0; j < nedgesoncell[i]; j++ ) {
			cellid = cellsoncell[ i * maxEdges + j ] - 1;
			edgeid = edgesoncell[ i * maxEdges + j ] - 1;

			if ( cellid >= 0 ) {
				connection.first = cellid;
				connection.second = dcedge[edgeid];
				cellList.push_back( connection );
			}
		}
		graph.push_back(cellList);
	}

	delete[] nedgesoncell;
	delete[] cellsoncell;
	delete[] edgesoncell;
	delete[] dcedge;
	cellList.clear();

	return graph;
}/*}}}*/
vector< vector< pair<int, double> > > readVertexGraph(const string inputFilename){/*{{{*/
	vector< vector< pair<int, double> > > graph;
	int *verticesonedge, *edgesonvertex;
	double *dvedge;
	bool add_vertid1, add_vertid2;
	int vertid1, vertid2, edgeid;
	pair<int, double> connection;
	vector< pair<int, double> > vertexList;

	graph.clear();
	verticesonedge = new int[nEdges * 2];
	edgesonvertex = new int[nVertices * vertexDegree];
	dvedge = new double[nEdges];
	netcdf_mpas_read_verticesonedge ( inputFilename, nEdges, verticesonedge );
	netcdf_mpas_read_edgesonvertex ( inputFilename, nVertices, vertexDegree, edgesonvertex );
	netcdf_mpas_read_dvedge ( inputFilename, nEdges, dvedge );

	for ( int i = 0; i < nVertices; i++ ) {
		vertexList.clear();

		for ( int j = 0; j < vertexDegree; j++ ) {
			edgeid = edgesonvertex[ i * vertexDegree + j ] - 1;

			if ( edgeid != -1 ) {
				vertid1 = verticesonedge[ edgeid * 2 ] - 1;
				vertid2 = verticesonedge[ edgeid * 2 + 1 ] - 1;

				add_vertid1 = true;
				add_vertid2 = true;
				if ( vertid1 == i ) add_vertid1 = false;
				if ( vertid2 == i ) add_vertid2 = false;

				if ( vertid1 >= nVertices ) add_vertid1 = false;
				if ( vertid2 >= nVertices ) add_vertid1 = false;

				for ( int k = 0; k < vertexList.size() && add_vertid1 && add_vertid2; k++ ) {
					if ( vertid1 == vertexList[k].first ) add_vertid1 = false;
					if ( vertid2 == vertexList[k].first ) add_vertid2 = false;
				}

				if ( add_vertid1 ) {
					connection.first = vertid1;
					connection.second = dvedge[ edgeid ];
					vertexList.push_back( connection );
				}

				if ( add_vertid2 ) {
					connection.first = vertid2;
					connection.second = dvedge[ edgeid ];
					vertexList.push_back( connection );
				}
			}
		}
		graph.push_back( vertexList );
	}

	delete[] edgesonvertex;
	delete[] verticesonedge;
	delete[] dvedge;
	vertexList.clear();

	return graph;
}/*}}}*/
vector< vector< pair<int, double> > > readEdgeGraph(const string inputFilename){/*{{{*/
	vector< vector< pair<int, double> > > graph;
	int *verticesonedge, *edgesonvertex;
	double *dvedge;
	int vertid, edgeid;
	bool add_edge;
	pair<int, double> connection;
	vector< pair<int, double> > edgeList;

	graph.clear();
	verticesonedge = new int[nEdges * 2];
	edgesonvertex = new int[nVertices * vertexDegree];
	dvedge = new double[nEdges];
	netcdf_mpas_read_verticesonedge ( inputFilename, nEdges, verticesonedge );
	netcdf_mpas_read_edgesonvertex ( inputFilename, nVertices, vertexDegree, edgesonvertex );
	netcdf_mpas_read_dvedge ( inputFilename, nEdges, dvedge );

	for ( int i = 0; i < nEdges; i++ ) {
		edgeList.clear();

		for ( int j = 0; j < 2; j++ ) {
			vertid = verticesonedge[ i * 2 + j ] - 1;

			for ( int k = 0; k < vertexDegree; k++ ) {
				edgeid = edgesonvertex[ vertid * vertexDegree + k ] - 1;

				if ( edgeid != -1 ) {
					add_edge = true;

					for ( int l = 0; l < edgeList.size(); l++ ) {
						if ( edgeid == edgeList[l].first ) add_edge = false;
					}

					if ( add_edge ) {
						connection.first = edgeid;
						connection.second = dvedge[i];
						edgeList.push_back( connection );
					}
				}
			}
		}
		graph.push_back( edgeList );
	}

	delete[] edgesonvertex;
	delete[] verticesonedge;
	delete[] dvedge;
	edgeList.clear();

	return graph;
}/*}}}*/
int resetFeatureInfo(){/*{{{*/
	connectivityGraph.clear();

	seedPoints.clear();

	regionPolygons.clear();
	regionNames.clear();
	regionsInGroup.clear();
	regionGroupNames.clear();

	transectPoints.clear();
	transectNames.clear();
	transectsInGroup.clear();
	transectGroupNames.clear();

	polygonConstants.clear();
	polygonMultiples.clear();
	pointLocations.clear();

	return 0;
}/*}}}*/
int getFeatureInfo(const string featureFilename){/*{{{*/
	ifstream json_file(featureFilename.c_str());
	Json::Value root;
	string groupName, tempGroupName;
	vector<int> pointIndices;
	vector<int> regionIndices;
	vector<int> transectIndices;
	vector<pnt> regionVertices;
	vector<pnt> transectVertices;
	vector< vector<pnt> > polygonList;
	vector< vector<pnt> > lineSegList;
	vector< tuple<string, string, int> > properties;
	bool addGroupToRegList = true;
	bool addGroupToTrnList = true;
	bool addGroupToPntList = true;

	json_file >> root;

	json_file.close();

	groupName = root["groupName"].asString();

	if ( groupName == "enterNameHere" ) {
		cout << " ** WARNING: Features file at " << featureFilename << " has an unset group name." << endl;
		cout << "             If you want to use groups to control sets of regions, define the group name, and rerun." << endl;
		cout << "             Using a deafult unique group name of '[feature]Group[N]' where [feature] is the feature type {region, point}" << endl;
		cout << "             and [N] is the current number of groups for that feature type." << endl;
	}
#ifdef _DEBUG
	cout << " Starting group: " << groupName << endl;
#endif

	pointIndices.clear();
	transectIndices.clear();
	regionIndices.clear();

	for ( int i = 0; i < root["features"].size(); i++ ){
		Json::Value feature = root["features"][i];
		string featureName = feature["properties"]["name"].asString();

		// Handle regions
		if ( feature["properties"]["object"].asString() == "region" ) {
			bool addRegToGroup;
			bool add_region;
			int regionIdx, idx;

			add_region = true;
			regionIdx = featureIndex( featureName, regionNames );
			if ( regionIdx >= 0 ) {
				add_region = false;
			}

			if ( add_region ) {
				polygonList.clear();
				regionIdx = regionNames.size();
#ifdef _DEBUG
				cout << "Adding region: " << featureName << " with index " << regionNames.size() << endl;
#endif
				regionNames.push_back(featureName);

				Json::Value geometry = feature["geometry"];
				properties.clear();
				properties = extractFeatureProperties( feature["properties"] );
				regionProperties.push_back(properties);
				Json::Value coordinates = geometry["coordinates"];

				if ( geometry["type"].asString() == "Polygon" ) {

					regionVertices.clear();
					for (int i = 0; i < coordinates[0].size(); i++){
						double lon = coordinates[0][i][0].asDouble() * M_PI/180.0;
						double lat = coordinates[0][i][1].asDouble() * M_PI/180.0;

#ifdef _DEBUG
						cout << " Added a region vertex: " << lat << ", " << lon << endl;
#endif
						pnt point = pntFromLatLon(lat, lon);
						point.setPositiveLonRange( lonRangePositive );
						regionVertices.push_back(point);
					}
#ifdef _DEBUG
					cout << "Added a polygon with: " << coordinates.size() << " vertices" << endl;
#endif
					polygonList.push_back(regionVertices);

				} else if (geometry["type"].asString() == "MultiPolygon" ) {
					for (int i = 0; i < coordinates.size(); i++){
						regionVertices.clear();
						for (int j = 0; j < coordinates[i][0].size(); j++){
							double lon = coordinates[i][0][j][0].asDouble() * M_PI/180.0;
							double lat = coordinates[i][0][j][1].asDouble() * M_PI/180.0;

#ifdef _DEBUG
							cout << " Added a region vertex: " << lat << ", " << lon << endl;
#endif
							pnt point = pntFromLatLon(lat, lon);
							point.setPositiveLonRange( lonRangePositive );
							regionVertices.push_back(point);
						}
#ifdef _DEBUG
						cout << "Added a polygon with: " << coordinates[i].size() << " vertices" << endl;
#endif
						polygonList.push_back(regionVertices);
						maxRegVertices = max(maxRegVertices, (int)regionVertices.size());
					}
#ifdef _DEBUG
					cout << "Added a total of: " << coordinates.size() << " polygons" << endl;
#endif
				}

				regionPolygons.push_back(polygonList);
				maxRegPolygons = max(maxRegPolygons, (int)polygonList.size());
				nRegions = nRegions + 1;

				// ** Check for region in group already
				addRegToGroup = true;
				for ( int_itr = regionIndices.begin(); int_itr != regionIndices.end(); int_itr++){
					if ( regionIdx == (*int_itr) ) {
						addRegToGroup = false;
					}
				}
				if ( addRegToGroup ) {
					regionIndices.push_back(regionIdx);
				}

				// * Add group name to list of groups for regions, and prevent it from being added with subsequent regions
				if ( addGroupToRegList ) {
					addGroupToRegList = false;
					if ( groupName == "enterNameHere" ){
						ostringstream convert;
						convert << regionGroupNames.size() + 1;
						tempGroupName = "regionGroup" + convert.str();
						regionGroupNames.push_back(tempGroupName);
					} else {
						regionGroupNames.push_back(groupName);
					}
				}

			}

		// Handle transects
		} else if ( feature["properties"]["object"].asString() == "transect" ) {
			bool addTransectToGroup;
			bool add_transect;
			int transectIdx, idx;

			add_transect = true;
			transectIdx = featureIndex( featureName, transectNames );
			if ( transectIdx != -1 ) {
				add_transect = false;
			}

			if ( add_transect ) {
				lineSegList.clear();
				transectIdx = transectNames.size();

#ifdef _DEBUG
				cout << "Adding transect: " << featureName << " with index " << transectIdx << endl;
#endif

				transectNames.push_back(featureName);

				Json::Value geometry = feature["geometry"];
				properties.clear();
				properties = extractFeatureProperties( feature["properties"] );
				transectProperties.push_back(properties);
				Json::Value coordinates = geometry["coordinates"];

				if ( geometry["type"].asString() == "LineString" ) {
					transectVertices.clear();

					for ( int i = 0; i < coordinates.size(); i++ ) {
						double lon = coordinates[i][0].asDouble() * M_PI/180.0;
						double lat = coordinates[i][1].asDouble() * M_PI/180.0;

#ifdef _DEBUG
						cout << " Added a transect vertex: " << lat << ", " << lon << endl;
#endif

						pnt point = pntFromLatLon(lat, lon);
						point.setPositiveLonRange( lonRangePositive );
						transectVertices.push_back(point);
					}

#ifdef _DEBUG
					cout << "Added a line segment with: " << coordinates.size() << " vertices" << endl;
#endif
					lineSegList.push_back(transectVertices);
					maxTrnVertices = max(maxTrnVertices, (int)transectVertices.size());
				} else if ( geometry["type"].asString() == "MultiLineString" ) {
					for ( int i = 0; i < coordinates.size(); i++ ){
						transectVertices.clear();
						for ( int j = 0; j < coordinates[i].size(); j++ ){
							double lon = coordinates[i][j][0].asDouble() * M_PI/180.0;
							double lat = coordinates[i][j][1].asDouble() * M_PI/180.0;

#ifdef _DEBUG
							cout << " Added a transect vertex: " << lat << ", " << lon << endl;
#endif
							pnt point = pntFromLatLon(lat, lon);
							point.setPositiveLonRange( lonRangePositive );
							transectVertices.push_back(point);
						}

#ifdef _DEBUG
						cout << "Added a transect with: " << coordinates[i].size() << " vertices" << endl;
#endif
						lineSegList.push_back(transectVertices);
						maxTrnVertices = max( maxTrnVertices, (int)transectVertices.size() );
					}
#ifdef _DEBUG
					cout << "Added a total of: " << coordinates.size() << " line segments" << endl;
#endif
				}

				transectPoints.push_back(lineSegList);
				maxTrnLines = max(maxTrnLines, (int)lineSegList.size());
				nTransects = nTransects + 1;

				addTransectToGroup = true;
				for ( int_itr = transectIndices.begin(); int_itr != transectIndices.end(); int_itr++ ) {
					if ( transectIdx == (*int_itr) ) {
						addTransectToGroup = false;
					}
				}

				if ( addTransectToGroup ) {
					transectIndices.push_back(transectIdx);
				}

				if ( addGroupToTrnList ) {
					addGroupToTrnList = false;
					if ( groupName == "enterNameHere" ){
						ostringstream convert;
						convert << transectGroupNames.size() + 1;
						tempGroupName = "transectGroup" + convert.str();
						transectGroupNames.push_back(tempGroupName);
					} else {
						transectGroupNames.push_back(groupName);
					}
				}
			}
		// Handle points
		} else if ( feature["properties"]["object"].asString() == "point" ) {
			bool addPointToGroup;
			bool add_point;
			int pointIdx, idx;

			add_point = true;
			pointIdx = featureIndex( featureName, pointNames );
			if ( pointIdx != -1 ) {
				add_point = false;
			}

			if ( add_point ) {
				double lon = feature["geometry"]["coordinates"][0].asDouble() * M_PI/180.0 ;
				double lat = feature["geometry"]["coordinates"][1].asDouble() * M_PI/180.0;

				pointIdx = pointNames.size();

#ifdef _DEBUG
				cout << "Adding point: " << featureName << " with index " << pointNames.size() << endl;
#endif

				pnt point = pntFromLatLon(lat, lon);
				point.setPositiveLonRange( lonRangePositive );
				point.idx = pointIdx;

				properties.clear();
				properties = extractFeatureProperties( feature["properties"] );
				pointProperties.push_back(properties);

				pointLocations.push_back(point);
				pointNames.push_back(featureName);
				nPoints = nPoints + 1;

				// ** Check for point in group already
				addPointToGroup = true;
				for ( int_itr = pointIndices.begin(); int_itr != pointIndices.end(); int_itr++){
					if ( pointIdx == (*int_itr) ) {
						addPointToGroup = false;
					}
				}

				if ( addPointToGroup ) {
					pointIndices.push_back(pointIdx);
				}

				// * Add group name to list of groups for points, and prevent it from being added with subsequent points
				if ( addGroupToPntList ) {
					addGroupToPntList = false;
					if ( groupName == "enterNameHere" ){
						ostringstream convert;
						convert << pointGroupNames.size() + 1;
						tempGroupName = "pointGroup" + convert.str();
						pointGroupNames.push_back(tempGroupName);
					} else {
						pointGroupNames.push_back(groupName);
					}
				}

			}
		} else {
			cout << "Skipping feature, can only process regions at this time..." << endl;
		}
	}


	if ( (int)regionIndices.size() > 0 ) {
		maxRegionsInGroup = max(maxRegionsInGroup, (int)regionIndices.size());
		regionsInGroup.push_back(regionIndices);
	}


	if ( (int)transectIndices.size() > 0 ) {
		maxTransectsInGroup = max(maxTransectsInGroup, (int)transectIndices.size());
		transectsInGroup.push_back(transectIndices);
	}

	if ( (int)pointIndices.size() > 0 ) {
		maxPointsInGroup = max(maxPointsInGroup, (int)pointIndices.size());
		pointsInGroup.push_back(pointIndices);
	}

	return 0;
}/*}}}*/
int getSeedInfo(const string seedFilename){/*{{{*/
	ifstream json_file(seedFilename.c_str());
	Json::Value root;

	json_file >> root;

	json_file.close();

	for ( int i = 0; i < root["features"].size(); i++ ){
		Json::Value feature = root["features"][i];

		if ( feature["properties"]["object"].asString() == "point" ) {
			double lon = feature["geometry"]["coordinates"][0].asDouble() * M_PI/180.0 ;
			double lat = feature["geometry"]["coordinates"][1].asDouble() * M_PI/180.0;

			pnt point = pntFromLatLon(lat, lon);
			point.setPositiveLonRange( lonRangePositive );
			point.idx = seedPoints.size();

#ifdef _DEBUG
			cout << "Adding seed point with index " << seedPoints.size() << endl;
#endif
			seedPoints.push_back(point);
		}
	}

	return 0;
}/*}}}*/
vector< vector< vector<int> > > buildClosestValuesForTransects(const vector< vector< vector<pnt> > > lineSegments, const vector<pnt> staticLocations){/*{{{*/
	vector< vector< vector<pnt> > >::const_iterator cseg_itr;
	vector< vector<pnt> >::const_iterator cline_itr;
	vector<pnt>::const_iterator cpnt_itr;
	vector< vector< vector<int> > > closestIndices;

	vector<int> minLocPoints;
	vector< vector<int> > minLocLines;

	double dist, min_dist;
	int min_loc, idx;;

	closestIndices.clear();

	for ( cseg_itr = lineSegments.begin(); cseg_itr != lineSegments.end(); cseg_itr++ ) {
		minLocLines.clear();

		for ( cline_itr = (*cseg_itr).begin(); cline_itr != (*cseg_itr).end(); cline_itr++ ) {
			minLocPoints.clear();
			for ( cpnt_itr = (*cline_itr).begin(); cpnt_itr != (*cline_itr).end(); cpnt_itr++ ){
				min_dist = 4.0 * M_PI;
				for ( idx = 0; idx < staticLocations.size(); idx++ ){
					dist = (*cpnt_itr).dotForAngle( staticLocations[idx] );

					if ( dist < min_dist ) {
						min_dist = dist;
						min_loc = idx;
					}
				}

				minLocPoints.push_back(min_loc);
			}

			minLocLines.push_back(minLocPoints);
		}
		closestIndices.push_back(minLocLines);
	}

	return closestIndices;
}/*}}}*/
vector< vector<int> > buildLinePaths( const vector< vector< vector<int> > > closestIndices, const vector< vector< pair<int, double> > > graph){/*{{{*/
	vector< vector< vector<int> > >::const_iterator cseg_itr;
	vector< vector<int> >::const_iterator cline_itr;
	vector< vector<int> > pathIndices ;

	vector<int> full_path, partial_path;

	int i, j;
	int start, end;

	partial_path.clear();
	pathIndices.clear();

	for ( cseg_itr = closestIndices.begin(); cseg_itr != closestIndices.end(); cseg_itr++ ) {
		full_path.clear();
		for ( cline_itr = (*cseg_itr).begin(); cline_itr != (*cseg_itr).end(); cline_itr++ ) {
			full_path.push_back( (*cline_itr).at(0) );
			for ( i = 0; i < (*cline_itr).size() - 1; i++ ){
				start = (*cline_itr).at( i );
				end = (*cline_itr).at( i + 1 );
				partial_path = walkGraph( graph, start, end );

				if ( partial_path.size() > 0 && ( full_path[ full_path.size() - 1 ] != partial_path[0] ) ) {
					full_path.push_back( partial_path[0] );
				}

				for ( j = 1; j < partial_path.size(); j++ ) {
					full_path.push_back( partial_path[j] );
				}

				partial_path.clear();
			}
		}
		pathIndices.push_back(full_path);
	}

	return pathIndices;
}/*}}}*/
vector< vector<int> > buildEdgePathSigns( const string inputFilename, const vector< vector<int> > edgePaths ){/*{{{*/
	int currEdge, nextEdge;
	int sharedVert, edgeSign;

	vector< vector<int> >::const_iterator paths_itr;
	vector< vector<int> > edgeSigns;
	vector<int> pathSigns;

	int *verticesonedge;

	pathSigns.clear();
	verticesonedge = new int [nEdges * 2];

	netcdf_mpas_read_verticesonedge ( inputFilename, nEdges, verticesonedge );

#ifdef _DEBUG
	cout << " Building signs for: " << edgePaths.size() << " paths." << endl;
#endif

	for ( paths_itr = edgePaths.begin(); paths_itr != edgePaths.end(); paths_itr++ ) {
		pathSigns.clear();
#ifdef _DEBUG
		cout << "    -- This path has: " << (*paths_itr).size() << " edges. " << endl;
#endif

		for ( int i = 0; i < (*paths_itr).size()-1; i++ ) {
			currEdge = (*paths_itr).at(i);
			nextEdge = (*paths_itr).at(i+1);

			sharedVert = -1;

			if ( verticesonedge[currEdge*2] == verticesonedge[nextEdge*2] || verticesonedge[currEdge*2] == verticesonedge[nextEdge*2 + 1] ){
				sharedVert = 1;
			}  else if ( verticesonedge[currEdge*2 + 1] == verticesonedge[nextEdge*2] || verticesonedge[currEdge*2 + 1] == verticesonedge[nextEdge*2 + 1] ){
				sharedVert = 2;
			}

			if ( sharedVert == 1 ) {
				edgeSign = -1;
			} else if ( sharedVert == 2 ) {
				edgeSign = 1;
			} else {
				cout << " ERROR: Edges " << currEdge << " and " << nextEdge << " don't share any vertices. Exiting..." << endl;
				exit(1);
			}

			pathSigns.push_back(edgeSign);
		}

		// Need to handle the last edge.
		if ( (*paths_itr).size() > 1 ) {
			currEdge = (*paths_itr).at( (*paths_itr).size() - 1 );
			nextEdge = (*paths_itr).at( (*paths_itr).size() - 2 ); // Really, edge before currEdge

			sharedVert = -1;

			if ( verticesonedge[currEdge*2] == verticesonedge[nextEdge*2] || verticesonedge[currEdge*2] == verticesonedge[nextEdge*2 + 1] ){
				sharedVert = 1;
			}  else if ( verticesonedge[currEdge*2 + 1] == verticesonedge[nextEdge*2] || verticesonedge[currEdge*2 + 1] == verticesonedge[nextEdge*2 + 1] ){
				sharedVert = 2;
			}

			if ( sharedVert == 1 ) {
				edgeSign = 1;
			} else if ( sharedVert == 2 ) {
				edgeSign = -1;
			} else {
				cout << " ERROR: Edges " << currEdge << " and " << nextEdge << " don't share any vertices. Exiting..." << endl;
				exit(1);
			}

			pathSigns.push_back(edgeSign);

		} else if ( (*paths_itr).size() == 1 ) {
			// If there's only one edge, arbitrarily give it a sign of 1, since direction doesn't matter.
			edgeSign = 1;
			pathSigns.push_back(edgeSign);
		}

		edgeSigns.push_back(pathSigns);
	}

	delete[] verticesonedge;

	return edgeSigns;
}/*}}}*/
int buildPolygonValues(){/*{{{*/
	int iReg, iPoly, i, j;
	double vert1Lat, vert1Lon;
	double vert2Lat, vert2Lon;

	// Build the constants and multiples arrays for each polygon, as is done in markCells()
	// These define the line segment of each edge, so we can test if a ray intersects the edges.

	polygonConstants.resize(regionPolygons.size());
	polygonMultiples.resize(regionPolygons.size());

	for (reg_itr = regionPolygons.begin(), iReg = 0; reg_itr != regionPolygons.end(); reg_itr++, iReg++){
		polygonConstants.at(iReg).resize((*reg_itr).size());
		polygonMultiples.at(iReg).resize((*reg_itr).size());

		for (poly_itr = (*reg_itr).begin(), iPoly = 0; poly_itr != (*reg_itr).end(); poly_itr++, iPoly++){
			polygonConstants.at(iReg).at(iPoly).resize((*poly_itr).size());
			polygonMultiples.at(iReg).at(iPoly).resize((*poly_itr).size());

			for ( i = 0, j = (*poly_itr).size()-1; i < (*poly_itr).size(); j=i, i++){
				vert1Lat = (*poly_itr).at(j).getLat();
				vert1Lon = (*poly_itr).at(j).getLon();

				vert2Lat = (*poly_itr).at(i).getLat();
				vert2Lon = (*poly_itr).at(i).getLon();

				if ( vert1Lat == vert2Lat ) {
					polygonConstants.at(iReg).at(iPoly).at(i) = vert1Lat;
					polygonMultiples.at(iReg).at(iPoly).at(i) = 0;
				} else {
					polygonConstants.at(iReg).at(iPoly).at(i) = vert1Lon - (vert1Lat * vert2Lon) / (vert2Lat - vert1Lat) + (vert1Lat * vert1Lon) / (vert2Lat - vert1Lat);
					polygonMultiples.at(iReg).at(iPoly).at(i) = (vert2Lon - vert1Lon) / (vert2Lat - vert1Lat);
				}
			}
		}
	}

	return 0;
}/*}}}*/
int buildMasks(vector<pnt> locations, int *masks){/*{{{*/
	bool inReg;
	pnt vert1, vert2;
	pnt vec1, vec2;
	pnt crossProd;
	double dot;

	double locLat, locLon;
	double vert1Lat, vert1Lon;
	double vert2Lat, vert2Lon;
	bool oddSides;

	int iReg, iPoly;
	int iLoc, i, j, idx;

	if ( regionPolygons.size() <= 0 ) {
		return 0;
	}

	#pragma omp parallel for default(shared) private(locLat, locLon, iReg, reg_itr, inReg, poly_itr, iPoly, oddSides, i, j, vert1Lat, vert1Lon, vert2Lat, vert2Lon)
	for ( iLoc = 0; iLoc < locations.size(); iLoc++ ){
#ifdef _DEBUG
		cout << " Masking location: " << iLoc << " of " << locations.size() << endl;
#endif

		locLat = locations.at(iLoc).getLat();
		locLon = locations.at(iLoc).getLon();
		for (reg_itr = regionPolygons.begin(), iReg = 0; reg_itr != regionPolygons.end(); reg_itr++, iReg++){
			inReg = false;
			for (poly_itr = (*reg_itr).begin(), iPoly = 0; poly_itr != (*reg_itr).end() && !inReg; poly_itr++, iPoly++){

				// Test for the number of edge intersections with a zonal ray.
				// If the number is odd, the point is inside the polygon
				// If the number is even, the point is outside the polygon
				oddSides = false;
				for ( i = 0, j = (*poly_itr).size()-1; i < (*poly_itr).size(); j=i, i++){
					vert1Lat = (*poly_itr).at(j).getLat();
					vert1Lon = (*poly_itr).at(j).getLon();

					vert2Lat = (*poly_itr).at(i).getLat();
					vert2Lon = (*poly_itr).at(i).getLon();

					if ( (vert1Lat < locLat && vert2Lat >= locLat) || (vert2Lat < locLat && vert1Lat >= locLat) ) {
						oddSides = oddSides^(locLat * polygonMultiples.at(iReg).at(iPoly).at(i) + polygonConstants.at(iReg).at(iPoly).at(i) < locLon);
					}
				}
				if ( oddSides ) {
					inReg = true;
				}
			}
			if ( inReg ) {
				masks[iLoc * nRegions + iReg] = 1;
			} else {
				masks[iLoc * nRegions + iReg] = 0;
			}
		}
	} // */
	return 0;
}/*}}}*/
int buildPointIndices(vector<pnt> testLocations, vector<pnt> staticLocations, int *indices){/*{{{*/
	int iTestLoc, iStaticLoc, idx;
	double minDist, dist;

	#pragma omp parallel for default(shared) private(idx, minDist, iStaticLoc, dist)
	for ( iTestLoc = 0; iTestLoc < testLocations.size(); iTestLoc++){
		idx = -1;
		minDist = HUGE_VAL;
		for ( iStaticLoc = 0; iStaticLoc < staticLocations.size(); iStaticLoc++){
			dist = testLocations[iTestLoc].dotForAngle(staticLocations[iStaticLoc]);

			if ( dist < minDist ) {
#ifdef _DEBUG
				cout << "    Old dist / idx: " << minDist << " " << idx << endl;
				cout << "    New dist / idx: " << dist << " " << iStaticLoc << endl;
				cout << "    Static Lat/Lon: " << staticLocations[iStaticLoc].getLat() << " " << staticLocations[iStaticLoc].getLon() << endl;
				cout << "    Test Lat/Lon: " << testLocations[iTestLoc].getLat() << " " << testLocations[iTestLoc].getLon() << endl;
#endif
				minDist = dist;
				idx = iStaticLoc;
			}
		}

		if ( idx == -1 ) {
			cout << " ERROR: No locations found for test location " << iTestLoc << endl;
			cout << testLocations[iTestLoc] << endl;
			indices[iTestLoc] = 0;
		} else {
			indices[iTestLoc] = idx+1;
		}
	}

	return 0;
}/*}}}*/
int buildAllFeatureGroups(){/*{{{*/
	vector<int> groupIndices;
	int featureIdx;

	// 'all' region group
	groupIndices.clear();

	for ( str_itr = regionNames.begin(), featureIdx = 0; str_itr != regionNames.end(); str_itr++, featureIdx++){
		groupIndices.push_back(featureIdx);
	}

	regionsInGroup.push_back(groupIndices);
	regionGroupNames.push_back("all");
	maxRegionsInGroup = max(maxRegionsInGroup, (int)groupIndices.size());

	// 'all' transect group
	groupIndices.clear();

	for ( str_itr = transectNames.begin(), featureIdx = 0; str_itr != transectNames.end(); str_itr++, featureIdx++){
		groupIndices.push_back(featureIdx);
	}

	transectsInGroup.push_back(groupIndices);
	transectGroupNames.push_back("all");
	maxTransectsInGroup = max(maxTransectsInGroup, (int)groupIndices.size());

	// 'all' point group
	groupIndices.clear();

	for ( str_itr = pointNames.begin(), featureIdx = 0; str_itr != pointNames.end(); str_itr++, featureIdx++){
		groupIndices.push_back(featureIdx);
	}

	pointsInGroup.push_back(groupIndices);
	pointGroupNames.push_back("all");
	maxPointsInGroup = max(maxPointsInGroup, (int)groupIndices.size());

	return 0;
}/*}}}*/
/*}}}*/

/* Output functions {{{*/
int outputMaskDimensions( const string outputFilename ){/*{{{*/
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

	// open the mpas file
	NcFile grid(outputFilename.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);

	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;

	// define dimensions
	NcDim *tempDim;

	// write dimensions
	if (!(tempDim = grid.add_dim("nCells", nCells))) return NC_ERR;
	if (!(tempDim = grid.add_dim("nVertices", nVertices))) return NC_ERR;
	if (!(tempDim = grid.add_dim("nEdges", nEdges))) return NC_ERR;

	if ( nRegions > 0 ) {
		if (!(tempDim = grid.add_dim("nRegions", nRegions))) return NC_ERR;
		if (!(tempDim = grid.add_dim("nRegionGroups", regionGroupNames.size()))) return NC_ERR;
		if (!(tempDim = grid.add_dim("maxRegionsInGroup", maxRegionsInGroup))) return NC_ERR;
	}

	if ( nTransects > 0 ) {
		for ( int i = 0; i < cellPaths.size(); i++ ) {
			maxCellsInTransect = max( maxCellsInTransect, (int)cellPaths[i].size() );
		}
		for ( int i = 0; i < vertexPaths.size(); i++ ) {
			maxVerticesInTransect = max( maxVerticesInTransect, (int)vertexPaths[i].size() );
		}
		for ( int i = 0; i < edgePaths.size(); i++ ) {
			maxEdgesInTransect = max( maxEdgesInTransect, (int)edgePaths[i].size() );
		}

		if (!(tempDim = grid.add_dim("nTransects", nTransects))) return NC_ERR;
		if (!(tempDim = grid.add_dim("nTransectGroups", transectGroupNames.size()))) return NC_ERR;
		if (!(tempDim = grid.add_dim("maxTransectsInGroup", maxTransectsInGroup))) return NC_ERR;
		if (!(tempDim = grid.add_dim("maxCellsInTransect", maxCellsInTransect))) return NC_ERR;
		if (!(tempDim = grid.add_dim("maxVerticesInTransect", maxVerticesInTransect))) return NC_ERR;
		if (!(tempDim = grid.add_dim("maxEdgesInTransect", maxEdgesInTransect))) return NC_ERR;
	}

	if ( nPoints > 0 ) {
		if (!(tempDim = grid.add_dim("nPoints", nPoints))) return NC_ERR;
		if (!(tempDim = grid.add_dim("nPointGroups", pointGroupNames.size()))) return NC_ERR;
		if (!(tempDim = grid.add_dim("maxPointsInGroup", maxPointsInGroup))) return NC_ERR;
	}

	if (!(tempDim = grid.add_dim("StrLen", 1024))) return NC_ERR;

	grid.close();

	// file closed when file obj goes out of scope
	return 0;
}/*}}}*/
int outputMaskAttributes( const string outputFilename, const string inputFilename ){/*{{{*/
	/************************************************************************
	 *
	 * This function writes the grid dimensions to the netcdf file named
	 * outputFilename
	 *
	 * **********************************************************************/
	// Return this code to the OS in case of failure.
	static const int NC_ERR = 2;
	char mesh_spec_str[1024];

	// set error behaviour (matches fortran behaviour)
	NcError err(NcError::verbose_nonfatal);

	// open the mpas file
	NcFile grid(outputFilename.c_str(), NcFile::Write);

	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;
	NcBool sphereAtt, radiusAtt, periodicAtt, xPeriodAtt, yPeriodAtt;
	NcBool history, id, spec, conventions, source, parent_id;
	string history_str = "";
	string id_str = "";
	string parent_str ="";

	history_str += "MpasMaskCreator.x ";
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

	sprintf(mesh_spec_str, "%2.1lf", (double)MESH_SPEC);

	if (!(history = grid.add_att(   "history", history_str.c_str() ))) return NC_ERR;
	if (!(spec = grid.add_att(   "mesh_spec", mesh_spec_str ))) return NC_ERR;
	if (!(conventions = grid.add_att(   "Conventions", "MPAS" ))) return NC_ERR;
	if (!(source = grid.add_att(   "source", "MpasMeshConverter.x" ))) return NC_ERR;
	if (!(id = grid.add_att(   "file_id", id_str.c_str() ))) return NC_ERR;

	grid.close();

	// file closed when file obj goes out of scope
	return 0;
}/*}}}*/
int outputMaskFields( const string outputFilename) {/*{{{*/
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

	// open the mpas file
	NcFile grid(outputFilename.c_str(), NcFile::Write);

	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;

	// fetch dimensions
	NcDim *nCellsDim = grid.get_dim( "nCells" );
	NcDim *nVerticesDim = grid.get_dim( "nVertices" );
	NcDim *nEdgesDim = grid.get_dim( "nEdges" );
	NcDim *StrLenDim = grid.get_dim( "StrLen" );

	int StrLen = StrLenDim->size();

	//Define nc variables
	NcVar *tempVar;

	string fieldName;
	double *x, *y, *z, *lat, *lon;
	int *indices;
	int *counts;
	char *names;
	double *double_vals;
	int *int_vals;
	char *char_vals;
	int i, j, idx;

	if ( seedPoints.size() > 0 ) {
		if (!(tempVar = grid.add_var("cellSeedMask", ncInt, nCellsDim))) return NC_ERR;
		if (!tempVar->put(&cellSeedMask[0], nCells)) return NC_ERR;
	}

	if ( nRegions > 0 ) {
		NcDim *nRegionsDim = grid.get_dim( "nRegions" );
		NcDim *nRegionGroupsDim = grid.get_dim( "nRegionGroups" );
		NcDim *maxRegionsInGroupDim = grid.get_dim( "maxRegionsInGroup" );

		int nRegionGroups = nRegionGroupsDim->size();

		// Build and write region masks and counts
		counts = new int[nRegions];

		for ( i = 0; i < nRegions; i++){
			counts[i] = 0;

			for ( j = 0; j < nCells; j++){
				counts[i] += cellMasks[j * nRegions + i];
			}
		}

		if (!(tempVar = grid.add_var("regionCellMasks", ncInt, nCellsDim, nRegionsDim))) return NC_ERR;
		if (!tempVar->put(cellMasks, nCells, nRegions)) return NC_ERR;
		if (!(tempVar = grid.add_var("nCellsInRegion", ncInt, nRegionsDim))) return NC_ERR;
		if (!tempVar->put(counts, nRegions)) return NC_ERR;
		delete[] cellMasks;

		for ( i = 0; i < nRegions; i++){
			counts[i] = 0;

			for ( j = 0; j < nVertices; j++){
				counts[i] += vertexMasks[j * nRegions + i];
			}
		}
		if (!(tempVar = grid.add_var("regionVertexMasks", ncInt, nVerticesDim, nRegionsDim))) return NC_ERR;
		if (!tempVar->put(vertexMasks, nVertices, nRegions)) return NC_ERR;
		if (!(tempVar = grid.add_var("nVerticesInRegion", ncInt, nRegionsDim))) return NC_ERR;
		if (!tempVar->put(counts, nRegions)) return NC_ERR;
		delete[] counts;
		delete[] vertexMasks;

		// Build and write region names
		names = new char[nRegions * StrLen];
		for (i = 0; i < nRegions; i++){
			for ( j = 0; j < StrLen; j++){
				names[i * StrLen + j ] = '\0';
			}
		}
		for (str_itr = regionNames.begin(), i=0; str_itr != regionNames.end(); str_itr++, i++){
			snprintf(&names[i*StrLen], StrLen, "%s", (*str_itr).c_str());
		}
		if (!(tempVar = grid.add_var("regionNames", ncChar, nRegionsDim, StrLenDim))) return NC_ERR;
		if (!tempVar->put(names, nRegions, StrLen)) return NC_ERR;
		delete[] names;

		// Build and write region groups and group names
		indices = new int[nRegionGroups * maxRegionsInGroup];
		counts = new int[nRegionGroups];
		names = new char[nRegionGroups * StrLen];

		for ( i = 0; i < nRegionGroups; i++ ) {
			for ( j = 0; j < maxRegionsInGroup; j++ ){
				indices[i * maxRegionsInGroup + j] = 0;
			}
		}

		for ( vec_int_itr = regionsInGroup.begin(), i = 0; vec_int_itr != regionsInGroup.end(); vec_int_itr++, i++){
			counts[i] = (*vec_int_itr).size();
			for ( int_itr = (*vec_int_itr).begin(), j = 0; int_itr != (*vec_int_itr).end(); int_itr++, j++){
				indices[i * maxRegionsInGroup + j] = (*int_itr)+1;
			}
		}

		for (i = 0; i < nRegionGroups; i++){
			for ( j = 0; j < StrLen; j++){
				names[i * StrLen + j ] = '\0';
			}
		}
		for (str_itr = regionGroupNames.begin(), i=0; str_itr != regionGroupNames.end(); str_itr++, i++){
			snprintf(&names[i*StrLen], StrLen, "%s", (*str_itr).c_str());
		}

		if (!(tempVar = grid.add_var("nRegionsInGroup", ncInt, nRegionGroupsDim))) return NC_ERR;
		if (!tempVar->put(counts, nRegionGroups)) return NC_ERR;
		if (!(tempVar = grid.add_var("regionsInGroup", ncInt, nRegionGroupsDim, maxRegionsInGroupDim))) return NC_ERR;
		if (!tempVar->put(indices, nRegionGroups, maxRegionsInGroup)) return NC_ERR;
		if (!(tempVar = grid.add_var("regionGroupNames", ncChar, nRegionGroupsDim, StrLenDim))) return NC_ERR;
		if (!tempVar->put(names, nRegionGroups, StrLen)) return NC_ERR;

		delete[] indices;
		delete[] counts;
		delete[] names;

		// Write out region properties:
		for ( i = 0; i < regionProperties[0].size(); i++ ) {
			fieldName = get<0>( regionProperties[0][i] ) + "Regions";

			if ( get<2>( regionProperties[0][i] ) == num_double ) {
				double_vals = new double[nRegions];

				for ( j = 0; j < nRegions; j++ ) {
					double_vals[j] = 0.0;
				}

				for ( j = 0; j < regionProperties.size(); j++ ) {
					// Find the target propery for this region:

					for ( int k = 0; k < regionProperties[j].size(); k++ ) {
						if ( get<0>( regionProperties[0][i] ) == get<0>( regionProperties[j][k] ) ) {
							stringstream convert( get<1>( regionProperties[j][k] ) );
							convert >> double_vals[j];
						}
					}
				}

				if (!(tempVar = grid.add_var(fieldName.c_str(), ncDouble, nRegionsDim))) return NC_ERR;
				if (!tempVar->put(double_vals, nRegions)) return NC_ERR;

				delete[] double_vals;
			} else if ( get<2>( regionProperties[0][i] ) == num_int ) {
				int_vals = new int[nRegions];

				for ( j = 0; j < nRegions; j++ ) {
					int_vals[j] = 0;
				}

				for ( j = 0; j < regionProperties.size(); j++ ) {
					// Find the target propery for this region:

					for ( int k = 0; k < regionProperties[j].size(); k++ ) {
						if ( get<0>( regionProperties[0][i] ) == get<0>( regionProperties[j][k] ) ) {
							stringstream convert( get<1>( regionProperties[j][k] ) );
							convert >> int_vals[j];
						}
					}
				}

				if (!(tempVar = grid.add_var(fieldName.c_str(), ncInt, nRegionsDim))) return NC_ERR;
				if (!tempVar->put(int_vals, nRegions)) return NC_ERR;

				delete[] int_vals;
			} else if ( get<2>( regionProperties[0][i] ) == text ) {
				char_vals = new char[nRegions * StrLen];

				for ( j = 0; j < nRegions * StrLen; j++ ) {
					char_vals[j] = '\0';
				}

				for ( j = 0; j < regionProperties.size(); j++ ) {
					// Find the target propery for this region:

					for ( int k = 0; k < regionProperties[j].size(); k++ ) {
						if ( get<0>( regionProperties[0][i] ) == get<0>( regionProperties[j][k] ) ) {
							snprintf(&char_vals[j * StrLen], StrLen, "%s", get<1>( regionProperties[j][k] ).c_str() );
						}
					}
				}

				if (!(tempVar = grid.add_var(fieldName.c_str(), ncChar, nRegionsDim, StrLenDim))) return NC_ERR;
				if (!tempVar->put(char_vals, nRegions, StrLen)) return NC_ERR;

				delete[] char_vals;
			}
		}
	}

	if ( nTransects > 0 ) {
		NcDim *nTransectsDim = grid.get_dim( "nTransects" );
		NcDim *nTransectGroupsDim = grid.get_dim( "nTransectGroups" );
		NcDim *maxTransectsInGroupDim = grid.get_dim( "maxTransectsInGroup" );
		NcDim *maxCellsInTransectDim = grid.get_dim( "maxCellsInTransect" );
		NcDim *maxVerticesInTransectDim = grid.get_dim( "maxVerticesInTransect" );
		NcDim *maxEdgesInTransectDim = grid.get_dim( "maxEdgesInTransect" );
		int nTransectGroups = nTransectGroupsDim->size();

		// Write out cell transect masks
		indices = new int[nCells * nTransects];
		for ( i = 0; i < nCells; i++ ){
			for (j = 0; j < nTransects; j++ ) {
				indices[i * nTransects + j ] = 0;
			}
		}

		for ( i = 0; i < cellPaths.size(); i++ ) {
			for ( int_itr = cellPaths[i].begin(); int_itr != cellPaths[i].end(); int_itr++ ) {
				indices[ (*int_itr) * nTransects + i ] = 1;
			}
		}

		if (!(tempVar = grid.add_var("transectCellMasks", ncInt, nCellsDim, nTransectsDim))) return NC_ERR;
		if (!tempVar->put(indices, nCells, nTransects)) return NC_ERR;

		delete[] indices;

		// Write out cell ID list
		indices = new int[nTransects * maxCellsInTransect];
		for ( i = 0; i < nTransects * maxCellsInTransect; i++ ) {
			indices[i] = 0;
		}

		for ( i = 0; i < cellPaths.size(); i++ ) {
			for (j = 0; j < cellPaths[i].size(); j++ ) {
				indices[i * maxCellsInTransect + j] = cellPaths[i][j] + 1;
			}
		}

		if (!(tempVar = grid.add_var("transectCellGlobalIDs", ncInt, nTransectsDim, maxCellsInTransectDim))) return NC_ERR;
		if (!tempVar->put(indices, nTransects, maxCellsInTransect)) return NC_ERR;

		delete[] indices;

		// Write out vertex transect masks
		indices = new int[nVertices * nTransects];
		for ( i = 0; i < nVertices; i++ ){
			for (j = 0; j < nTransects; j++ ) {
				indices[i * nTransects + j ] = 0;
			}
		}

		for ( i = 0; i < vertexPaths.size(); i++ ) {
			for ( int_itr = vertexPaths[i].begin(); int_itr != vertexPaths[i].end(); int_itr++ ) {
				indices[ (*int_itr) * nTransects + i ] = 1;
			}
		}

		if (!(tempVar = grid.add_var("transectVertexMasks", ncInt, nVerticesDim, nTransectsDim))) return NC_ERR;
		if (!tempVar->put(indices, nVertices, nTransects)) return NC_ERR;

		delete[] indices;

		// Write out vertex ID list
		indices = new int[nTransects * maxVerticesInTransect];
		for ( i = 0; i < nTransects * maxVerticesInTransect; i++ ) {
			indices[i] = 0;
		}

		for ( i = 0; i < vertexPaths.size(); i++ ) {
			for (j = 0; j < vertexPaths[i].size(); j++ ) {
				indices[i * maxVerticesInTransect + j] = vertexPaths[i][j] + 1;
			}
		}

		if (!(tempVar = grid.add_var("transectVertexGlobalIDs", ncInt, nTransectsDim, maxVerticesInTransectDim))) return NC_ERR;
		if (!tempVar->put(indices, nTransects, maxVerticesInTransect)) return NC_ERR;

		delete[] indices;

		// Write out edge transect masks
		indices = new int[nEdges * nTransects];
		for ( i = 0; i < nEdges; i++ ){
			for (j = 0; j < nTransects; j++ ) {
				indices[i * nTransects + j ] = 0;
			}
		}

		for ( i = 0; i < edgePaths.size(); i++ ) {
			for ( int_itr = edgePaths[i].begin(); int_itr != edgePaths[i].end(); int_itr++ ) {
				indices[ (*int_itr) * nTransects + i ] = 1;
			}
		}

		if (!(tempVar = grid.add_var("transectEdgeMasks", ncInt, nEdgesDim, nTransectsDim))) return NC_ERR;
		if (!tempVar->put(indices, nEdges, nTransects)) return NC_ERR;

		delete[] indices;

		// Build and write transect names
		names = new char[nTransects * StrLen];
		for (i = 0; i < nTransects; i++){
			for ( j = 0; j < StrLen; j++){
				names[i * StrLen + j ] = '\0';
			}
		}
		for (str_itr = transectNames.begin(), i=0; str_itr != transectNames.end(); str_itr++, i++){
			snprintf(&names[i*StrLen], StrLen, "%s", (*str_itr).c_str());
		}
		if (!(tempVar = grid.add_var("transectNames", ncChar, nTransectsDim, StrLenDim))) return NC_ERR;
		if (!tempVar->put(names, nTransects, StrLen)) return NC_ERR;
		delete[] names;

		// Build and write region groups and group names
		indices = new int[nTransectGroups * maxTransectsInGroup];
		counts = new int[nTransectGroups];
		names = new char[nTransectGroups * StrLen];

		for ( i = 0; i < nTransectGroups; i++ ) {
			for ( j = 0; j < maxTransectsInGroup; j++ ){
				indices[i * maxTransectsInGroup + j] = 0;
			}
		}

		for ( vec_int_itr = transectsInGroup.begin(), i = 0; vec_int_itr != transectsInGroup.end(); vec_int_itr++, i++){
			counts[i] = (*vec_int_itr).size();
			for ( int_itr = (*vec_int_itr).begin(), j = 0; int_itr != (*vec_int_itr).end(); int_itr++, j++){
				indices[i * maxTransectsInGroup + j] = (*int_itr)+1;
			}
		}

		for (i = 0; i < nTransectGroups; i++){
			for ( j = 0; j < StrLen; j++){
				names[i * StrLen + j ] = '\0';
			}
		}
		for (str_itr = transectGroupNames.begin(), i=0; str_itr != transectGroupNames.end(); str_itr++, i++){
			snprintf(&names[i*StrLen], StrLen, "%s", (*str_itr).c_str());
		}

		if (!(tempVar = grid.add_var("nTransectsInGroup", ncInt, nTransectGroupsDim))) return NC_ERR;
		if (!tempVar->put(counts, nTransectGroups)) return NC_ERR;
		if (!(tempVar = grid.add_var("transectsInGroup", ncInt, nTransectGroupsDim, maxTransectsInGroupDim))) return NC_ERR;
		if (!tempVar->put(indices, nTransectGroups, maxTransectsInGroup)) return NC_ERR;
		if (!(tempVar = grid.add_var("transectGroupNames", ncChar, nTransectGroupsDim, StrLenDim))) return NC_ERR;
		if (!tempVar->put(names, nTransectGroups, StrLen)) return NC_ERR;

		delete[] indices;
		delete[] counts;
		delete[] names;

		// Write out signed edge transect masks
		indices = new int[nEdges * nTransects];
		for ( i = 0; i < nEdges; i++ ){
			for (j = 0; j < nTransects; j++ ) {
				indices[i * nTransects + j ] = 0;
			}
		}

		if ( edgePaths.size() != edgePathSigns.size() ) {
			cout << "ERROR: Number of paths are not equal in edgePaths and edgePathSigns." << endl;
			cout << "       " << edgePaths.size() << " != " << edgePathSigns.size() << endl;
			cout << "Exiting..." << endl;
			exit(1);
		}

		for ( i = 0; i < edgePaths.size(); i++ ) {
			if ( edgePaths[i].size() != edgePathSigns[i].size() ) {
				cout << "ERROR: Number of edges are not equal in path " << i << endl;
				cout << "       " << edgePaths[i].size() << " != " << edgePathSigns[i].size() << endl;
				cout << "Exiting..." << endl;
				exit(1);
			}
			for ( j = 0; j < edgePaths[i].size(); j++ ) {
				indices[ edgePaths[i][j] * nTransects + i ] = edgePathSigns[i][j];
			}
		}

		if (!(tempVar = grid.add_var("transectEdgeMaskSigns", ncInt, nEdgesDim, nTransectsDim))) return NC_ERR;
		if (!tempVar->put(indices, nEdges, nTransects)) return NC_ERR;

		delete[] indices;

		// Write out edge ID list
		indices = new int[nTransects * maxEdgesInTransect];
		for ( i = 0; i < nTransects * maxEdgesInTransect; i++ ) {
			indices[i] = 0;
		}

		for ( i = 0; i < edgePaths.size(); i++ ) {
			for (j = 0; j < edgePaths[i].size(); j++ ) {
				indices[i * maxEdgesInTransect + j] = edgePaths[i][j] + 1;
			}
		}

		if (!(tempVar = grid.add_var("transectEdgeGlobalIDs", ncInt, nTransectsDim, maxEdgesInTransectDim))) return NC_ERR;
		if (!tempVar->put(indices, nTransects, maxEdgesInTransect)) return NC_ERR;

		delete[] indices;

		// Write out transect properties:
		for ( i = 0; i < transectProperties[0].size(); i++ ) {
			fieldName = get<0>( transectProperties[0][i] ) + "Transects";

			if ( get<2>( transectProperties[0][i] ) == num_double ) {
				double_vals = new double[nTransects];

				for ( j = 0; j < nTransects; j++ ) {
					double_vals[j] = 0.0;
				}

				for ( j = 0; j < transectProperties.size(); j++ ) {
					// Find the target propery for this region:

					for ( int k = 0; k < transectProperties[j].size(); k++ ) {
						if ( get<0>( transectProperties[0][i] ) == get<0>( transectProperties[j][k] ) ) {
							stringstream convert( get<1>( transectProperties[j][k] ) );
							convert >> double_vals[j];
						}
					}
				}

				if (!(tempVar = grid.add_var(fieldName.c_str(), ncDouble, nTransectsDim))) return NC_ERR;
				if (!tempVar->put(double_vals, nTransects)) return NC_ERR;

				delete[] double_vals;
			} else if ( get<2>( transectProperties[0][i] ) == num_int ) {
				int_vals = new int[nTransects];

				for ( j = 0; j < nTransects; j++ ) {
					int_vals[j] = 0;
				}

				for ( j = 0; j < transectProperties.size(); j++ ) {
					// Find the target propery for this region:

					for ( int k = 0; k < transectProperties[j].size(); k++ ) {
						if ( get<0>( transectProperties[0][i] ) == get<0>( transectProperties[j][k] ) ) {
							stringstream convert( get<1>( transectProperties[j][k] ) );
							convert >> int_vals[j];
						}
					}
				}

				if (!(tempVar = grid.add_var(fieldName.c_str(), ncInt, nTransectsDim))) return NC_ERR;
				if (!tempVar->put(int_vals, nTransects)) return NC_ERR;

				delete[] int_vals;
			} else if ( get<2>( transectProperties[0][i] ) == text ) {
				char_vals = new char[nTransects * StrLen];

				for ( j = 0; j < nTransects * StrLen; j++ ) {
					char_vals[j] = '\0';
				}

				for ( j = 0; j < transectProperties.size(); j++ ) {
					// Find the target propery for this region:

					for ( int k = 0; k < transectProperties[j].size(); k++ ) {
						if ( get<0>( transectProperties[0][i] ) == get<0>( transectProperties[j][k] ) ) {
							snprintf(&char_vals[j * StrLen], StrLen, "%s", get<1>( transectProperties[j][k] ).c_str() );
						}
					}
				}

				if (!(tempVar = grid.add_var(fieldName.c_str(), ncChar, nTransectsDim, StrLenDim))) return NC_ERR;
				if (!tempVar->put(char_vals, nTransects, StrLen)) return NC_ERR;

				delete[] char_vals;
			}
		}
	}

	if ( nPoints > 0 ) {
		NcDim *nPointsDim = grid.get_dim( "nPoints" );
		NcDim *nPointGroupsDim = grid.get_dim( "nPointGroups" );
		NcDim *maxPointsInGroupDim = grid.get_dim( "maxPointsInGroup" );
		int nPointGroups = nPointGroupsDim->size();

#ifdef _WRITE_POINT_MASK
		// Build and write point mask for cells, to see where the points live.
		cellPointMasks = new int[nCells];
		for ( i = 0; i < nCells; i++ ) {
			cellPointMasks[i] = 0;
		}
		for (i = 0; i < nPoints; i++){
			idx = pointCellIndices[i];
			cellPointMasks[ idx - 1 ] = 1;
		}

		if (!(tempVar = grid.add_var("cellPointMask", ncInt, nCellsDim))) return NC_ERR;
		if (!tempVar->put(cellPointMasks, nCells)) return NC_ERR;
#endif

		// Build and write point index information
		if (!(tempVar = grid.add_var("pointCellGlobalID", ncInt, nPointsDim))) return NC_ERR;
		if (!tempVar->put(pointCellIndices, nPoints)) return NC_ERR;
		delete[] pointCellIndices;

		if (!(tempVar = grid.add_var("pointVertexGlobalID", ncInt, nPointsDim))) return NC_ERR;
		if (!tempVar->put(pointVertexIndices, nPoints)) return NC_ERR;
		delete[] pointVertexIndices;

		// Build and write point names
		names = new char[nPoints * StrLen];
		for (i = 0; i < nPoints; i++){
			for ( j = 0; j < StrLen; j++){
				names[i * StrLen + j ] = '\0';
			}
		}
		for (str_itr = pointNames.begin(), i=0; str_itr != pointNames.end(); str_itr++, i++){
			snprintf(&names[i*StrLen], StrLen, "%s", (*str_itr).c_str());
		}
		if (!(tempVar = grid.add_var("pointNames", ncChar, nPointsDim, StrLenDim))) return NC_ERR;
		if (!tempVar->put(names, nPoints, StrLen)) return NC_ERR;
		delete[] names;

		// Build and write point groups and group names
		indices = new int[nPointGroups * maxPointsInGroup];
		counts = new int[nPointGroups];
		names = new char[nPointGroups * StrLen];

		for ( i = 0; i < nPointGroups; i++ ) {
			for ( j = 0; j < maxPointsInGroup; j++ ){
				indices[i * maxPointsInGroup + j] = 0;
			}
		}

		for ( vec_int_itr = pointsInGroup.begin(), i = 0; vec_int_itr != pointsInGroup.end(); vec_int_itr++, i++){
			counts[i] = (*vec_int_itr).size();
			for ( int_itr = (*vec_int_itr).begin(), j = 0; int_itr != (*vec_int_itr).end(); int_itr++, j++){
				indices[i * maxPointsInGroup + j] = (*int_itr)+1;
			}
		}

		for (i = 0; i < nPointGroups; i++){
			for ( j = 0; j < StrLen; j++){
				names[i * StrLen + j ] = '\0';
			}
		}
		for (str_itr = pointGroupNames.begin(), i=0; str_itr != pointGroupNames.end(); str_itr++, i++){
			snprintf(&names[i*StrLen], StrLen, "%s", (*str_itr).c_str());
		}

		if (!(tempVar = grid.add_var("nPointsInGroup", ncInt, nPointGroupsDim))) return NC_ERR;
		if (!tempVar->put(counts, nPointGroups)) return NC_ERR;
		if (!(tempVar = grid.add_var("pointsInGroup", ncInt, nPointGroupsDim, maxPointsInGroupDim))) return NC_ERR;
		if (!tempVar->put(indices, nPointGroups, maxPointsInGroup)) return NC_ERR;
		if (!(tempVar = grid.add_var("pointGroupNames", ncChar, nPointGroupsDim, StrLenDim))) return NC_ERR;
		if (!tempVar->put(names, nPointGroups, StrLen)) return NC_ERR;

		delete[] indices;
		delete[] counts;
		delete[] names;

		// Build point coordinates and output
		x = new double[nPoints];
		y = new double[nPoints];
		z = new double[nPoints];
		lat = new double[nPoints];
		lon = new double[nPoints];

		for (pnt_itr = pointLocations.begin(), i = 0; pnt_itr != pointLocations.end(); pnt_itr++, i++){
			x[i] = (*pnt_itr).x;
			y[i] = (*pnt_itr).y;
			z[i] = (*pnt_itr).z;
			if ( !lonRangePositive ) (*pnt_itr).setPositiveLonRange(true);
			lat[i] = (*pnt_itr).getLat();
			lon[i] = (*pnt_itr).getLon();
			if ( !lonRangePositive ) (*pnt_itr).setPositiveLonRange(lonRangePositive);
		}

		if (!(tempVar = grid.add_var("xPoint", ncDouble, nPointsDim))) return NC_ERR;
		if (!tempVar->put(x, nPoints)) return NC_ERR;
		if (!(tempVar = grid.add_var("yPoint", ncDouble, nPointsDim))) return NC_ERR;
		if (!tempVar->put(y, nPoints)) return NC_ERR;
		if (!(tempVar = grid.add_var("zPoint", ncDouble, nPointsDim))) return NC_ERR;
		if (!tempVar->put(z, nPoints)) return NC_ERR;
		if (!(tempVar = grid.add_var("latPoint", ncDouble, nPointsDim))) return NC_ERR;
		if (!tempVar->put(lat, nPoints)) return NC_ERR;
		if (!(tempVar = grid.add_var("lonPoint", ncDouble, nPointsDim))) return NC_ERR;
		if (!tempVar->put(lon, nPoints)) return NC_ERR;
		delete[] x;
		delete[] y;
		delete[] z;
		delete[] lat;
		delete[] lon;

		// Write out point properties:
		for ( i = 0; i < pointProperties[0].size(); i++ ) {
			fieldName = get<0>( pointProperties[0][i] ) + "Points";

			if ( get<2>( pointProperties[0][i] ) == num_double ) {
				double_vals = new double[nPoints];

				for ( j = 0; j < nPoints; j++ ) {
					double_vals[j] = 0.0;
				}

				for ( j = 0; j < pointProperties.size(); j++ ) {
					// Find the target propery for this region:

					for ( int k = 0; k < pointProperties[j].size(); k++ ) {
						if ( get<0>( pointProperties[0][i] ) == get<0>( pointProperties[j][k] ) ) {
							stringstream convert( get<1>( pointProperties[j][k] ) );
							convert >> double_vals[j];
						}
					}
				}

				if (!(tempVar = grid.add_var(fieldName.c_str(), ncDouble, nPointsDim))) return NC_ERR;
				if (!tempVar->put(double_vals, nPoints)) return NC_ERR;

				delete[] double_vals;
			} else if ( get<2>( pointProperties[0][i] ) == num_int ) {
				int_vals = new int[nPoints];

				for ( j = 0; j < nPoints; j++ ) {
					int_vals[j] = 0;
				}

				for ( j = 0; j < pointProperties.size(); j++ ) {
					// Find the target propery for this region:

					for ( int k = 0; k < pointProperties[j].size(); k++ ) {
						if ( get<0>( pointProperties[0][i] ) == get<0>( pointProperties[j][k] ) ) {
							stringstream convert( get<1>( pointProperties[j][k] ) );
							convert >> int_vals[j];
						}
					}
				}

				if (!(tempVar = grid.add_var(fieldName.c_str(), ncInt, nPointsDim))) return NC_ERR;
				if (!tempVar->put(int_vals, nPoints)) return NC_ERR;

				delete[] int_vals;
			} else if ( get<2>( pointProperties[0][i] ) == text ) {
				char_vals = new char[nPoints * StrLen];

				for ( j = 0; j < nPoints * StrLen; j++ ) {
					char_vals[j] = '\0';
				}

				for ( j = 0; j < pointProperties.size(); j++ ) {
					// Find the target propery for this region:

					for ( int k = 0; k < pointProperties[j].size(); k++ ) {
						if ( get<0>( pointProperties[0][i] ) == get<0>( pointProperties[j][k] ) ) {
							snprintf(&char_vals[j * StrLen], StrLen, "%s", get<1>( pointProperties[j][k] ).c_str() );
						}
					}
				}

				if (!(tempVar = grid.add_var(fieldName.c_str(), ncChar, nPointsDim, StrLenDim))) return NC_ERR;
				if (!tempVar->put(char_vals, nPoints, StrLen)) return NC_ERR;

				delete[] char_vals;
			}
		}
	}

	grid.close();

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

