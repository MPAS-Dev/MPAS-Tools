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

#include "netcdf_utils.h"
#include "pnt.h"
#include "edge.h"

#define MESH_SPEC 1.0
#define ID_LEN 10

//#define _WRITE_POINT_MASK

using namespace std;

int nCells, nVertices;
int nRegions = 0;
int nPoints = 0;
int maxRegPolygons = 0;
int maxRegVertices = 0;
int maxRegionsInGroup = 0;
int maxPointsInGroup = 0;
int maxEdges;
double sphereRadius;
bool spherical, periodic;
string in_history = "";
string in_file_id = "";
string in_parent_id = "";

// Mask and location information {{{

int *cellPointMasks;
int *cellMasks, *vertexMasks, *pointCellIndices, *pointVertexIndices;
vector<pnt> cells;
vector<pnt> vertices;
vector<pnt> pointLocations;
vector< vector< vector<pnt> > > regionPolygons;
vector<string> regionNames;
vector<string> pointNames;
vector< vector<int> > regionsInGroup;
vector< vector<int> > pointsInGroup;
vector<string> regionGroupNames;
vector<string> pointGroupNames;
vector< vector< vector<double> > > polygonConstants;
vector< vector< vector<double> > > polygonMultiples;

// }}}

// Iterators {{{
vector<pnt>::iterator pnt_itr, vert_itr;
vector<string>::iterator str_itr;
vector< vector<int> >::iterator mask_itr, vec_int_itr;
vector< vector< vector<pnt> > >::iterator reg_itr;
vector< vector<pnt> >::iterator poly_itr;
vector<int>::iterator int_itr;
// }}}

/* Building/Reading functions {{{ */
int readGridInfo(const string inputFilename);
int readCells(const string inputFilename);
int readVertices(const string inputFilename);
int resetFeatureInfo();
int getFeatureInfo(const string featureFilename);
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

string gen_random(const int len);

int main ( int argc, char *argv[] ) {
	int error;
	string masks_name = "features.geojson";
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
	//
	//  If the input file was not specified, get it now.
	//
	if ( argc <= 1 )
	{
		cout << "\n";
		cout << "MPAS_MASK_CREATOR:\n";
		cout << "  Please enter the MPAS NetCDF input filename.\n";

		cin >> in_name;

		cout << "\n";
		cout << "MPAS_MASK_CREATOR:\n";
		cout << "  Please enter the MPAS masks NetCDF output filename.\n";

		cin >> out_name;

		cout << "\n";
		cout << "MPAS_MASK_CREATOR:\n";
		cout << "  Please enter a file containing a set of regions to create masks from.\n";

		cin >> masks_name;
	}
	else if (argc == 2)
	{
		in_name = argv[1];

		cout << "\n";
		cout << "MPAS_MESH_CONVERTER:\n";
		cout << "  Output name not specified. Using default of " << out_name << endl;
		cout << "  No features file specified. Using default of " << masks_name << endl;
	}
	else if (argc == 3)
	{
		in_name = argv[1];
		out_name = argv[2];
		cout << "\n";
		cout << "MPAS_MESH_CONVERTER:\n";
		cout << "  No features file specified. Using default of " << masks_name << endl;
	} else if (argc > 3)
	{
		in_name = argv[1];
		out_name = argv[2];
	}

	if(in_name == out_name){
		cout << "   ERROR: Input and output names are the same." << endl;
		return 1;
	}


	srand(time(NULL));

	cout << "Reading input grid." << endl;
	error = readGridInfo(in_name);
	if(error) return 1;

	error = resetFeatureInfo();

	cout << "Building feature information." << endl;
	if (argc > 3) {
		for ( int i = 3; i < argc; i++){
			masks_name = argv[i];
			error = getFeatureInfo(masks_name);
			if(error) return 1;
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

	cout << "Deleting cell center information" << endl;
	cells.clear();

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

	cout << "Deleting vertex information" << endl;
	vertices.clear();

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

/* Building/Ordering functions {{{ */
int readGridInfo(const string inputFilename){/*{{{*/
#ifdef _DEBUG
	cout << endl << endl << "Begin function: readGridInput" << endl << endl;
#endif

	nCells = netcdf_mpas_read_dim(inputFilename, "nCells");
	nVertices = netcdf_mpas_read_dim(inputFilename, "nVertices");
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
int resetFeatureInfo(){/*{{{*/
	regionPolygons.clear();
	regionNames.clear();
	regionsInGroup.clear();
	regionGroupNames.clear();
	polygonConstants.clear();
	polygonMultiples.clear();
	pointLocations.clear();

	return 0;
}/*}}}*/
int getFeatureInfo(const string featureFilename){/*{{{*/
	ifstream json_file(featureFilename);	
	Json::Value root;
	string groupName, tempGroupName;
	vector<int> pointIndices;
	vector<int> regionIndices;
	vector<pnt> regionVertices;
	vector< vector<pnt> > polygonList;
	bool addGroupToRegList = true;
	bool addGroupToPntList = true;
	bool addRegToGroup;
	bool createAllRegGroup;

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
	regionIndices.clear();

	for ( int i = 0; i < root["features"].size(); i++ ){
		Json::Value feature = root["features"][i];
		string featureName = feature["properties"]["name"].asString();

		if ( feature["properties"]["object"].asString() == "region" ) {
			bool add_region = true;
			int regionIdx, idx;

			for ( str_itr = regionNames.begin(), idx = 0; str_itr != regionNames.end(); str_itr++, idx++){
				if ( (*str_itr) == featureName ) {
					add_region = false;
					regionIdx = idx;
				}
			}

			if ( add_region ) {
				polygonList.clear();
				regionIdx = regionNames.size();
#ifdef _DEBUG
				cout << "Adding region: " << featureName << " with index " << regionNames.size() << endl;
#endif
				regionNames.push_back(featureName);

				Json::Value geometry = feature["geometry"];
				Json::Value coordinates = geometry["coordinates"];

				if ( geometry["type"].asString() == "Polygon" ) {

					regionVertices.clear();
					for (int i = 0; i < coordinates[0].size(); i++){
						double lon = coordinates[0][i][0].asDouble() * M_PI/180.0;
						double lat = coordinates[0][i][1].asDouble() * M_PI/180.0;

						if ( lon < 0.0 ) {
							lon = lon + (2.0 * M_PI);
						}

#ifdef _DEBUG
						cout << " Added a region vertex: " << lat << ", " << lon << endl;
#endif
						pnt point = pntFromLatLon(lat, lon);
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

							if ( lon < 0.0 ) {
								lon = lon + (2.0 * M_PI);
							}

#ifdef _DEBUG
							cout << " Added a region vertex: " << lat << ", " << lon << endl;
#endif
							pnt point = pntFromLatLon(lat, lon);
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
						tempGroupName = "regionGroup" + to_string( regionGroupNames.size() + 1 );
						regionGroupNames.push_back(tempGroupName);
					} else {
						regionGroupNames.push_back(groupName);
					}
				}

			}

		} else if ( feature["properties"]["object"].asString() == "point" ) {
			bool addPointToGroup;
			bool add_point = true;
			int pointIdx, idx;

			for ( str_itr = pointNames.begin(), idx = 0; str_itr != pointNames.end(); str_itr++, idx++ ) {
				if ( (*str_itr) == featureName ) {
					add_point = false;
					pointIdx = idx;
				}
			}

			if ( add_point ) {
				double lon = feature["geometry"]["coordinates"][0].asDouble() * M_PI/180.0 ;
				double lat = feature["geometry"]["coordinates"][1].asDouble() * M_PI/180.0;

				if ( lon < 0.0 ) {
					lon = lon + ( 2.0 * M_PI );
				}

				pointIdx = pointNames.size();

#ifdef _DEBUG
				cout << "Adding point: " << featureName << " with index " << pointNames.size() << endl;
#endif

				pnt point = pntFromLatLon(lat, lon);
				point.idx = pointIdx;

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
						tempGroupName = "pointGroup" + to_string( pointGroupNames.size() + 1 );
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

	if ( (int)pointIndices.size() > 0 ) {
		maxPointsInGroup = max(maxPointsInGroup, (int)pointIndices.size());
		pointsInGroup.push_back(pointIndices);
	}

	return 0;
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
	
	// open the scvtmesh file
	NcFile grid(outputFilename.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);

	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;
	
	// define dimensions
	NcDim *tempDim;

	// write dimensions
	if (!(tempDim = grid.add_dim("nCells", nCells))) return NC_ERR;
	if (!(tempDim = grid.add_dim("nVertices", nVertices))) return NC_ERR;

	if ( nRegions > 0 ) {
		if (!(tempDim = grid.add_dim("nRegions", nRegions))) return NC_ERR;
		if (!(tempDim = grid.add_dim("nRegionGroups", regionGroupNames.size()))) return NC_ERR;
		if (!(tempDim = grid.add_dim("maxRegionsInGroup", maxRegionsInGroup))) return NC_ERR;
	}

	if ( nPoints > 0 ) {
		if (!(tempDim = grid.add_dim("nPoints", nPoints))) return NC_ERR;
		if (!(tempDim = grid.add_dim("nPointGroups", pointGroupNames.size()))) return NC_ERR;
		if (!(tempDim = grid.add_dim("maxPointsInGroup", maxPointsInGroup))) return NC_ERR;
	}

	if (!(tempDim = grid.add_dim("StrLen", 64))) return NC_ERR;

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
	
	// open the scvtmesh file
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
	
	// open the scvtmesh file
	NcFile grid(outputFilename.c_str(), NcFile::Write);
	
	// check to see if the file was opened
	if(!grid.is_valid()) return NC_ERR;
	
	// fetch dimensions
	NcDim *nCellsDim = grid.get_dim( "nCells" );
	NcDim *nVerticesDim = grid.get_dim( "nVertices" );
	NcDim *StrLenDim = grid.get_dim( "StrLen" );

	int StrLen = StrLenDim->size();

	//Define nc variables
	NcVar *tempVar;

	double *x, *y, *z, *lat, *lon;
	int *indices;
	int *counts;
	char *names;
	int i, j, idx;


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
			lat[i] = (*pnt_itr).getLat();
			lon[i] = (*pnt_itr).getLon();
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

