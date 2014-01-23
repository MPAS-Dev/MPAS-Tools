#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#ifdef _LINUX
	#include <unistd.h>
#endif

#include "constants.h"
#include "vec_utils.h"
#include "netcdf_utils.h"

using namespace std;

#ifdef _MACOS
	#include <GLUT/glut.h>
#elif _LINUX
	#include <GL/glut.h>
#endif

int main ( int argc, char *argv[] );
void display ( );
void mouse ( int btn, int state, int x, int y );
void gl_init ( );
void myReshape ( int w, int h );
void single_screenshot ( );
void timestamp ( );

void build_connectivity();
void setup_ranges();
void build_range(int id);

void rescale_cells_and_vertices();

void draw_cells();
void draw_triangles();
void draw_edges();

void draw_cell_lines();
void draw_triangle_lines();
void draw_edge_lines();

void build_regions();
void draw_regions();

void color_mesh();
void color_cells();
void color_triangles();
void color_edges();

void arrowKeys( int a_keys, int x, int y );
void keyPressed( unsigned char key, int x, int y );
void translateView ( double updown, double leftright);
void polarView( double distance, double twist, double elevation, double azimuth );
void hsv_to_rgb(float h, float s, float v, float& r, float& g, float& b);

double getLat(double x, double y, double z);
double getLon(double x, double y, double z);
double getX(double lat, double lon);
double getY(double lat, double lon);
double getZ(double lat, double lon);

void drawSphere(double r, int lats, int longs);

void screenshot (char filename[160],int x, int y) ;
void control_sequence();

//
//  Global data.
//
string filename;
static GLint axis = 1;
GLint window;
int drawing = 1;
int draw_lines = 1;
bool on_sphere;
bool draw_sphere;
double sphere_radius;
int color_bar = 0;
double missing_value = -1e34;

double line_factor = 1.008;
double region_line_factor = 1.01;
double region_center_factor = 1.013;
double range_factor = 0.80;
double zNear = 0.1;
double zFar = 2.5;

double projUpDown = 0.0;
double projLeftRight = 0.0;
double projDistance	= 3.0;
double projTwist		= 0.0;
double projElevation	= 0.0;
double projAzimuth		= 0.0;

bool single_ss = false;
int loop_it = 0;
static GLfloat theta[3] = { 0.0, 0.0, 0.0 };
double theta_speed = 0.020;

int ntime;
int nvertlevels;
int ncells;
int nvertices;
int nedges;
int maxedges;
int pixel_height;
int pixel_width;

int edge_field = -1;
int cell_field = -1;
int vertex_field = -1;
int cur_level = 0;
int cur_time = 0;
int cur_screenshot = 1;

double *xcell;
double *ycell;
double *zcell;
double *xvertex;
double *yvertex;
double *zvertex;

int *verticesoncell;
int *verticesonedge;
int *cellsonvertex;
int *cellsonedge;
int *nedgesoncell;

double *cell_values;
double *triangle_values;
double *edge_values;

vector<bool> cell_cells_save;
vector<GLfloat> cell_cells;
vector<GLfloat> vertex_cells;
vector<GLfloat> edge_cells;
vector<GLfloat> cell_lines;
vector<GLfloat> vertex_lines;
vector<GLfloat> edge_lines;
vector<GLfloat> cell_colors;
vector<GLfloat> vertex_colors;
vector<GLfloat> edge_colors;

vector< vector<double> > ranges; // 0 - min, 1 - max
vector<double> hard_ranges;

vector<GLfloat> region_centers;
vector<GLfloat> region_lines;
vector<double> region_radii;
vector<GLfloat> region_angles;
bool regions_built = false;
bool region_draw = false;
int region_line_div = 180;

double xyz_center[3];
double xyz_max[3];
double xyz_min[3];
double xyz_range[3];
double xyz_scale;

int main ( int argc, char *argv[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    MAIN is the main program for TRIANGULATION_FACES.
	//
	//  Discussion:
	//
	//    This program reads certain information from an MPAS NETCDF grid file,
	//    and displays the faces of the triangulation.
	//
	//  Usage:
	//
	//    triangulation_faces file.nc
	//
	//    where
	//
	//    * file.nc is an MPAS NETCDF grid file.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    01 January 2011
	//
	//  Author:
	//
	//    John Burkardt, Geoff Womeldorff, Doug Jacobsen
	//
	//  Reference:
	//
	//    Edward Angel,
	//    Interactive Computer Graphics:
	//    A Top-Down Approach with OpenGL,
	//    Second Edition,
	//    Addison Wesley, 2000.
	//
	int cell;
	int edge;
	int i;
	int v;

	cout << "\n";
	timestamp ( );

	cout << "\n";
	cout << "MPAS_DRAW:\n";
	cout << "  C++ version\n";
	cout << "  Read an MPAS NETCDF grid file\n";
	cout << "  Visualize the mpas grid/output file.\n";
	cout << "\n";
	cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
	//
	//  If the input file was not specified, get it now.
	//
	if ( argc <= 1 )
	{
		cout << "\n";
		cout << "MPAS_DRAW:\n";
		cout << "  Please enter the MPAS NETCDF grid filename.\n";

		cin >> filename;
	}
	else if (argc == 2)
	{
		filename = argv[1];
	}
	else if (argc == 3)
	{
		filename = argv[1];
		single_ss = true;
	}

	build_connectivity();
	setup_ranges();
	color_cells();
	color_triangles();
	color_edges();

	//
	//  Hand things over to OpenGL.
	//
	glutInit ( &argc, argv );
	glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
	glutInitWindowSize ( 800, 800 );
	glutInitWindowPosition ( 0, 0 );
	window = glutCreateWindow ( filename.c_str ( ) );
	glutReshapeFunc ( myReshape );
	glutDisplayFunc ( display );
	glutIdleFunc ( single_screenshot );
	glutMouseFunc ( mouse );
	glutKeyboardFunc( keyPressed );
	glutSpecialFunc( arrowKeys );

	//
	//  Enable hidden surface removal.
	//
	glEnable ( GL_DEPTH_TEST );

	gl_init ( );

	glutMainLoop ( );


	//
	//  Things that won't actually happen because we never return from glutMainLoop:
	//
	//
	//  Terminate.
	//
	cout << "\n";
	cout << "TRIANGULATION_FACES::\n";
	cout << "  Normal end of execution.\n";

	cout << "\n";
	timestamp ( );

	return 0;
}/*}}}*/
void display ( ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    DISPLAY generates the graphics output.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    John Burkardt, Geoff Womeldorff, Doug Jacobsen
	//
	//  Reference:
	//
	//    Edward Angel,
	//    Interactive Computer Graphics:
	//    A Top-Down Approach with OpenGL,
	//    Second Edition,
	//    Addison Wesley, 2000.

	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();						// moves to center of screen
	translateView ( projUpDown, projLeftRight );
	polarView( projDistance, projTwist, projElevation, projAzimuth );

	if(on_sphere && draw_sphere){
		drawSphere(0.98, 40, 40);
	}

	switch(drawing){
		case 0:
			draw_triangles();
			break;
		case 1:
			draw_cells();
			break;
		case 2:
			draw_edges();
			break;
	}

	switch(draw_lines){
		case 0:
			draw_triangle_lines();
			break;
		case 1:
			draw_cell_lines();
			break;
		case 2:
			draw_edge_lines();
			break;
	}

	switch(region_draw){
		case true:
			draw_regions();
			break;
		default:
			break;
	}

	//
	//  Clear all the buffers.
	//
	glFlush ( );
	//
	//  Switch between the two buffers for fast animation.
	//
	glutSwapBuffers ( );

	return;
}/*}}}*/
void mouse ( int btn, int state, int x, int y ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    MOUSE determines the response to mouse input.
	//
	//  Discussion:
	//
	//    The original routine assumed the user had a three button mouse, and
	//    dedicated one axis to each.
	//
	//    Since Apple prefers the esthetics of a one button mouse, we're forced
	//    to live with that choice.  This routine alternately pauses rotation,
	//    or increments the rotation axis by 1, no matter which button is pushed.
	//
	//  Modified:
	//
	//    30 December 2008
	//
	//  Author:
	//
	//    John Burkardt
	//
	//  Reference:
	//
	//    Edward Angel,
	//    Interactive Computer Graphics:
	//    A Top-Down Approach with OpenGL,
	//    Second Edition,
	//    Addison Wesley, 2000.
	//
	if ( ( btn == GLUT_LEFT_BUTTON   && state == GLUT_DOWN ) ||
			( btn == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN ) ||
			( btn == GLUT_RIGHT_BUTTON  && state == GLUT_DOWN ) ) {
	}

	axis = axis % 3;

	return;
}/*}}}*/
void gl_init ( ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    gl_init initializes OpenGL state variables dealing with viewing and attributes.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    01 January 2010
	//
	//  Author:
	//
	//    John Burkardt, Geoff Womeldorff, Doug Jacobsen
	//
	//  Reference:
	//
	//    Edward Angel,
	//    Interactive Computer Graphics:
	//    A Top-Down Approach with OpenGL,
	//    Second Edition,
	//    Addison Wesley, 2000.
	GLfloat line_width;
	GLfloat point_size;
	//
	//  Set the background to WHITE.
	//
	glClearColor ( 1.0, 1.0, 1.0, 1.0 );
	//
	//  Antialiasing.
	//
	glDepthFunc( GL_LESS );
	glEnable(GL_DEPTH_TEST);
	glBlendFunc ( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glHint ( GL_LINE_SMOOTH_HINT, GL_DONT_CARE );
	glShadeModel( GL_SMOOTH );

	glEnableClientState( GL_VERTEX_ARRAY );
	glEnableClientState( GL_COLOR_ARRAY );


	glEnable ( GL_POINT_SMOOTH );
	glEnable ( GL_LINE_SMOOTH );
	//
	//  The default point size is 1.0.
	//
	if ( ncells <= 400 ) {
		point_size = 16.0;
	} else if ( ncells <= 800 ) {
		point_size = 8.0;
	} else if ( ncells <= 1600 ) {
		point_size = 4.0;
	} else if ( ncells <= 3200 ) {
		point_size = 2.0;
	} else {
		point_size = 1.0;
	}

	point_size = 8.0;
	glPointSize ( point_size );
	//
	//  The default line width is 1.0.
	//
	if ( ncells <= 1600 ){
		line_width = 1.0;
	} else {
		line_width = 0.1;
	}
	glLineWidth ( line_width );

	return;
}/*}}}*/
void myReshape ( int w, int h ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    MYRESHAPE determines the window mapping.
	//
	//  Modified:
	//
	//    30 December 2008
	//
	//  Author:
	//
	//    John Burkardt, Geoff Womeldorff, Doug Jacobsen
	//
	//  Reference:
	//
	//    Edward Angel,
	//    Interactive Computer Graphics:
	//    A Top-Down Approach with OpenGL,
	//    Second Edition,
	//    Addison Wesley, 2000.
	//
	//
	
	if ( w <= h )
	{
		glOrtho (
				-1.05, +1.05,
				- 1.05 * ( GLfloat ) h / ( GLfloat ) w, + 1.05 * ( GLfloat ) h / ( GLfloat ) w,
				-10.0, 10.0 );
	}
	else
	{
		glOrtho (
				- 1.05 * ( GLfloat ) h / ( GLfloat ) w, + 1.05 * ( GLfloat ) h / ( GLfloat ) w,
				- 1.05, + 1.05,
				-10.0, 10.0 );
	}

	glViewport ( 0, 0, w, h );
	glMatrixMode ( GL_PROJECTION );
	glLoadIdentity ( );

	gluPerspective( 45.0f, (GLfloat)w / (GLfloat)h, 0.1f, 5.0f );
	//gluPerspective( 45.0f, (GLfloat)w / (GLfloat)h, (GLfloat)zNear, (GLfloat)zFar);
	//glMatrixMode ( GL_MODELVIEW );

	return;
}/*}}}*/
void single_screenshot ( ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    SINGLE_SCREENSHOT takes a single screenshot of the grid, and exits mpas_draw
	//
	//  Modified:
	//
	//    02/08/13
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//
	//color_mesh();
	glutPostRedisplay ( );
	
	if(single_ss){

		if(loop_it == 3){
			char filename[100];
			sprintf(filename, "ss.%04d.tga", cur_screenshot);
			screenshot(filename, 800, 800);
			exit(0);
		}
		loop_it++;
	}

	return;
}/*}}}*/
void timestamp ( ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    TIMESTAMP prints the current YMDHMS date as a time stamp.
	//
	//  Example:
	//
	//    31 May 2001 09:45:54 AM
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    08 July 2009
	//
	//  Author:
	//
	//    John Burkardt
	//
	//  Parameters:
	//
	//    None
	//
# define TIME_SIZE 40

	static char time_buffer[TIME_SIZE];
	const struct std::tm *tm_ptr;
	size_t len;
	std::time_t now;

	now = std::time ( NULL );
	tm_ptr = std::localtime ( &now );

	len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

	std::cout << time_buffer << "\n";

	return;
# undef TIME_SIZE
}/*}}}*/

void build_connectivity(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    BUILD_CONNECTIVITY builds the connectivity arrays for Display to draw
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    John Burkardt, Geoff Womeldorff, Doug Jacobsen
	//

	int i, j;
	int v1, v2;
	int c1, c2, c3, dc;
	int n_lilgons; // Smaller than 5 sides
	int n_pentagons; // 5 Sides
	int n_heptagons; // 7 Sides
	int n_bigagons; // Larger than 7 sides
	double xcell_min, xcell_max, ycell_min, ycell_max, zcell_min, zcell_max;
	double xvertex_min, xvertex_max, yvertex_min, yvertex_max, zvertex_min, zvertex_max;
	double xcell_range, ycell_range, zcell_range;
	double xvertex_range, yvertex_range, zvertex_range;
	double distance1, distance2, distance3;
	double distcell_max, distvertex_max;

	double max_distance, min_distance;
	bool keep_shape;


	//
	//  Get sizes.
	//
	ntime = netcdf_mpas_read_dim(filename, "Time" );
	nvertlevels = netcdf_mpas_read_dim(filename, "nVertLevels");
	ncells = netcdf_mpas_read_dim (filename, "nCells");
	nvertices = netcdf_mpas_read_dim(filename, "nVertices");
	nedges = netcdf_mpas_read_dim(filename, "nEdges");
	maxedges = netcdf_mpas_read_dim(filename, "maxEdges");
	on_sphere = netcdf_mpas_read_onsphere(filename);
	sphere_radius = netcdf_mpas_read_sphereradius(filename);
	draw_sphere = false;

	cell_cells.reserve((ncells+1)*maxedges*3*3);
	cell_lines.reserve((ncells+1)*maxedges*3*2);
	vertex_cells.reserve((nvertices+1)*3*3);
	vertex_lines.reserve((nvertices+1)*3*3*2);
	edge_cells.reserve((nedges+1)*4*3);
	edge_lines.reserve((nedges+1)*4*2*3);

	cout << "\n";
	if(!on_sphere){
		line_factor = 0.002;
		cout << "  Points are NOT on a sphere. " << endl;
	} else {
		cout << "  Points ARE on a sphere. " << endl;
	}

	cout << "  The number of time steps NTIME   = " << ntime << "\n";
	cout << "  The number of vertical levels NVERTLEVELS = " << nvertlevels << "\n";
	cout << "  The number of cells NCELLS       = " << ncells << "\n";
	cout << "  The number of vertices NVERTICES =  " << nvertices << "\n";
	cout << "  The number of edges NEDGES =  " << nedges << "\n";

	xcell = new double[ncells];
	ycell = new double[ncells];
	zcell = new double[ncells];

	xvertex = new double[nvertices];
	yvertex = new double[nvertices];
	zvertex = new double[nvertices];

	verticesoncell = new int[maxedges*ncells];
	verticesonedge = new int[2*nedges];
	cellsonvertex = new int[3*nvertices];
	cellsonedge = new int[2*nedges];
	nedgesoncell = new int[ncells];

	netcdf_mpas_read_xyzcell ( filename, ncells, xcell, ycell, zcell );
	netcdf_mpas_read_xyzvertex ( filename, nvertices, xvertex, yvertex, zvertex );

	rescale_cells_and_vertices();

	netcdf_mpas_read_verticesoncell ( filename, maxedges, ncells, verticesoncell );
	netcdf_mpas_read_verticesonedge ( filename, nedges, verticesonedge );
	netcdf_mpas_read_cellsonvertex ( filename, nvertices, cellsonvertex );
	netcdf_mpas_read_cellsonedge ( filename, nedges, cellsonedge );
	netcdf_mpas_read_nedgesoncell ( filename, ncells, nedgesoncell );

	n_lilgons = 0;
	n_pentagons = 0;
	n_heptagons = 0;
	n_bigagons = 0;

	r8vec_min_max(ncells, xcell, xcell_min, xcell_max, missing_value);
	r8vec_min_max(ncells, ycell, ycell_min, ycell_max, missing_value);
	r8vec_min_max(ncells, zcell, zcell_min, zcell_max, missing_value);
	r8vec_min_max(nvertices, xvertex, xvertex_min, xvertex_max, missing_value);
	r8vec_min_max(nvertices, yvertex, yvertex_min, yvertex_max, missing_value);
	r8vec_min_max(nvertices, zvertex, zvertex_min, zvertex_max, missing_value);

	xcell_range = (xcell_max - xcell_min);
	ycell_range = (ycell_max - ycell_min);
	zcell_range = (zcell_max - zcell_min);
	xvertex_range = (xvertex_max - xvertex_min);
	yvertex_range = (yvertex_max - yvertex_min);
	zvertex_range = (zvertex_max - zvertex_min);

	distcell_max = sqrt( pow(xcell_range/2.0,2) + pow(ycell_range/2.0, 2) + pow(zcell_range/2.0, 2))/2.0;
	distvertex_max = sqrt( pow(xvertex_range/2.0,2) + pow(yvertex_range/2.0, 2) + pow(zvertex_range/2.0, 2))/2.0;

	max_distance = 0.0;

	//
	//  connectivity is 1 based.  Fix that.
	//
	for ( i = 0; i < nvertices; i++ )
	{
		for ( j = 0; j < 3; j++ )
		{
			cellsonvertex[i*3+j]--;
		}
	}

	for( i = 0; i < ncells; i++){
		for ( j = 0; j < maxedges; j++){
			verticesoncell[i*maxedges + j]--;
		}
		if(nedgesoncell[i] < 5){
			n_lilgons++;
		} else if(nedgesoncell[i] == 5){
			n_pentagons++;
		} else if(nedgesoncell[i] == 7){
			n_heptagons++;
		} else if(nedgesoncell[i] > 7){
			n_bigagons++;
		}
	}

	cout << " Number of Lilgons (< 5 sides) " << n_lilgons << endl;
	cout << " Number of Pentagons " << n_pentagons << endl;
	cout << " Number of Heptagons " << n_heptagons << endl;
	cout << " Number of Bigagons (> 7 sides) " << n_bigagons << endl;

	for( i = 0; i < nedges; i++){
		for( j = 0; j < 2; j++){
			cellsonedge[i*2 + j]--;
			verticesonedge[i*2 + j]--;
		}
	}


	for(i = 0; i < nvertices; i++){
		c1 = cellsonvertex[i*3];
		c2 = cellsonvertex[i*3 + 1];
		c3 = cellsonvertex[i*3 + 2];

		keep_shape = true;

		if(c1 == -1){
			if(c2 == -1){
				dc = c3;
			} else {
				dc = c2;
			}
		} else {
			dc = c1;
		}

		if(c1 == -1){
			c1 = dc;
		}
		if(c2 == -1){
			c2 = dc;
		}
		if(c3 == -1){
			c3 = dc;
		}

		distance1 = sqrt( max( pow(xcell[c1] - xcell[c2], 2), max( pow(xcell[c2] - xcell[c3], 2), pow(xcell[c3] - xcell[c1], 2))));
		distance2 = sqrt( max( pow(ycell[c1] - ycell[c2], 2), max( pow(ycell[c2] - ycell[c3], 2), pow(ycell[c3] - ycell[c1], 2))));
		distance3 = sqrt( max( pow(zcell[c1] - zcell[c2], 2), max( pow(zcell[c2] - zcell[c3], 2), pow(zcell[c3] - zcell[c1], 2))));

		if(distance1 > xcell_range/2.0 || distance2 > ycell_range/2.0 || distance3 > zcell_range/2.0){
			keep_shape = false;
		}

		if(keep_shape || on_sphere){
			vertex_cells.push_back(xcell[c1]);
			vertex_cells.push_back(ycell[c1]);
			vertex_cells.push_back(zcell[c1]);

			vertex_cells.push_back(xcell[c2]);
			vertex_cells.push_back(ycell[c2]);
			vertex_cells.push_back(zcell[c2]);

			vertex_cells.push_back(xcell[c3]);
			vertex_cells.push_back(ycell[c3]);
			vertex_cells.push_back(zcell[c3]);

			if(on_sphere){
				vertex_lines.push_back(xcell[c1]*line_factor);
				vertex_lines.push_back(ycell[c1]*line_factor);
				vertex_lines.push_back(zcell[c1]*line_factor);
				vertex_lines.push_back(xcell[c2]*line_factor);
				vertex_lines.push_back(ycell[c2]*line_factor);
				vertex_lines.push_back(zcell[c2]*line_factor);

				vertex_lines.push_back(xcell[c2]*line_factor);
				vertex_lines.push_back(ycell[c2]*line_factor);
				vertex_lines.push_back(zcell[c2]*line_factor);
				vertex_lines.push_back(xcell[c3]*line_factor);
				vertex_lines.push_back(ycell[c3]*line_factor);
				vertex_lines.push_back(zcell[c3]*line_factor);

				vertex_lines.push_back(xcell[c3]*line_factor);
				vertex_lines.push_back(ycell[c3]*line_factor);
				vertex_lines.push_back(zcell[c3]*line_factor);
				vertex_lines.push_back(xcell[c1]*line_factor);
				vertex_lines.push_back(ycell[c1]*line_factor);
				vertex_lines.push_back(zcell[c1]*line_factor);
			} else {
				vertex_lines.push_back(xcell[c1]);
				vertex_lines.push_back(ycell[c1]);
				vertex_lines.push_back(line_factor);
				vertex_lines.push_back(xcell[c2]);
				vertex_lines.push_back(ycell[c2]);
				vertex_lines.push_back(line_factor);

				vertex_lines.push_back(xcell[c2]);
				vertex_lines.push_back(ycell[c2]);
				vertex_lines.push_back(line_factor);
				vertex_lines.push_back(xcell[c3]);
				vertex_lines.push_back(ycell[c3]);
				vertex_lines.push_back(line_factor);

				vertex_lines.push_back(xcell[c3]);
				vertex_lines.push_back(ycell[c3]);
				vertex_lines.push_back(line_factor);
				vertex_lines.push_back(xcell[c1]);
				vertex_lines.push_back(ycell[c1]);
				vertex_lines.push_back(line_factor);
			}
		} else {
			vertex_cells.push_back(0.0);
			vertex_cells.push_back(0.0);
			vertex_cells.push_back(0.0);
			vertex_cells.push_back(0.0);
			vertex_cells.push_back(0.0);
			vertex_cells.push_back(0.0);
			vertex_cells.push_back(0.0);
			vertex_cells.push_back(0.0);
			vertex_cells.push_back(0.0);
		}
	}

	delete(cellsonvertex);

	for(i = 0; i < ncells; i++){
		for(j = 0; j < nedgesoncell[i]; j++){
			v1 = verticesoncell[i*maxedges + j%nedgesoncell[i]];
			v2 = verticesoncell[i*maxedges + (j+1)%nedgesoncell[i]];

			keep_shape = true;

			distance1 = sqrt(max( pow(xvertex[v1] - xvertex[v2], 2), max( pow(xvertex[v1] - xcell[i], 2), pow(xvertex[v2] - xcell[i], 2))));
			distance2 = sqrt(max( pow(yvertex[v1] - yvertex[v2], 2), max( pow(yvertex[v1] - ycell[i], 2), pow(yvertex[v2] - ycell[i], 2))));
			distance3 = sqrt(max( pow(zvertex[v1] - zvertex[v2], 2), max( pow(zvertex[v1] - zcell[i], 2), pow(zvertex[v2] - zcell[i], 2))));

			if(distance1 > xvertex_range/2.0 || distance2 > yvertex_range/2.0 || distance3 > zvertex_range/2.0){
				keep_shape = false;
			}

			if(keep_shape || on_sphere){
				cell_cells.push_back(xcell[i]);
				cell_cells.push_back(ycell[i]);
				cell_cells.push_back(zcell[i]);
				cell_cells.push_back(xvertex[v1]);
				cell_cells.push_back(yvertex[v1]);
				cell_cells.push_back(zvertex[v1]);
				cell_cells.push_back(xvertex[v2]);
				cell_cells.push_back(yvertex[v2]);
				cell_cells.push_back(zvertex[v2]);

				if(on_sphere){
					cell_lines.push_back(xvertex[v1]*line_factor);
					cell_lines.push_back(yvertex[v1]*line_factor);
					cell_lines.push_back(zvertex[v1]*line_factor);

					cell_lines.push_back(xvertex[v2]*line_factor);
					cell_lines.push_back(yvertex[v2]*line_factor);
					cell_lines.push_back(zvertex[v2]*line_factor);
				} else {
					cell_lines.push_back(xvertex[v1]);
					cell_lines.push_back(yvertex[v1]);
					cell_lines.push_back(line_factor);

					cell_lines.push_back(xvertex[v2]);
					cell_lines.push_back(yvertex[v2]);
					cell_lines.push_back(line_factor);
				}
			} else {
				cell_cells.push_back(0.0);
				cell_cells.push_back(0.0);
				cell_cells.push_back(0.0);
				cell_cells.push_back(0.0);
				cell_cells.push_back(0.0);
				cell_cells.push_back(0.0);
				cell_cells.push_back(0.0);
				cell_cells.push_back(0.0);
				cell_cells.push_back(0.0);
			}
		}
	}

	delete(verticesoncell);

	for(i = 0; i < nedges; i++){
		c1 = cellsonedge[i*2];
		c2 = cellsonedge[i*2+1];
		v1 = verticesonedge[i*2];
		v2 = verticesonedge[i*2+1];

		keep_shape = true;

		if(c1 == -1){
			c1 = c2;
		} else if(c2 == -1){
			c2 = c1;
		}

		if(v1 == -1){
			v1 = v2;
		} else if (v2 == -1){
			v2 = v1;
		}

		distance1 = pow(xcell[c1] - xvertex[v1], 2);
		distance2 = pow(ycell[c1] - yvertex[v1], 2);
		distance3 = pow(zcell[c1] - zvertex[v1], 2);

		if(c1 != c2){
			distance1 = max(distance1, pow(xcell[c2] - xvertex[v1], 2));
			distance2 = max(distance2, pow(ycell[c2] - yvertex[v1], 2));
			distance3 = max(distance3, pow(zcell[c2] - zvertex[v1], 2));

			distance1 = max(distance1, pow(xcell[c2] - xvertex[v2], 2));
			distance2 = max(distance2, pow(ycell[c2] - yvertex[v2], 2));
			distance3 = max(distance3, pow(zcell[c2] - zvertex[v2], 2));

			distance1 = max(distance1, pow(xcell[c1] - xvertex[v2], 2));
			distance2 = max(distance2, pow(ycell[c1] - yvertex[v2], 2));
			distance3 = max(distance3, pow(zcell[c1] - zvertex[v2], 2));
		} else {
			distance1 = max(distance1, pow(xvertex[v1] - xvertex[v2], 2));
			distance2 = max(distance2, pow(yvertex[v1] - yvertex[v2], 2));
			distance3 = max(distance3, pow(zvertex[v1] - zvertex[v2], 2));

			distance1 = max(distance1, pow(xcell[c1] - xvertex[v2], 2));
			distance2 = max(distance2, pow(ycell[c1] - yvertex[v2], 2));
			distance3 = max(distance3, pow(zcell[c1] - zvertex[v2], 2));
		}

		if(distance1 > xvertex_range/2.0 || distance2 > yvertex_range/2.0 || distance3 > zvertex_range/2.0){
			keep_shape = false;
		}

		if(keep_shape || on_sphere){ 
			edge_cells.push_back(xcell[c1]);	
			edge_cells.push_back(ycell[c1]);	
			edge_cells.push_back(zcell[c1]);	

			edge_cells.push_back(xvertex[v1]);
			edge_cells.push_back(yvertex[v1]);
			edge_cells.push_back(zvertex[v1]);

			edge_cells.push_back(xvertex[v2]);
			edge_cells.push_back(yvertex[v2]);
			edge_cells.push_back(zvertex[v2]);

			edge_cells.push_back(xcell[c2]);	
			edge_cells.push_back(ycell[c2]);	
			edge_cells.push_back(zcell[c2]);	

			edge_cells.push_back(xvertex[v2]);
			edge_cells.push_back(yvertex[v2]);
			edge_cells.push_back(zvertex[v2]);

			edge_cells.push_back(xvertex[v1]);
			edge_cells.push_back(yvertex[v1]);
			edge_cells.push_back(zvertex[v1]);

			if(on_sphere){
				edge_lines.push_back(xcell[c1]*line_factor);
				edge_lines.push_back(ycell[c1]*line_factor);
				edge_lines.push_back(zcell[c1]*line_factor);
				edge_lines.push_back(xvertex[v1]*line_factor);
				edge_lines.push_back(yvertex[v1]*line_factor);
				edge_lines.push_back(zvertex[v1]*line_factor);

				edge_lines.push_back(xvertex[v1]*line_factor);
				edge_lines.push_back(yvertex[v1]*line_factor);
				edge_lines.push_back(zvertex[v1]*line_factor);
				if(c1 == c2){
					edge_lines.push_back(xvertex[v2]*line_factor);
					edge_lines.push_back(yvertex[v2]*line_factor);
					edge_lines.push_back(zvertex[v2]*line_factor);
				} else {
					edge_lines.push_back(xcell[c2]*line_factor);
					edge_lines.push_back(ycell[c2]*line_factor);
					edge_lines.push_back(zcell[c2]*line_factor);

					edge_lines.push_back(xcell[c2]*line_factor);
					edge_lines.push_back(ycell[c2]*line_factor);
					edge_lines.push_back(zcell[c2]*line_factor);
					edge_lines.push_back(xvertex[v2]*line_factor);
					edge_lines.push_back(yvertex[v2]*line_factor);
					edge_lines.push_back(zvertex[v2]*line_factor);
				}

				edge_lines.push_back(xvertex[v2]*line_factor);
				edge_lines.push_back(yvertex[v2]*line_factor);
				edge_lines.push_back(zvertex[v2]*line_factor);
				edge_lines.push_back(xcell[c1]*line_factor);
				edge_lines.push_back(ycell[c1]*line_factor);
				edge_lines.push_back(zcell[c1]*line_factor);
			} else {
				edge_lines.push_back(xcell[c1]);
				edge_lines.push_back(ycell[c1]);
				edge_lines.push_back(line_factor);
				edge_lines.push_back(xvertex[v1]);
				edge_lines.push_back(yvertex[v1]);
				edge_lines.push_back(line_factor);

				edge_lines.push_back(xvertex[v1]);
				edge_lines.push_back(yvertex[v1]);
				edge_lines.push_back(line_factor);
				if(c1 == c2){
					edge_lines.push_back(xvertex[v2]);
					edge_lines.push_back(yvertex[v2]);
					edge_lines.push_back(line_factor);
				} else {
					edge_lines.push_back(xcell[c2]);
					edge_lines.push_back(ycell[c2]);
					edge_lines.push_back(line_factor);

					edge_lines.push_back(xcell[c2]);
					edge_lines.push_back(ycell[c2]);
					edge_lines.push_back(line_factor);
					edge_lines.push_back(xvertex[v2]);
					edge_lines.push_back(yvertex[v2]);
					edge_lines.push_back(line_factor);
				}

				edge_lines.push_back(xvertex[v2]);
				edge_lines.push_back(yvertex[v2]);
				edge_lines.push_back(line_factor);
				edge_lines.push_back(xcell[c1]);
				edge_lines.push_back(ycell[c1]);
				edge_lines.push_back(line_factor);
			}
		} else {
			edge_cells.push_back(0.0);
			edge_cells.push_back(0.0);
			edge_cells.push_back(0.0);

			edge_cells.push_back(0.0);
			edge_cells.push_back(0.0);
			edge_cells.push_back(0.0);

			edge_cells.push_back(0.0);
			edge_cells.push_back(0.0);
			edge_cells.push_back(0.0);

			edge_cells.push_back(0.0);
			edge_cells.push_back(0.0);
			edge_cells.push_back(0.0);

			edge_cells.push_back(0.0);
			edge_cells.push_back(0.0);
			edge_cells.push_back(0.0);

			edge_cells.push_back(0.0);
			edge_cells.push_back(0.0);
			edge_cells.push_back(0.0);
		}
	}

	delete(verticesonedge);
	delete(cellsonedge);

	cell_values = new double[ncells];
	triangle_values = new double[nvertices];
	edge_values = new double[nedges];

	color_cells();
	color_triangles();
	color_edges();

	delete [] xcell;
	delete [] ycell;
	delete [] zcell;
	delete [] xvertex;
	delete [] yvertex;
	delete [] zvertex;

}/*}}}*/
void setup_ranges(){/*{{{*/
	int num_vars;

	num_vars = netcdf_mpas_read_num_vars(filename);

	ranges.resize(num_vars);

	return;
}/*}}}*/
void build_range(int id){/*{{{*/
	int num_items;
	double max, min;
	vector<double> temp_data;
	bool failed;

	if(color_bar == 2){
		if(ranges[id].size() == 0){
			cout << "Building min-max range for coloring." << endl;
			num_items = netcdf_mpas_field_num_items(filename, id);

			temp_data.resize(num_items);

			netcdf_mpas_read_full_field(filename, id, &temp_data[0]);

			r8vec_min_max(num_items, &temp_data[0], min, max, missing_value);

			ranges[id].push_back(min);
			ranges[id].push_back(max);
			cout << "Range build complete." << endl;
		}
	} else if(color_bar == 1){
		hard_ranges.clear();
		cout << endl;
		cout << "Input minimum for color bar:" << endl;

		do{
			cin >> min;
			if(cin.fail()){
				failed = true;
				cin.clear();
				cin.ignore(1024,'\n');
				cout << "Invalid input or the input buffer needed to be cleared. Please try again." << endl;
			} else {
				failed = false;
			}
		}while(failed);

		cout << "Input maximum for color bar:" << endl;
		do{
			cin >> max;
			if(cin.fail()){
				failed = true;
				cin.clear();
				cin.ignore(1024,'\n');
				cout << "Invalid input or the input buffer needed to be cleared. Please try again." << endl;
			} else {
				failed = false;
			}
		}while(failed);

		hard_ranges.push_back(min);
		hard_ranges.push_back(max);
	} else {
		switch(drawing){
			case 0:
				num_items = nvertices;
				break;
			case 1:
				num_items = ncells;
				break;
			case 2:
				num_items = nedges;
				break;
			default:
				return;
				break;
		}

		temp_data.resize(num_items);

		netcdf_mpas_read_field(filename, id, &temp_data.at(0), cur_time, cur_level);
	}
}/*}}}*/

void rescale_cells_and_vertices(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    RESCALE_CELLS_AND_VERTICES scales xcell, ycell, and zcell arrays to have a norm of 1
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    John Burkardt, Doug Jacobsen
	//


	int i;
	int cell, vertex;
	double norm;
	xyz_min[0] = r8vec_min ( ncells, xcell, missing_value);
	xyz_max[0] = r8vec_max ( ncells, xcell, missing_value);
	xyz_min[1] = r8vec_min ( ncells, ycell, missing_value);
	xyz_max[1] = r8vec_max ( ncells, ycell, missing_value);
	xyz_min[2] = r8vec_min ( ncells, zcell, missing_value);
	xyz_max[2] = r8vec_max ( ncells, zcell, missing_value);

	xyz_range[0] = xyz_max[0] - xyz_min[0];
	xyz_range[1] = xyz_max[1] - xyz_min[1];
	xyz_range[2] = xyz_max[2] - xyz_min[2];

	if ( xyz_range[0] == 0.0 ){
		cout << "\n";
		cout << "rescale_cells_and_vertices(): - Fatal error!\n";
		cout << "  The X data range is 0.\n";
		exit ( 1 );
	}

	if ( xyz_range[1] == 0.0 ){
		cout << "\n";
		cout << "rescale_cells_and_vertices(): - Fatal error!\n";
		cout << "  The Y data range is 0.\n";
		exit ( 1 );
	}

	if(on_sphere){
		if ( xyz_range[2] == 0.0 ){
			cout << "\n";
			cout << "rescale_cells_and_vertices(): - Fatal error!\n";
			cout << "  The Z data range is 0.\n";
			exit ( 1 );
		}
	}

	xyz_scale = 0.0;
	for (i = 0; i < 3; i++ )
	{
		xyz_center[i] = ( xyz_min[i] + xyz_max[i] ) / 2.0;
		xyz_scale = std::max ( xyz_scale, ( xyz_max[i] - xyz_min[i] ) / 2.0 );
	}

	for ( cell = 0; cell < ncells; cell++ )
	{
		norm = xyz_scale;
		xcell[cell] = (xcell[cell] - xyz_center[0]) / norm;
		ycell[cell] = (ycell[cell] - xyz_center[1]) / norm;
		zcell[cell] = (zcell[cell] - xyz_center[2]) / norm;
	}
/*
	// Do Vertices
	xyz_min[0] = r8vec_min ( nvertices, xvertex, missing_value);
	xyz_max[0] = r8vec_max ( nvertices, xvertex, missing_value);
	xyz_min[1] = r8vec_min ( nvertices, yvertex, missing_value);
	xyz_max[1] = r8vec_max ( nvertices, yvertex, missing_value);
	xyz_min[2] = r8vec_min ( nvertices, zvertex, missing_value);
	xyz_max[2] = r8vec_max ( nvertices, zvertex, missing_value);

	xyz_range[0] = xyz_max[0] - xyz_min[0];
	xyz_range[1] = xyz_max[1] - xyz_min[1];
	xyz_range[2] = xyz_max[2] - xyz_min[2];

	if ( xyz_range[0] == 0.0 ){
		cout << "\n";
		cout << "rescale_cells_and_vertices(): - Fatal error!\n";
		cout << "  The X data range is 0.\n";
		exit ( 1 );
	}

	if ( xyz_range[1] == 0.0 ){
		cout << "\n";
		cout << "rescale_cells_and_vertices(): - Fatal error!\n";
		cout << "  The Y data range is 0.\n";
		exit ( 1 );
	}

	if(on_sphere){
		if ( xyz_range[2] == 0.0 ){
			cout << "\n";
			cout << "rescale_cells_and_vertices(): - Fatal error!\n";
			cout << "  The Z data range is 0.\n";
			exit ( 1 );
		}
	}

	xyz_scale = 0.0;
	for (i = 0; i < 3; i++ )
	{
		xyz_center[i] = ( xyz_min[i] + xyz_max[i] ) / 2.0;
		xyz_scale = std::max ( xyz_scale, ( xyz_max[i] - xyz_min[i] ) / 2.0 );
	}
	*/

	for ( vertex = 0; vertex < nvertices; vertex++ )
	{
		norm = xyz_scale;
		xvertex[vertex] = (xvertex[vertex] - xyz_center[0]) / norm;
		yvertex[vertex] = (yvertex[vertex] - xyz_center[1]) / norm;
		zvertex[vertex] = (zvertex[vertex] - xyz_center[2]) / norm;
	}

}/*}}}*/

void draw_triangles ( ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    DRAW_TRIANGLES draws the Delaunay triangles
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    Geoff Womeldorff, Doug Jacobsen
	//

	glColorPointer( 3, GL_FLOAT, 0, &vertex_colors[0] );
	glVertexPointer( 3, GL_FLOAT, 0, &vertex_cells[0] );
	glDrawArrays( GL_TRIANGLES, 0, vertex_cells.size()/3);

	return;
}/*}}}*/
void draw_cells ( ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    DRAW_CELLS draws the Voronoi cells
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    Geoff Womeldorff, Doug Jacobsen
	//
	
	glColorPointer( 3, GL_FLOAT, 0, &cell_colors[0] );
	glVertexPointer( 3, GL_FLOAT, 0, &cell_cells[0] );
	glDrawArrays( GL_TRIANGLES, 0, cell_cells.size()/3);

	return;
}/*}}}*/
void draw_edges ( ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    DRAW_EDGES draws quadrilaterals to represent edge values
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    Geoff Womeldorff, Doug Jacobsen
	//
	
	glColorPointer( 3, GL_FLOAT, 0, &edge_colors[0] );
	glVertexPointer( 3, GL_FLOAT, 0, &edge_cells[0] );
	glDrawArrays( GL_TRIANGLES, 0, edge_cells.size()/3);

	return;
}/*}}}*/

void draw_triangle_lines(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    DRAW_TRIANGLE_LINES draws the lines for the Delaunay triangle edges
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    Geoff Womeldorff, Doug Jacobsen
	//

	glDisableClientState( GL_COLOR_ARRAY );
	glColor3f ( LINE_R, LINE_G, LINE_B );
	glVertexPointer( 3, GL_FLOAT, 0, &vertex_lines[0] );
	glDrawArrays( GL_LINES, 0, vertex_lines.size()/3 );
	glEnableClientState( GL_COLOR_ARRAY);

	return;
}/*}}}*/
void draw_cell_lines(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    DRAW_CELL_LINES draws the lines for the Voronoi cell edges.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    Geoff Womeldorff, Doug Jacobsen
	//

	glDisableClientState( GL_COLOR_ARRAY );
	glColor3f ( LINE_R, LINE_G, LINE_B );
	glVertexPointer( 3, GL_FLOAT, 0, &cell_lines[0] );
	glDrawArrays( GL_LINES, 0, cell_lines.size()/3 );
	glEnableClientState( GL_COLOR_ARRAY);

	return;
}/*}}}*/
void draw_edge_lines(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    DRAW_EDGE_LINES draws the lines for the edge "edges".
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    Geoff Womeldorff, Doug Jacobsen
	//

	glDisableClientState( GL_COLOR_ARRAY );
	glColor3f ( LINE_R, LINE_G, LINE_B );
	glVertexPointer( 3, GL_FLOAT, 0, &edge_lines[0] );
	glDrawArrays( GL_LINES, 0, edge_lines.size()/3 );
	glEnableClientState( GL_COLOR_ARRAY);

	return;
}/*}}}*/

void build_regions(){/*{{{*/
	ifstream reg_in("RegionCenters");
	ifstream rad_in("RegionRadii");
	double t;
	double x, y, z, r;
	double scale;
	double planar_r;
	int i, j;

	if(!reg_in){
		cout << "File RegionCenters not found. Skipping regions." << endl;
		return;
	} else {
		while(!reg_in.eof()){
			reg_in >> x;
			reg_in >> y;
			reg_in >> z;

			if(reg_in.good()){
				region_centers.push_back(x);
				region_centers.push_back(y);
				region_centers.push_back(z);
			}
		}
	}

	reg_in.close();

	if(!rad_in){
		cout << "File RegionRadii not found. Skipping regions." << endl;
		return;
	} else {
		while(!rad_in.eof()){
			rad_in >> r;

			if(rad_in.good()){
				region_radii.push_back(r);
			}
		}
	}

	rad_in.close();

	for(i = 0; i < region_radii.size(); i++){
		x = region_centers.at(i*3);
		y = region_centers.at(i*3+1);
		z = region_centers.at(i*3+2);

		scale = sin(region_radii.at(i) + M_PI/2.0);	
		
		region_lines.push_back(x*scale*region_line_factor);
		region_lines.push_back(y*scale*region_line_factor);
		region_lines.push_back(z*scale*region_line_factor);

		region_radii.at(i) = sin(region_radii.at(i));

		region_angles.push_back(atan2(region_centers.at(i*3+1),region_centers.at(i*3))*180.0/M_PI);
		region_angles.push_back((asin(region_centers.at(i*3+2))+M_PI/2.0)*180.0/M_PI);
	}
	

	for(i = 0; i < region_centers.size(); i++){
		region_centers.at(i) = region_centers.at(i)*region_center_factor;
	}

	regions_built = true;
}/*}}}*/
void draw_regions(){/*{{{*/
	double theta, phi;
	float h, s, v;
	float r, g, b;
	glDisableClientState( GL_COLOR_ARRAY );

	s = 0.8;
	v = 0.0;

	for(int i = 0; i < region_radii.size(); i++){
		phi = region_angles.at(i*2);
		theta = region_angles.at(i*2+1);
		h = (1.0*i/region_radii.size()) * 0.8;

		hsv_to_rgb(h,s,v,r,g,b);

		glPushMatrix();
			glColor3f( r, g, b);
			glTranslated(region_lines.at(i*3), region_lines.at(i*3+1), region_lines.at(i*3+2));

			glRotated(phi, 0.0, 0.0, 1.0);
			glRotated(-theta, 0.0, 1.0, 0.0);

			gluCylinder(gluNewQuadric(), region_radii.at(i)*1.015, region_radii.at(i)*1.015, 0.015, 360, 5);
			gluDisk(gluNewQuadric(), region_radii.at(i), region_radii.at(i)*1.015, 360, 5);
		glPopMatrix();
	}

	glColor3f( REG_R, REG_G, REG_B );
	glVertexPointer( 3, GL_FLOAT, 0, &region_centers[0] );
	glDrawArrays( GL_POINTS, 0, region_centers.size()/3.0 );

	glEnableClientState( GL_COLOR_ARRAY );
}/*}}}*/

void color_mesh(){/*{{{*/
	switch(drawing){
		case 0:
			color_triangles();
			break;
		case 1:
			color_cells();
			break;
		case 2:
			color_edges();
			break;
		default:
			break;
	}
	return;
}/*}}}*/
void color_cells(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    COLOR_CELLS builds the color array used to display the Voronoi cells
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    Doug Jacobsen
	//

	int i, j, k, o;
	double max, min;
	long *dims;
	int num_dims = 0;
	int cell_dim;
	int vert_dim;
	int time_dim;
	int num_items;
	float h, s, v;
	float r, g, b;

	cell_colors.clear();

	s = 1.0;
	v = 1.0;
	if(cell_field == -1){
		for(i = 0; i < ncells; i++){
			o = nedgesoncell[i];

			for(j = 0; j < 3*o; j++){
				cell_colors.push_back(0.8);
				cell_colors.push_back(0.8);
				cell_colors.push_back(0.8);
			} 
		}
	} else {
		netcdf_mpas_read_field(filename, cell_field, cell_values, cur_time, cur_level);

		if(color_bar == 2){
			min = ranges[cell_field].at(0);
			max = ranges[cell_field].at(1);
		} else if (color_bar == 1){
			min = hard_ranges.at(0);
			max = hard_ranges.at(1);
		} else {
			r8vec_min_max(ncells, cell_values, min, max, missing_value);
		}

		cout << "Min: " << min << " Max: " << max << endl;

		for(i = 0; i < ncells; i++){

			o = nedgesoncell[i];

			if(cell_values[i] == missing_value){
				r = 0.5;
				g = 0.5;
				b = 0.5;
			} else {
				if((max-min) != 0.0){
					if(cell_values[i] >= max){
						h = range_factor;
					} else if (cell_values[i] <= min){
						h = 0.0;
					} else {
						h = (cell_values[i] - min)/(max-min) * range_factor;
					}
				} else {
					h = (cell_values[i] - min) * range_factor;
				}

				hsv_to_rgb(h, s, v, b, g, r);
			}

			for(j = 0; j < o*3; j++){
				cell_colors.push_back(r);
				cell_colors.push_back(g);
				cell_colors.push_back(b);
			}

		}
	}

	return;
}/*}}}*/
void color_triangles(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    COLOR_TRIANGLES builds the color array used to display the Delaunay triangles
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    Doug Jacobsen
	//

	int i, j;
	double max, min;
	float h, s, v;
	float r, g, b;

	vertex_colors.clear();

	s = 1.0;
	v = 1.0;
	if(vertex_field == -1){
		for(i = 0; i < nvertices; i++){
			for(j = 0; j < 3; j++){
				vertex_colors.push_back(0.8);
				vertex_colors.push_back(0.8);
				vertex_colors.push_back(0.8);
			}
		}
	} else {
		netcdf_mpas_read_field(filename, vertex_field, triangle_values, cur_time, cur_level);

		if(color_bar == 2){
			max = ranges[vertex_field][0];
			min = ranges[vertex_field][1];
		} else if(color_bar == 1) {
			min = hard_ranges.at(0);
			max = hard_ranges.at(1);
		} else {
			r8vec_min_max(nvertices, triangle_values, min, max, missing_value);
		}

		cout << "Min: " << min << " Max: " << max << endl;

		for(i = 0; i < nvertices; i++){
			if(triangle_values[i] == missing_value){
				r = 0.5;
				g = 0.5;
				b = 0.5;
			} else { 
				if(max-min != 0.0){
					if(triangle_values[i] >= max){
						h = range_factor;
					} else if (triangle_values[i] <= min){
						h = 0.0;
					} else {
						h = (triangle_values[i] - min)/(max-min)*range_factor;
					}
				}else{
					h = (triangle_values[i] - min)/1.0 * range_factor;
				}

				hsv_to_rgb(h, s, v, b, g, r);
			}

			for(j = 0; j < 3; j++){
				vertex_colors.push_back(r);
				vertex_colors.push_back(g);
				vertex_colors.push_back(b);
			}
		}
	}

	return;
}/*}}}*/
void color_edges(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    COLOR_EDGES builds the color array used to display the edge quads
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	int i, j, o;
	double max, min;
	float h, s, v;
	float r, g, b;

	edge_colors.clear();

	s = 1.0;
	v = 1.0;
	if(edge_field == -1){
		for(i = 0; i < nedges; i++){
			for(j = 0; j < 6; j++){
				edge_colors.push_back(0.8);
				edge_colors.push_back(0.8);
				edge_colors.push_back(0.8);
			}
		}
	} else {
		netcdf_mpas_read_field(filename, edge_field, edge_values, cur_time, cur_level);

		if(color_bar == 2){
			min = ranges[edge_field][0];
			max = ranges[edge_field][1];
		} else if(color_bar == 1) {
			min = hard_ranges.at(0);
			max = hard_ranges.at(1);
		} else {
			r8vec_min_max(nedges, edge_values, min, max, missing_value);
		}

		cout << "Min: " << min << " Max: " << max << endl;

		for(i = 0; i < nedges; i++){
			if(edge_values[i] == missing_value){
				r = 0.5;
				g = 0.5;
				b = 0.5;
			} else {
				if(max-min != 0.0){
					if(edge_values[i] >= max){
						h = range_factor;
					} else if (edge_values[i] <= min){
						h = 0.0;
					} else {
						h = (edge_values[i] - min)/(max-min)*range_factor;
					}
				} else {
					h = (edge_values[i] - min)/1.0 * range_factor;
				}

				hsv_to_rgb(h, s, v, b, g, r);
			}

			for(j = 0; j < 6; j++){
				edge_colors.push_back(r);
				edge_colors.push_back(g);
				edge_colors.push_back(b);
			}
		}
	}
}/*}}}*/

void arrowKeys( int a_keys, int x, int y ) {/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    ARROWKEYS processes special keys, for example the arrow keys, and the f-keys
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    Geoff Womeldorff
	//

	switch( a_keys ) {
		case GLUT_KEY_UP:
			projElevation = (double)( ( (int)projElevation + 1 ) % 360 );
			break;

		case GLUT_KEY_DOWN:
			projElevation = (double) ( ( (int)projElevation - 1 ) % 360 );
			break;

		case GLUT_KEY_RIGHT:
			projAzimuth = (double) ( ( (int)projAzimuth + 1 ) % 360 );
			break;

		case GLUT_KEY_LEFT:
			projAzimuth = (double) ( ( (int)projAzimuth - 1 ) % 360 );
			break;

		case GLUT_KEY_F1:
			glutFullScreen();
			break;

		case GLUT_KEY_F2:
			glutReshapeWindow( kWindowWidth, kWindowHeight );
			break;

		default:
			break;

	}

	return;
}/*}}}*/
void keyPressed( unsigned char key, int x, int y ) {/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    KEYPRESSED processes non-special keys, for example the letter keys
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    Geoff Womeldorff, Doug Jacobsen
	//


	usleep( 100 );

	// Used to get codes of keys you don't know.
//	cout << "Key pressed: " << key << endl;
//	cout << "Code is: " << (int)key << endl;

	switch( key ) {

		case KEY_ESCAPE:
			glutDestroyWindow( window );
			glDisableClientState( GL_VERTEX_ARRAY );
			glDisableClientState( GL_COLOR_ARRAY );
			//glDisableClientState( GL_NORMAL_ARRAY );
			delete [] nedgesoncell;

			cout << "End of normal execution. Exiting." << endl;

			exit( 0 );
		case KEY_l:
			cur_level = (cur_level + 1) % nvertlevels;
			cout << "Current vertical level: " << cur_level+1 << " out of " << nvertlevels << endl;
			color_mesh();
			break;
		case KEY_L:
			cur_level = (cur_level - 1) % nvertlevels;
			if(cur_level < 0)
				cur_level = nvertlevels - 1;
			cout << "Current vertical level: " << cur_level+1 << " out of " << nvertlevels << endl;
			color_mesh();
			break;
		case KEY_t:
			cur_time = (cur_time + 1) % ntime;
			cout << "Current time level: " << cur_time+1 << " out of " << ntime << endl;
			color_mesh();
			break;
		case KEY_T:
			cur_time = (cur_time - 1) % ntime;
			if(cur_time < 0)
				cur_time = ntime-1;
			cout << "Current time level: " << cur_time+1 << " out of " << ntime << endl;
			color_mesh();
			break;
		case KEY_s:
			projLeftRight -= 0.1;
			break;
		case KEY_f:
			projLeftRight += 0.1;
			break;
		case KEY_d:
			projUpDown -= 0.1;
			break;
		case KEY_e:
			projUpDown += 0.1;
			break;
		case KEY_S:
			projAzimuth = (double) ( ( (int)projAzimuth - 3 ) % 360 );
			break;
		case KEY_F:
			projAzimuth = (double) ( ( (int)projAzimuth + 3 ) % 360 );
			break;
		case KEY_D:
			projElevation = (double) ( ( (int)projElevation - 3 ) % 360 );
			break;
		case KEY_E:
			projElevation = (double)( ( (int)projElevation + 3 ) % 360 );
			break;
		case KEY_B:
			if(color_bar != 2){
				color_bar = 2;
			} else {
				color_bar = 0;
			}

			switch(color_bar){
				case 0:
					cout << "Coloring based on current time/vertlevel slice" << endl;
					break;
				case 2:
					cout << "Coloring base on entire time range" << endl;
					break;
				default:
					break;
			}
			if(color_bar == 2){
				switch(drawing){
					case 0:
						build_range(vertex_field);
						cout << " Range of values: Min = " << ranges.at(vertex_field).at(0) << ", Max = " << ranges.at(vertex_field).at(1) << endl;
						break;
					case 1:
						build_range(cell_field);
						cout << " Range of values: Min = " << ranges.at(cell_field).at(0) << ", Max = " << ranges.at(cell_field).at(1) << endl;
						break;
					case 2:
						build_range(edge_field);
						cout << " Range of values: Min = " << ranges.at(edge_field).at(0) << ", Max = " << ranges.at(edge_field).at(1) << endl;
						break;
					default:
						break;
				}
			}

			break;
		case KEY_b:
			color_bar = (color_bar + 1) % 2;
			switch(color_bar){
				case 0:
					cout << "Coloring based on current time/vertlevel slice" << endl;
					break;
				case 1:
					cout << "Coloring base on hard coded ranges" << endl;
					break;
				default:
					break;
			}
			if(color_bar == 1){
				switch(drawing){
					case 0:
						build_range(vertex_field);
						cout << " Range of values: Min = " << hard_ranges.at(0) << ", Max = " << hard_ranges.at(1) << endl;
						break;
					case 1:
						build_range(cell_field);
						cout << " Range of values: Min = " << hard_ranges.at(0) << ", Max = " << hard_ranges.at(1) << endl;
						break;
					case 2:
						build_range(edge_field);
						cout << " Range of values: Min = " << hard_ranges.at(0) << ", Max = " << hard_ranges.at(1) << endl;
						break;
					default:
						break;
				}
			}

			break;
		case KEY_w:
			draw_lines = (draw_lines + 1) % 4;
			break;
		case KEY_c:
			drawing = ( drawing + 1 ) % 3;
			break;
		case KEY_COMMA:
			projDistance += 0.05;
			zFar += 0.05;
			break;
		case KEY_PERIOD:
			projDistance -= 0.05;
			zFar -= 0.05;
			break;
		case KEY_v:
			switch(drawing){
				case 0:
					vertex_field = netcdf_mpas_list_nvertex_fields(filename);
					netcdf_mpas_print_field_info(filename, vertex_field);
					if(vertex_field >= 0 && color_bar == 2){
						build_range(vertex_field);
						cout << " Range of values: Min = " << ranges.at(vertex_field).at(0) << ", Max = " << ranges.at(vertex_field).at(1) << endl;
					}
					break;
				case 1:
					cell_field = netcdf_mpas_list_ncell_fields(filename);
					netcdf_mpas_print_field_info(filename, cell_field);
					if(cell_field >= 0 && color_bar == 2){
						build_range(cell_field);
						cout << " Range of values: Min = " << ranges.at(cell_field).at(0) << ", Max = " << ranges.at(cell_field).at(1) << endl;
					}
					break;
				case 2:
					edge_field = netcdf_mpas_list_nedge_fields(filename);
					netcdf_mpas_print_field_info(filename, edge_field);
					if(edge_field >= 0 && color_bar == 2){
						build_range(edge_field);
						cout << " Range of values: Min = " << ranges.at(edge_field).at(0) << ", Max = " << ranges.at(edge_field).at(1) << endl;
					}
					break;
				default:
					break;
			}
			color_mesh();
			break;
		case KEY_V:
			cout << "playing video." << endl;
			cur_time = 0;
			for(int i = 0; i < ntime; i++){
				cur_time = i;
				color_mesh();
				display();
				usleep(3.3e4);
			}
			break;
		case KEY_R:
			cout << " Trying regions. " << endl;
			if(!regions_built){
				build_regions();
			}

			if(regions_built){
				region_draw = !region_draw;
			}
			break;
		case KEY_r:
			cout << "Resetting to default parameters." << endl;
			cur_time = 0;
			cur_level = 0;
			cell_field = 0;
			edge_field = 0;
			vertex_field = 0;
			color_bar = 0;
			projDistance = 3.0;
			projUpDown = 0.0;
			projLeftRight = 0.0;
			projTwist = 0.0;
			projElevation = 0.0;
			projAzimuth = 0.0;

			color_mesh();
			break;
		case KEY_0:
			range_factor = 1.0;
			color_mesh();
			break;
		case KEY_1:
			range_factor = 0.1;
			color_mesh();
			break;
		case KEY_2:
			range_factor = 0.2;
			color_mesh();
			break;
		case KEY_3:
			range_factor = 0.3;
			color_mesh();
			break;
		case KEY_4:
			range_factor = 0.4;
			color_mesh();
			break;
		case KEY_5:
			range_factor = 0.5;
			color_mesh();
			break;
		case KEY_6:
			range_factor = 0.6;
			color_mesh();
			break;
		case KEY_7:
			range_factor = 0.7;
			color_mesh();
			break;
		case KEY_8:
			range_factor = 0.8;
			color_mesh();
			break;
		case KEY_9:
			range_factor = 0.9;
			color_mesh();
			break;
		case KEY_O:
			draw_sphere = !draw_sphere;
			break;
		case KEY_q:
			char filename[100];
			sprintf(filename, "ss.%04d.tga", cur_screenshot);
			screenshot(filename, 800, 800);
			cur_screenshot++;
			break;
		case KEY_Q:
			control_sequence();
			break;
		default:
			break;
	}

	return;
}/*}}}*/

void translateView ( double updown, double leftright){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    TRANSLATEVIEW Translates the drawing left, right, up, or down
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    Doug Jacobsen
	//

	glTranslated ( leftright, updown, 0.0);
}/*}}}*/
void polarView( double distance, double twist, double elevation, double azimuth ) {/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    POLARVIEW Rotates the drawing (in two directions), and translates in and out
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    Geoff Womeldorff
	//

	glTranslated( 0.0, 0.0, -1.0 * distance );
	glRotated( -1.0 * twist, 0.0, 0.0, 1.0 );
	glRotated( -1.0 * elevation, 1.0, 0.0, 0.0 );
	glRotated( azimuth, 0.0, 0.0, 1.0 );
}/*}}}*/

void hsv_to_rgb(float h, float s, float v, float& r, float& g, float& b){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    HSV_TO_RGB convertes a hue-saturation-value triplet to a red-green-blue triplet
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    30 December 2010
	//
	//  Author:
	//
	//    Geoff Womeldorff, Doug Jacobsen
	//

	int i;
	float f;
	float m;
	float n;

	h *= 6.0;

	i = floor( h );
	f = h - i;

	if(!(i & 1)) f = 1 - f; // if i is even
	m = v * (1 - s);
	n = v * (1 - s * f);

	switch (i) {
		case 6:
		case 0: 
			r = v; g = n; b = m;
			break;
		case 1:
			r = n; g = v; b = m;
			break;
		case 2:
			r = m; g = v; b = n;
			break;
		case 3:
			r = m; g = n; b = v;
			break;
		case 4: 
			r = n; g = m; b = v;
			break;
		case 5:
			r = v; g = m; b = n;
			break;	
	}

	return;
}/*}}}*/

double getLat(double x, double y, double z){/*{{{*/
	return asin(z); 
}/*}}}*/
double getLon(double x, double y, double z){/*{{{*/
			double lon;

			lon = atan2(y,x);

			if(lon < 0){
				return 2.0 * M_PI + lon;
			} else {
				return lon;
			}
}/*}}}*/
double getX(double lat, double lon){/*{{{*/
	return cos(lon) * cos(lat);
}/*}}}*/
double getY(double lat, double lon){/*{{{*/
	return sin(lon) * cos(lat);
}/*}}}*/
double getZ(double lat, double lon){/*{{{*/
	return sin(lat);
}/*}}}*/

void drawSphere(double r, int lats, int longs) {/*{{{*/
	int i, j;
	double z_c;

	z_c = -0.01;
	//z_c = 0.0;
	for(i = 0; i < lats; i++) {
		double lat0 = M_PI * (-0.5 + (double) i  / lats);
		double z0  = sin(lat0);
		double zr0 =  cos(lat0);

		double lat1 = M_PI * (-0.5 + (double) (i+1) / lats);
		double z1 = sin(lat1);
		double zr1 = cos(lat1);

		glBegin(GL_QUAD_STRIP);
		for(j = 0; j <= longs; j++) {
			double lng = 2 * M_PI * (double) (j) / longs;
			double x = cos(lng);
			double y = sin(lng);

			glNormal3f(x * zr0, y * zr0, z0);
			glVertex3f(r * x * zr0, r * y * zr0, r * (z0+z_c));
			glColor3f(205.0/255.0, 190.0/255.0, 112.0/255.0);
			glNormal3f(x * zr1, y * zr1, z1);
			glVertex3f(r * x * zr1, r * y * zr1, r * (z1+z_c));
			glColor3f(205.0/255.0, 190.0/255.0, 112.0/255.0);
			//glColor3f(205.0/255.0, 190.0/255.0, 112.0/255.0); // Light golden rod3
			//glColor3f(154.0/255.0, 205.0/255.0, 50.0/255.0); // Yellow Green
		}
		glEnd();
	}
}/*}}}*/

void screenshot (char filename[160],int x, int y) /*{{{*/
{// get the image data
	long imageSize = x * y * 3;
	unsigned char *data = new unsigned char[imageSize];
	glReadPixels(0,0,x,y, GL_BGR,GL_UNSIGNED_BYTE,data);// split x and y sizes into bytes
	int xa= x % 256;
	int xb= (x-xa)/256;
	int ya= y % 256;
	int yb= (y-ya)/256;

	//assemble the header
	unsigned char header[18]={0,0,2,0,0,0,0,0,0,0,0,0,(char)xa,(char)xb,(char)ya,(char)yb,24,0};

	// write header and data to file
	fstream File(filename, ios::out | ios::binary);
	File.write (reinterpret_cast<char *>(header), sizeof (char)*18);
	File.write (reinterpret_cast<char *>(data), sizeof (char)*imageSize);
	File.close();

	delete[] data;
	data=NULL;

	cout << "Screenshot written to: " << filename << endl;
}/*}}}*/
void control_sequence(){/*{{{*/
	ifstream ctrl_seq("ControlFile");
	string seq;
	string junk;
	
	cout << "Trying to read file: ControlFile" << endl;

	if(!ctrl_seq){
		cout << "  -- FAILED -- Error reading file" << endl;
		return;
	}

	getline(ctrl_seq, seq);
	const char *p = seq.c_str();

	while(*p != '\0'){
		keyPressed(*p++, 0, 0);
		display();
	}

	ctrl_seq.close();
}/*}}}*/
