#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <unistd.h>

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
int color_map = 1;
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
int 	nvertdimension;  // this is the size of the dimension being currently treated as the vertical dimension
string vertdimensionname; // this is the name of the dimension being currently treated as the vertical dimension
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

		distance1 = fabs(xcell[c1] - xvertex[v1]);
		distance2 = fabs(ycell[c1] - yvertex[v1]);
		distance3 = fabs(zcell[c1] - zvertex[v1]);

		if(c1 != c2){
			distance1 = max(distance1, fabs(xcell[c2] - xvertex[v1]));
			distance2 = max(distance2, fabs(ycell[c2] - yvertex[v1]));
			distance3 = max(distance3, fabs(zcell[c2] - zvertex[v1]));

			distance1 = max(distance1, fabs(xcell[c2] - xvertex[v2]));
			distance2 = max(distance2, fabs(ycell[c2] - yvertex[v2]));
			distance3 = max(distance3, fabs(zcell[c2] - zvertex[v2]));

			distance1 = max(distance1, fabs(xcell[c1] - xvertex[v2]));
			distance2 = max(distance2, fabs(ycell[c1] - yvertex[v2]));
			distance3 = max(distance3, fabs(zcell[c1] - zvertex[v2]));
		} else {
			distance1 = max(distance1, fabs(xvertex[v1] - xvertex[v2]));
			distance2 = max(distance2, fabs(yvertex[v1] - yvertex[v2]));
			distance3 = max(distance3, fabs(zvertex[v1] - zvertex[v2]));

			distance1 = max(distance1, fabs(xcell[c1] - xvertex[v2]));
			distance2 = max(distance2, fabs(ycell[c1] - yvertex[v2]));
			distance3 = max(distance3, fabs(zcell[c1] - zvertex[v2]));
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


void color_mapping_standard(float h, float& r, float& g, float& b){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    COLOR_MAP_STANDARD converts an array value to a r,g,b colour with a fixed 
	//    saturation and brightness and varying hue
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    6 Oct 2014
	//
	//  Author:
	//
	//    Adrian Turner
	//
	float s, v;

	s = 1.0;
	v = 1.0;

	hsv_to_rgb(h, s, v, b, g, r);

}/*}}}*/

void color_mapping_viridis(float h, float& r, float& g, float& b){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    COLOR_MAP_VIRIDIS converts an array value to a r,g,b colour using  
	//    the matplotlib viridis colourmap
	//    https://github.com/matplotlib/matplotlib/blob/196f3446a3d5178c58144cee796fa8e8aa8d2917/lib/matplotlib/_cm_listed.py
	//
	//  Licensing:
	//
	//    This code is distributed under the matplotlib license: MATPLOTLIB_LICENSE
	//
	//  Modified:
	//
	//    April 11, 2016
	//
	//  Author:
	//
	//    Phillip J. Wolfram
	//

	int value;
	/* colormap def {{{*/
	float viridis_data[256][3] = {{0.267004, 0.004874, 0.329415},
																		{0.268510, 0.009605, 0.335427},
																		{0.269944, 0.014625, 0.341379},
																		{0.271305, 0.019942, 0.347269},
																		{0.272594, 0.025563, 0.353093},
																		{0.273809, 0.031497, 0.358853},
																		{0.274952, 0.037752, 0.364543},
																		{0.276022, 0.044167, 0.370164},
																		{0.277018, 0.050344, 0.375715},
																		{0.277941, 0.056324, 0.381191},
																		{0.278791, 0.062145, 0.386592},
																		{0.279566, 0.067836, 0.391917},
																		{0.280267, 0.073417, 0.397163},
																		{0.280894, 0.078907, 0.402329},
																		{0.281446, 0.084320, 0.407414},
																		{0.281924, 0.089666, 0.412415},
																		{0.282327, 0.094955, 0.417331},
																		{0.282656, 0.100196, 0.422160},
																		{0.282910, 0.105393, 0.426902},
																		{0.283091, 0.110553, 0.431554},
																		{0.283197, 0.115680, 0.436115},
																		{0.283229, 0.120777, 0.440584},
																		{0.283187, 0.125848, 0.444960},
																		{0.283072, 0.130895, 0.449241},
																		{0.282884, 0.135920, 0.453427},
																		{0.282623, 0.140926, 0.457517},
																		{0.282290, 0.145912, 0.461510},
																		{0.281887, 0.150881, 0.465405},
																		{0.281412, 0.155834, 0.469201},
																		{0.280868, 0.160771, 0.472899},
																		{0.280255, 0.165693, 0.476498},
																		{0.279574, 0.170599, 0.479997},
																		{0.278826, 0.175490, 0.483397},
																		{0.278012, 0.180367, 0.486697},
																		{0.277134, 0.185228, 0.489898},
																		{0.276194, 0.190074, 0.493001},
																		{0.275191, 0.194905, 0.496005},
																		{0.274128, 0.199721, 0.498911},
																		{0.273006, 0.204520, 0.501721},
																		{0.271828, 0.209303, 0.504434},
																		{0.270595, 0.214069, 0.507052},
																		{0.269308, 0.218818, 0.509577},
																		{0.267968, 0.223549, 0.512008},
																		{0.266580, 0.228262, 0.514349},
																		{0.265145, 0.232956, 0.516599},
																		{0.263663, 0.237631, 0.518762},
																		{0.262138, 0.242286, 0.520837},
																		{0.260571, 0.246922, 0.522828},
																		{0.258965, 0.251537, 0.524736},
																		{0.257322, 0.256130, 0.526563},
																		{0.255645, 0.260703, 0.528312},
																		{0.253935, 0.265254, 0.529983},
																		{0.252194, 0.269783, 0.531579},
																		{0.250425, 0.274290, 0.533103},
																		{0.248629, 0.278775, 0.534556},
																		{0.246811, 0.283237, 0.535941},
																		{0.244972, 0.287675, 0.537260},
																		{0.243113, 0.292092, 0.538516},
																		{0.241237, 0.296485, 0.539709},
																		{0.239346, 0.300855, 0.540844},
																		{0.237441, 0.305202, 0.541921},
																		{0.235526, 0.309527, 0.542944},
																		{0.233603, 0.313828, 0.543914},
																		{0.231674, 0.318106, 0.544834},
																		{0.229739, 0.322361, 0.545706},
																		{0.227802, 0.326594, 0.546532},
																		{0.225863, 0.330805, 0.547314},
																		{0.223925, 0.334994, 0.548053},
																		{0.221989, 0.339161, 0.548752},
																		{0.220057, 0.343307, 0.549413},
																		{0.218130, 0.347432, 0.550038},
																		{0.216210, 0.351535, 0.550627},
																		{0.214298, 0.355619, 0.551184},
																		{0.212395, 0.359683, 0.551710},
																		{0.210503, 0.363727, 0.552206},
																		{0.208623, 0.367752, 0.552675},
																		{0.206756, 0.371758, 0.553117},
																		{0.204903, 0.375746, 0.553533},
																		{0.203063, 0.379716, 0.553925},
																		{0.201239, 0.383670, 0.554294},
																		{0.199430, 0.387607, 0.554642},
																		{0.197636, 0.391528, 0.554969},
																		{0.195860, 0.395433, 0.555276},
																		{0.194100, 0.399323, 0.555565},
																		{0.192357, 0.403199, 0.555836},
																		{0.190631, 0.407061, 0.556089},
																		{0.188923, 0.410910, 0.556326},
																		{0.187231, 0.414746, 0.556547},
																		{0.185556, 0.418570, 0.556753},
																		{0.183898, 0.422383, 0.556944},
																		{0.182256, 0.426184, 0.557120},
																		{0.180629, 0.429975, 0.557282},
																		{0.179019, 0.433756, 0.557430},
																		{0.177423, 0.437527, 0.557565},
																		{0.175841, 0.441290, 0.557685},
																		{0.174274, 0.445044, 0.557792},
																		{0.172719, 0.448791, 0.557885},
																		{0.171176, 0.452530, 0.557965},
																		{0.169646, 0.456262, 0.558030},
																		{0.168126, 0.459988, 0.558082},
																		{0.166617, 0.463708, 0.558119},
																		{0.165117, 0.467423, 0.558141},
																		{0.163625, 0.471133, 0.558148},
																		{0.162142, 0.474838, 0.558140},
																		{0.160665, 0.478540, 0.558115},
																		{0.159194, 0.482237, 0.558073},
																		{0.157729, 0.485932, 0.558013},
																		{0.156270, 0.489624, 0.557936},
																		{0.154815, 0.493313, 0.557840},
																		{0.153364, 0.497000, 0.557724},
																		{0.151918, 0.500685, 0.557587},
																		{0.150476, 0.504369, 0.557430},
																		{0.149039, 0.508051, 0.557250},
																		{0.147607, 0.511733, 0.557049},
																		{0.146180, 0.515413, 0.556823},
																		{0.144759, 0.519093, 0.556572},
																		{0.143343, 0.522773, 0.556295},
																		{0.141935, 0.526453, 0.555991},
																		{0.140536, 0.530132, 0.555659},
																		{0.139147, 0.533812, 0.555298},
																		{0.137770, 0.537492, 0.554906},
																		{0.136408, 0.541173, 0.554483},
																		{0.135066, 0.544853, 0.554029},
																		{0.133743, 0.548535, 0.553541},
																		{0.132444, 0.552216, 0.553018},
																		{0.131172, 0.555899, 0.552459},
																		{0.129933, 0.559582, 0.551864},
																		{0.128729, 0.563265, 0.551229},
																		{0.127568, 0.566949, 0.550556},
																		{0.126453, 0.570633, 0.549841},
																		{0.125394, 0.574318, 0.549086},
																		{0.124395, 0.578002, 0.548287},
																		{0.123463, 0.581687, 0.547445},
																		{0.122606, 0.585371, 0.546557},
																		{0.121831, 0.589055, 0.545623},
																		{0.121148, 0.592739, 0.544641},
																		{0.120565, 0.596422, 0.543611},
																		{0.120092, 0.600104, 0.542530},
																		{0.119738, 0.603785, 0.541400},
																		{0.119512, 0.607464, 0.540218},
																		{0.119423, 0.611141, 0.538982},
																		{0.119483, 0.614817, 0.537692},
																		{0.119699, 0.618490, 0.536347},
																		{0.120081, 0.622161, 0.534946},
																		{0.120638, 0.625828, 0.533488},
																		{0.121380, 0.629492, 0.531973},
																		{0.122312, 0.633153, 0.530398},
																		{0.123444, 0.636809, 0.528763},
																		{0.124780, 0.640461, 0.527068},
																		{0.126326, 0.644107, 0.525311},
																		{0.128087, 0.647749, 0.523491},
																		{0.130067, 0.651384, 0.521608},
																		{0.132268, 0.655014, 0.519661},
																		{0.134692, 0.658636, 0.517649},
																		{0.137339, 0.662252, 0.515571},
																		{0.140210, 0.665859, 0.513427},
																		{0.143303, 0.669459, 0.511215},
																		{0.146616, 0.673050, 0.508936},
																		{0.150148, 0.676631, 0.506589},
																		{0.153894, 0.680203, 0.504172},
																		{0.157851, 0.683765, 0.501686},
																		{0.162016, 0.687316, 0.499129},
																		{0.166383, 0.690856, 0.496502},
																		{0.170948, 0.694384, 0.493803},
																		{0.175707, 0.697900, 0.491033},
																		{0.180653, 0.701402, 0.488189},
																		{0.185783, 0.704891, 0.485273},
																		{0.191090, 0.708366, 0.482284},
																		{0.196571, 0.711827, 0.479221},
																		{0.202219, 0.715272, 0.476084},
																		{0.208030, 0.718701, 0.472873},
																		{0.214000, 0.722114, 0.469588},
																		{0.220124, 0.725509, 0.466226},
																		{0.226397, 0.728888, 0.462789},
																		{0.232815, 0.732247, 0.459277},
																		{0.239374, 0.735588, 0.455688},
																		{0.246070, 0.738910, 0.452024},
																		{0.252899, 0.742211, 0.448284},
																		{0.259857, 0.745492, 0.444467},
																		{0.266941, 0.748751, 0.440573},
																		{0.274149, 0.751988, 0.436601},
																		{0.281477, 0.755203, 0.432552},
																		{0.288921, 0.758394, 0.428426},
																		{0.296479, 0.761561, 0.424223},
																		{0.304148, 0.764704, 0.419943},
																		{0.311925, 0.767822, 0.415586},
																		{0.319809, 0.770914, 0.411152},
																		{0.327796, 0.773980, 0.406640},
																		{0.335885, 0.777018, 0.402049},
																		{0.344074, 0.780029, 0.397381},
																		{0.352360, 0.783011, 0.392636},
																		{0.360741, 0.785964, 0.387814},
																		{0.369214, 0.788888, 0.382914},
																		{0.377779, 0.791781, 0.377939},
																		{0.386433, 0.794644, 0.372886},
																		{0.395174, 0.797475, 0.367757},
																		{0.404001, 0.800275, 0.362552},
																		{0.412913, 0.803041, 0.357269},
																		{0.421908, 0.805774, 0.351910},
																		{0.430983, 0.808473, 0.346476},
																		{0.440137, 0.811138, 0.340967},
																		{0.449368, 0.813768, 0.335384},
																		{0.458674, 0.816363, 0.329727},
																		{0.468053, 0.818921, 0.323998},
																		{0.477504, 0.821444, 0.318195},
																		{0.487026, 0.823929, 0.312321},
																		{0.496615, 0.826376, 0.306377},
																		{0.506271, 0.828786, 0.300362},
																		{0.515992, 0.831158, 0.294279},
																		{0.525776, 0.833491, 0.288127},
																		{0.535621, 0.835785, 0.281908},
																		{0.545524, 0.838039, 0.275626},
																		{0.555484, 0.840254, 0.269281},
																		{0.565498, 0.842430, 0.262877},
																		{0.575563, 0.844566, 0.256415},
																		{0.585678, 0.846661, 0.249897},
																		{0.595839, 0.848717, 0.243329},
																		{0.606045, 0.850733, 0.236712},
																		{0.616293, 0.852709, 0.230052},
																		{0.626579, 0.854645, 0.223353},
																		{0.636902, 0.856542, 0.216620},
																		{0.647257, 0.858400, 0.209861},
																		{0.657642, 0.860219, 0.203082},
																		{0.668054, 0.861999, 0.196293},
																		{0.678489, 0.863742, 0.189503},
																		{0.688944, 0.865448, 0.182725},
																		{0.699415, 0.867117, 0.175971},
																		{0.709898, 0.868751, 0.169257},
																		{0.720391, 0.870350, 0.162603},
																		{0.730889, 0.871916, 0.156029},
																		{0.741388, 0.873449, 0.149561},
																		{0.751884, 0.874951, 0.143228},
																		{0.762373, 0.876424, 0.137064},
																		{0.772852, 0.877868, 0.131109},
																		{0.783315, 0.879285, 0.125405},
																		{0.793760, 0.880678, 0.120005},
																		{0.804182, 0.882046, 0.114965},
																		{0.814576, 0.883393, 0.110347},
																		{0.824940, 0.884720, 0.106217},
																		{0.835270, 0.886029, 0.102646},
																		{0.845561, 0.887322, 0.099702},
																		{0.855810, 0.888601, 0.097452},
																		{0.866013, 0.889868, 0.095953},
																		{0.876168, 0.891125, 0.095250},
																		{0.886271, 0.892374, 0.095374},
																		{0.896320, 0.893616, 0.096335},
																		{0.906311, 0.894855, 0.098125},
																		{0.916242, 0.896091, 0.100717},
																		{0.926106, 0.897330, 0.104071},
																		{0.935904, 0.898570, 0.108131},
																		{0.945636, 0.899815, 0.112838},
																		{0.955300, 0.901065, 0.118128},
																		{0.964894, 0.902323, 0.123941},
																		{0.974417, 0.903590, 0.130215},
																		{0.983868, 0.904867, 0.136897},
																		{0.993248, 0.906157, 0.143936}};
	/*}}}*/

	// interpolate	color based on h
	value = fmin(fmax((int) (h*255), 0), 255);
	r = viridis_data[value][0];
	g = viridis_data[value][1];
	b = viridis_data[value][2];

	//cout << "value = " << value << " h = " << h << " r=" << r << " g=" << g << " b=" << b << endl;

}/*}}}*/

void color_mapping_jet(float h, float& r, float& g, float& b){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    COLOR_MAP_JET converts an array value to a r,g,b colour using  
	//    the matlab jet colourmap
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    6 Oct 2014
	//
	//  Author:
	//
	//    Adrian Turner
	//

	b = fmax(fmin(fmin(4.0 * h + 0.5, -4.0 * h + 2.5),1.0),0.0);

	g = fmax(fmin(fmin(4.0 * h - 0.5, -4.0 * h + 3.5),1.0),0.0);

	r = fmax(fmin(fmin(4.0 * h - 1.5, -4.0 * h + 4.5),1.0),0.0);

}/*}}}*/

void color_mapping(double value, double min, double max, double missing_value, float& r, float& g, float& b){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    COLOR_MAP converts an array value to a r,g,b colour
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    6 Oct 2014
	//
	//  Author:
	//
	//    Adrian Turner
	//
	float h;

	if(value == missing_value){
		r = 0.5;
		g = 0.5;
		b = 0.5;
	} else {
		if(max-min != 0.0){
			if(value >= max){
				h = 1.0;
			} else if (value <= min){
				h = 0.0;
			} else {
				h = (value - min)/(max - min);
			}
		} else {
			h = (value - min)/1.0;
		}

		switch(color_map){
			case 0:
				h = h * range_factor;
				color_mapping_standard(h, r, g, b);
				break;
			case 1:
				color_mapping_jet(h, r, g, b);
				break;
			case 2:
				color_mapping_viridis(h, r, g, b);
				break;
			default:
				h = h * range_factor;
				color_mapping_standard(h, r, g, b);
				break;
		}

	}

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
        //     6 October  2014
	//
	//  Author:
	//
	//    Doug Jacobsen
        //    Adrian Turner (color map)
	//

	int i, j, k, o;
	double max, min;
	long *dims;
	int num_dims = 0;
	int cell_dim;
	int time_dim;
	int num_items;
	float r, g, b;

	cell_colors.clear();

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

			color_mapping(cell_values[i], min, max, missing_value, r, g, b);

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
        //     6 October  2014
	//
	//  Author:
	//
	//    Doug Jacobsen
        //    Adrian Turner (color map)
	//

	int i, j;
	double max, min;
	float r, g, b;

	vertex_colors.clear();

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

			color_mapping(triangle_values[i], min, max, missing_value, r, g, b);

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
        //     6 October  2014
	//
	//  Author:
	//
	//    Doug Jacobsen
        //    Adrian Turner (color map)
	//
	int i, j, o;
	double max, min;
	float r, g, b;

	edge_colors.clear();

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

			color_mapping(edge_values[i], min, max, missing_value, r, g, b);

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

	char xtime[64];
	int xtime_found;

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
			cur_level = (cur_level + 1) % nvertdimension;
			cout << "l/L keys control dimension: " << vertdimensionname << "; Current dimension value: " << cur_level+1 << " out of " << nvertdimension << endl;
			color_mesh();
			break;
		case KEY_L:
			cur_level = (cur_level - 1) % nvertdimension;
			if(cur_level < 0)
				cur_level = nvertdimension - 1;
			cout << "l/L keys control dimension: " << vertdimensionname << "; Current dimension value: " << cur_level+1 << " out of " << nvertdimension << endl;
			color_mesh();
			break;
		case KEY_t:
			if ( ntime > 0 ) {
				cur_time = (cur_time + 1) % ntime;
			}
			xtime_found = netcdf_mpas_get_xtime(filename, cur_time, &xtime[0]);
			if ( xtime_found ) {
				cout << "Current time level: " << cur_time+1 << " out of " << ntime << ".  xtime=" << xtime << endl;
			} else {
				cout << "Current time level: " << cur_time+1 << " out of " << ntime << "." << endl;
			}
			color_mesh();
			break;
		case KEY_T:
			if ( ntime > 0 ) {
				cur_time = (cur_time - 1) % ntime;
				if(cur_time < 0)
					cur_time = ntime-1;
			}
			xtime_found = netcdf_mpas_get_xtime(filename, cur_time, &xtime[0]);
			if ( xtime_found ) {
				cout << "Current time level: " << cur_time+1 << " out of " << ntime << ".  xtime=" << xtime << endl;
			} else {
				cout << "Current time level: " << cur_time+1 << " out of " << ntime << "." << endl;
			}
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
		case KEY_m:
			color_map = (color_map + 1) % 3;
			switch(color_map){
				case 0:
					cout << "Color mapping based on variable hue" << endl;
					break;
				case 1:
					cout << "Color mapping with matlab jet" << endl;
					break;
				case 2:
					cout << "Color mapping with matplotlib viridis" << endl;
					break;
				default:
					break;
			}
			color_mesh();
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
					netcdf_mpas_get_vert_dim_info(filename, vertex_field, nvertdimension, vertdimensionname); // Update vertical dimension info
					cur_level = cur_level % nvertdimension;
					cout << "l/L keys control dimension: " << vertdimensionname << "; Current dimension value: " << cur_level+1 << " out of " << nvertdimension << endl;
					if(vertex_field >= 0 && color_bar == 2){
						build_range(vertex_field);
						cout << " Range of values: Min = " << ranges.at(vertex_field).at(0) << ", Max = " << ranges.at(vertex_field).at(1) << endl;
					}
					break;
				case 1:
					cell_field = netcdf_mpas_list_ncell_fields(filename);
					netcdf_mpas_print_field_info(filename, cell_field);
					netcdf_mpas_get_vert_dim_info(filename, cell_field, nvertdimension, vertdimensionname); // Update vertical dimension info
					cur_level = cur_level % nvertdimension;
					cout << "l/L keys control dimension: " << vertdimensionname << "; Current dimension value: " << cur_level+1 << " out of " << nvertdimension << endl;
					if(cell_field >= 0 && color_bar == 2){
						build_range(cell_field);
						cout << " Range of values: Min = " << ranges.at(cell_field).at(0) << ", Max = " << ranges.at(cell_field).at(1) << endl;
					}
					break;
				case 2:
					edge_field = netcdf_mpas_list_nedge_fields(filename);
					netcdf_mpas_print_field_info(filename, edge_field);
					netcdf_mpas_get_vert_dim_info(filename, edge_field, nvertdimension, vertdimensionname); // Update vertical dimension info
					cur_level = cur_level % nvertdimension;
					cout << "l/L keys control dimension: " << vertdimensionname << "; Current dimension value: " << cur_level+1 << " out of " << nvertdimension << endl;
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
