#include <string>
#include <netcdfcpp.h>

using namespace std;

int netcdf_mpas_read_num_vars(string filename);

/* Attribute reading functions {{{*/
bool netcdf_mpas_read_onsphere(string filename);
double netcdf_mpas_read_sphereradius(string filename);
bool netcdf_mpas_read_isperiodic(string filename);
double netcdf_mpas_read_xperiod(string filename);
double netcdf_mpas_read_yperiod(string filename);
string netcdf_mpas_read_history(string filename);
string netcdf_mpas_read_fileid(string filename);
string netcdf_mpas_read_parentid(string filename);
double netcdf_mpas_read_meshspec(string filename);
/*}}}*/

/* Dimension reading functions {{{*/
int netcdf_mpas_read_dim ( string filename, string dim_name );
/*}}}*/

/* Cell Reading Functions {{{*/
void netcdf_mpas_read_xyzcell ( string filename, int ncells, double xcell[], double ycell[], double zcell[] );
void netcdf_mpas_read_latloncell ( string filename, int ncells, double latcell[], double loncell[] );
void netcdf_mpas_read_areacell ( string filename, int ncells, double areacell[] );
void netcdf_mpas_read_nedgesoncell ( string filename, int ncells, int edgesoncell[] );
void netcdf_mpas_read_cellsoncell ( string filename, int ncells, int maxedges, int cellsoncell[] );
void netcdf_mpas_read_edgesoncell ( string filename, int ncells, int maxedges, int edgesoncell[] );
void netcdf_mpas_read_verticesoncell ( string filename, int ncells, int maxedges, int verticesoncell[] );
void netcdf_mpas_read_mesh_density ( string filename, int ncells, double mesh_density[] );
void netcdf_mpas_read_cullcell ( string filename, int ncells, int cullcell[] );
void netcdf_mpas_read_regioncellmasks ( string filename, int ncells, int nregions, int regioncellmasks[] );
void netcdf_mpas_read_transectcellmasks ( string filename, int ncells, int ntransects, int transectcellmasks[] );
void netcdf_mpas_read_cellseedmask ( string filename, int ncells, int cellseedmask[] );
/* }}} */

/* Vertex Reading Functions {{{*/
void netcdf_mpas_read_xyzvertex ( string filename, int nvertices, double xvertex[], double yvertex[], double zvertex[] );
void netcdf_mpas_read_latlonvertex ( string filename, int nvertices, double latvertex[], double lonvertex[] );
void netcdf_mpas_read_areatriangle ( string filename, int nvertices, double areatriangle[] );
void netcdf_mpas_read_cellsonvertex ( string filename, int nvertices, int vertexdegree, int cellsonvertex[] );
void netcdf_mpas_read_kiteareasonvertex ( string filename, int nvertices, int vertexdegree, double kiteareasonvertex[] );
void netcdf_mpas_read_edgesonvertex ( string filename, int nvertices, int vertexdegree, int edgesonvertex[] );
/* }}} */

/* Edge Reading Functions {{{*/
void netcdf_mpas_read_xyzedge ( string filename, int nedges, double xedge[], double yedge[], double zedge[] );
void netcdf_mpas_read_latlonedge ( string filename, int nedges, double latedge[], double lonedge[]);
void netcdf_mpas_read_verticesonedge ( string filename, int nedges, int verticesonedge[] );
void netcdf_mpas_read_cellsonedge ( string filename, int nedges, int cellsonedge[] );
void netcdf_mpas_read_dvedge ( string filename, int nedges, double dvedge[] );
void netcdf_mpas_read_dcedge ( string filename, int nedges, double dcedge[] );
void netcdf_mpas_read_angleedge ( string filename, int nedges, double angleedge[] );
void netcdf_mpas_read_weightsonedge ( string filename, int nedges, int maxedges2, double weightsonedge[] );
void netcdf_mpas_read_edgesonedge ( string filename, int nedges, int maxedges2, int edgesonedge[] );
void netcdf_mpas_read_nedgesonedge ( string filename, int nedges, int nedgesonedge[] );
/* }}} */

