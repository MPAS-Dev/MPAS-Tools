#include <string>
#include <string.h>
#include <netcdfcpp.h>
#include <cstdlib>
#include <vector>

using namespace std;

int netcdf_mpas_read_num_vars(string filename){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_NUM_VARS gets the number of variables in the netcdf file
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    11 February 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Output, int NETCDF_MPAS_READ_NUM_VARS, the number of variables
	//
	long int var_size;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get Ntime, which is a NETCDF dimension.
	//
	var_size = ncid.num_vars();

	ncid.close();

	return var_size;
}/*}}}*/

/* Attribute reading functions {{{*/
bool netcdf_mpas_read_onsphere(string filename){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_ONSPHERE determines if points are on a sphere, or in a plane
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    11 February 2014
	//
	//  Author:
	//
	//     Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Output, int NETCDF_MPAS_READ_ONSPHERE, true or false based on the global
	//    										 attribute, on_a_sphere
	//
	NcAtt *att_id;
	NcValues *vals;
	bool valid;
	string tmp_name;
	string sph_name = "on_a_sphere";
	//
	//  Open the file.
	//
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//  Get Ntime, which is a NETCDF dimension.
	//
	valid = false;

	for(int i = 0; i < ncid.num_atts() && !valid; i++){
		tmp_name = ncid.get_att(i)->name();

		if(sph_name == tmp_name){
			valid = true;
		}
	}

	if(valid){
		att_id = ncid.get_att(sph_name.c_str());
		
		vals = att_id -> values();

		sph_name = "YES";

		for(int i = 0; i < vals -> num(); i++){
			tmp_name = vals -> as_string(i);
			if(tmp_name.find(sph_name) != string::npos){
				return true;
			}
		}
	} else {
		return true;
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return false;

}/*}}}*/
//****************************************************************************80
double netcdf_mpas_read_sphereradius(string filename){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_SPHERERADIUS determines if the radius of the sphere.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    11 February 2014
	//
	//  Author:
	//
	//     Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Output, double NETCDF_MPAS_READ_SPHERERADIUS, true or false based on the global
	//    										 attribute, on_a_sphere
	//
	NcAtt *att_id;
	NcValues *vals;
	bool valid;
	string tmp_name;
	string sph_name = "sphere_radius";
	//
	//  Open the file.
	//
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//  Get Ntime, which is a NETCDF dimension.
	//
	valid = false;

	for(int i = 0; i < ncid.num_atts() && !valid; i++){
		tmp_name = ncid.get_att(i)->name();

		if(sph_name == tmp_name){
			valid = true;
		}
	}

	if(valid){
		att_id = ncid.get_att(sph_name.c_str());
		
		vals = att_id -> values();

		sph_name = "YES             ";

		for(int i = 0; i < vals -> num(); i++){
			tmp_name = vals -> as_string(i);

			return atof(vals -> as_string(i));
		}
	} else {
		return 0.0;
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return 0.0;

}/*}}}*/
//****************************************************************************80
bool netcdf_mpas_read_isperiodic(string filename){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_ISPERIODIC determines if a mesh is intended to be periodic.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    7 March 2014
	//
	//  Author:
	//
	//     Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Output, int NETCDF_MPAS_READ_ISPERIODIC, true or false based on the global
	//    										 attribute, is_periodic. Default is false
	//
	NcAtt *att_id;
	NcValues *vals;
	bool valid;
	string tmp_name;
	string sph_name = "is_periodic";
	//
	//  Open the file.
	//
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//  Get Ntime, which is a NETCDF dimension.
	//
	valid = false;

	for(int i = 0; i < ncid.num_atts() && !valid; i++){
		tmp_name = ncid.get_att(i)->name();

		if(sph_name == tmp_name){
			valid = true;
		}
	}

	if(valid){
		att_id = ncid.get_att(sph_name.c_str());
		
		vals = att_id -> values();

		sph_name = "YES";

		for(int i = 0; i < vals -> num(); i++){
			tmp_name = vals -> as_string(i);
			if(tmp_name.find(sph_name) != string::npos){
				return true;
			}
		}
	} else {
		return false;
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return false;

}/*}}}*/
//****************************************************************************80
double netcdf_mpas_read_xperiod(string filename){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_XPERIOD reads the global attribute x_period
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    11 February 2014
	//
	//  Author:
	//
	//     Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Output, double NETCDF_MPAS_READ_XPERIOD, the value of the global
	//    										 attribute, x_period
	//
	NcAtt *att_id;
	NcValues *vals;
	bool valid;
	string tmp_name;
	string sph_name = "x_period";
	//
	//  Open the file.
	//
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//  Get Ntime, which is a NETCDF dimension.
	//
	valid = false;

	for(int i = 0; i < ncid.num_atts() && !valid; i++){
		tmp_name = ncid.get_att(i)->name();

		if(sph_name == tmp_name){
			valid = true;
		}
	}

	if(valid){
		att_id = ncid.get_att(sph_name.c_str());
		
		vals = att_id -> values();

		for(int i = 0; i < vals -> num(); i++){
			tmp_name = vals -> as_string(i);

			return atof(vals -> as_string(i));
		}
	} else {
		return -1.0;
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return -1.0;

}/*}}}*/
//****************************************************************************80
double netcdf_mpas_read_yperiod(string filename){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_YPERIOD reads the global attribute y_period
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    11 February 2014
	//
	//  Author:
	//
	//     Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Output, double NETCDF_MPAS_READ_YPERIOD, the value of the global
	//    										 attribute, y_period
	//
	NcAtt *att_id;
	NcValues *vals;
	bool valid;
	string tmp_name;
	string sph_name = "y_period";
	//
	//  Open the file.
	//
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//  Get Ntime, which is a NETCDF dimension.
	//
	valid = false;

	for(int i = 0; i < ncid.num_atts() && !valid; i++){
		tmp_name = ncid.get_att(i)->name();

		if(sph_name == tmp_name){
			valid = true;
		}
	}

	if(valid){
		att_id = ncid.get_att(sph_name.c_str());
		
		vals = att_id -> values();

		for(int i = 0; i < vals -> num(); i++){
			tmp_name = vals -> as_string(i);

			return atof(vals -> as_string(i));
		}
	} else {
		return -1.0;
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return -1.0;

}/*}}}*/
//****************************************************************************80
string netcdf_mpas_read_history(string filename){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_HISTORY reads the history attribute from a file.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    11 February 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Output, string NETCDF_MPAS_READ_HISTORY, the value of the history attribute
	//
	NcAtt *att_id;
	NcValues *vals;
	bool valid;
	string tmp_name;
	string sph_name = "history";
	string history = "";
	//
	//  Open the file.
	//
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//  Get Ntime, which is a NETCDF dimension.
	//
	valid = false;

	for(int i = 0; i < ncid.num_atts() && !valid; i++){
		tmp_name = ncid.get_att(i)->name();

		if(sph_name == tmp_name){
			valid = true;
		}
	}

	if(valid){
		att_id = ncid.get_att(sph_name.c_str());
		
		vals = att_id -> values();

		for(int i = 0; i < vals -> num(); i++){
			history += vals -> as_string(i);
			return history;
		}
	} else {
		return history;
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return history;

}/*}}}*/
//****************************************************************************80
string netcdf_mpas_read_fileid(string filename){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_FILEID reads the file_id attribute from a file.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    11 February 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Output, string NETCDF_MPAS_READ_FILEID, the value of the file_id attribute
	//
	NcAtt *att_id;
	NcValues *vals;
	bool valid;
	string tmp_name;
	string sph_name = "file_id";
	string id_str = "";
	//
	//  Open the file.
	//
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//  Get Ntime, which is a NETCDF dimension.
	//
	valid = false;

	for(int i = 0; i < ncid.num_atts() && !valid; i++){
		tmp_name = ncid.get_att(i)->name();

		if(sph_name == tmp_name){
			valid = true;
		}
	}

	if(valid){
		att_id = ncid.get_att(sph_name.c_str());
		
		vals = att_id -> values();

		for(int i = 0; i < vals -> num(); i++){
			id_str += vals -> as_string(i);
			return id_str;
		}
	} else {
		return id_str;
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return id_str;

}/*}}}*/
//****************************************************************************80
string netcdf_mpas_read_parentid(string filename){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_PARENTID reads the parent_id attribute from a file.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    11 February 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Output, string NETCDF_MPAS_READ_PARENTID, the value of the parent_id attribute
	//
	NcAtt *att_id;
	NcValues *vals;
	bool valid;
	string tmp_name;
	string sph_name = "parent_id";
	string id_str = "";
	//
	//  Open the file.
	//
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif

//	NcError err(NcError::silent_nonfatal); // Don't error if the variable isn't found.
	//
	//  Get Ntime, which is a NETCDF dimension.
	//
	valid = false;

	for(int i = 0; i < ncid.num_atts() && !valid; i++){
		tmp_name = ncid.get_att(i)->name();

		if(sph_name == tmp_name){
			valid = true;
		}
	}

	if(valid){
		att_id = ncid.get_att(sph_name.c_str());
		
		vals = att_id -> values();

		for(int i = 0; i < vals -> num(); i++){
			id_str += vals -> as_string(i);
			return id_str;
		}
	} else {
		return id_str;
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return id_str;

}/*}}}*/
//****************************************************************************80
double netcdf_mpas_read_meshspec(string filename){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_MESHSPEC returns the mesh_spec attribute.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    18 February 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Output, double NETCDF_MPAS_READ_MESHSPEC, the value of the mesh_spec attribute
	//
	NcAtt *att_id;
	NcValues *vals;
	bool valid;
	string tmp_name;
	string sph_name = "mesh_spec";
	//
	//  Open the file.
	//
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//  Get Ntime, which is a NETCDF dimension.
	//
	valid = false;

	for(int i = 0; i < ncid.num_atts() && !valid; i++){
		tmp_name = ncid.get_att(i)->name();

		if(sph_name == tmp_name){
			valid = true;
		}
	}

	if(valid){
		att_id = ncid.get_att(sph_name.c_str());
		
		vals = att_id -> values();

		sph_name = "YES             ";

		for(int i = 0; i < vals -> num(); i++){
			tmp_name = vals -> as_string(i);

			return atof(vals -> as_string(i));
		}
	} else {
		return 0.0;
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return 0.0;

}/*}}}*/
/*}}}*/

/* Dimension reading functions {{{*/
//****************************************************************************80
int netcdf_mpas_read_dim ( string filename, string dim_name ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_DIM gets the size of the dimension with name dim_name
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    11 February 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, string DIM_NAME, the name of the dimension in question
	//
	//    Output, int NETCDF_MPAS_READ_DIM, the value of the dimension.
	//
	int ntime;
	long int dim_size;
	NcDim *dim_id;
	bool valid;
	string tmp_name;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get Ntime, which is a NETCDF dimension.
	//
	valid = false;

	for(int i = 0; i < ncid.num_dims() && !valid; i++){
		tmp_name = ncid.get_dim(i)->name();

		if(dim_name == tmp_name){
			valid = true;
		}
	}

	if(valid){
		dim_id = ncid.get_dim(dim_name.c_str());
		dim_size = (*dim_id).size ( );
	} else {
		dim_size = 1;
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	ntime = ( int ) dim_size;

	return ntime;
}/*}}}*/
/*}}}*/

/* Cell Reading Functions {{{*/
//****************************************************************************80
void netcdf_mpas_read_xyzcell ( string filename, int ncells, double xcell[], double ycell[], double zcell[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_XYZCELL reads xCell, yCell, zCell.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    11 February 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NCELLS, the number of nodes.
	//
	//    Output, double XCELL[NCELLS], YCELL[NCELLS], ZCELL[NCELLS], the
	//    coordinates of the nodes.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "    Reading xCell" << endl;
#endif
	var_id = ncid.get_var ( "xCell" );
	(*var_id).get ( &xcell[0], ncells );
#ifdef _DEBUG
	cout << "    Reading yCell" << endl;
#endif
	var_id = ncid.get_var ( "yCell" );
	(*var_id).get ( &ycell[0], ncells );
#ifdef _DEBUG
	cout << "    Reading zCell" << endl;
#endif
	var_id = ncid.get_var ( "zCell" );
	(*var_id).get ( &zcell[0], ncells );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_latloncell ( string filename, int ncells, double latcell[], double loncell[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_LATLONCELL reads latCell, lonCell
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    11 February 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NCELLS, the number of nodes.
	//
	//    Output, double LATCELL[NCELLS], LONCELL[NCELLS], the coordinates of the cell centers.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading latCell" << endl;
#endif
	var_id = ncid.get_var ( "latCell" );
	(*var_id).get ( &latcell[0], ncells );
#ifdef _DEBUG
	cout << "   Reading lonCell" << endl;
#endif
	var_id = ncid.get_var ( "lonCell" );
	(*var_id).get ( &loncell[0], ncells );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_areacell ( string filename, int ncells, double areacell[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_AREACELL reads areaCell
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    11 February 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NCELLS, the number of nodes.
	//
	//    Output, double AREACELL[NCELLS] the area of each cell
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading areaCell" << endl;
#endif
	var_id = ncid.get_var ( "areaCell" );
	(*var_id).get ( &areacell[0], ncells );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_nedgesoncell ( string filename, int ncells, int nedgesoncell[] ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_NEDGESONCELL gets the nEdgesOnCell information.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    08/08/2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NCELLS, the number of cells.
	//
	//    Output, int NEDGESONCELL[NCELLS];
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading nEdgesOnCell" << endl;
#endif
	var_id = ncid.get_var ( "nEdgesOnCell" );
	(*var_id).get ( &nedgesoncell[0], ncells );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_cellsoncell ( string filename, int ncells, int maxedges, int cellsoncell[] ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_CELLSONCELL gets the cellsOnCell information.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    02/12/2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NCELLS, the number of cells.
	//
	//    Input, int MAXEDGES, the maximum number of edges around each cell.
	//
	//    Output, int CELLSONCELL[MAXEDGES*NCELLS];
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading cellsOnCell" << endl;
#endif
	var_id = ncid.get_var ( "cellsOnCell" );
	(*var_id).get ( &cellsoncell[0], ncells, maxedges );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_edgesoncell ( string filename, int ncells, int maxedges, int edgesoncell[] ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_EDGESONCELL gets the edgesOnCell information.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    02/12/2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NCELLS, the number of cells.
	//
	//    Input, int MAXEDGES, the maximum number of edges around each cell.
	//
	//    Output, int EDGESONCELL[MAXEDGES*NCELLS];
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading edgesOnCell" << endl;
#endif
	var_id = ncid.get_var ( "edgesOnCell" );
	(*var_id).get ( &edgesoncell[0], ncells, maxedges );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_verticesoncell ( string filename, int ncells, int maxedges, int verticesoncell[] ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_VERTICESONCELL gets the verticesOnCell information.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    02/12/2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NCELLS, the number of cells.
	//
	//    Input, int MAXEDGES, the maximum number of edges around each cell.
	//
	//    Output, int VERTICESONCELL[MAXEDGES*NCELLS];
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading verticesOnCell" << endl;
#endif
	var_id = ncid.get_var ( "verticesOnCell" );
	(*var_id).get ( &verticesoncell[0], ncells, maxedges );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_mesh_density ( string filename, int ncells, double mesh_density[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_MESH_DENSITY reads meshDensity
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    11 February 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NCELLS, the number of nodes.
	//
	//    Output, double mesh_density[NCELLS] value of the density function at each cell center.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
	NcError err(NcError::silent_nonfatal); // Don't error if the variable isn't found.
	
#ifdef _DEBUG
	cout << "   Reading meshDensity" << endl;
#endif
	var_id = ncid.get_var ( "meshDensity" );
	if(var_id == NULL){
		for(int i = 0; i < ncells; i++){
			mesh_density[i] = 1.0;
		}
	} else {
		(*var_id).get ( &mesh_density[0], ncells );
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_cullcell ( string filename, int ncells, int cullcell[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_MESH_CULLCELL reads cullCell
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    11 February 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NCELLS, the number of nodes.
	//
	//    Output, int cullcell[NCELLS] a mask on cells, 1 if the cell is kept, 0 if it's removed.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	NcError err(NcError::silent_nonfatal); // Don't error if the variable isn't found.
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading cullCell" << endl;
#endif
	var_id = ncid.get_var ( "cullCell" );
	if(var_id == NULL){
		for(int i = 0; i < ncells; i++){
			cullcell[i] = 0;
		}
	} else {
		(*var_id).get ( &cullcell[0], ncells );
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_regioncellmasks ( string filename, int ncells, int nregions, int regioncellmasks[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_REGIONCELLMASKS reads regionCellMasks
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    21 January 2016
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NCELLS, the number of nodes.
	//
	//    Input, int NREGIONS, the number of regions.
	//
	//    Output, int REGIONCELLMASKS[NCELLS] a mask on cells where each region has a mask of 1 or 0 if the cell is within the region.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	NcError err(NcError::silent_nonfatal); // Don't error if the variable isn't found.
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading regionCellMasks" << endl;
#endif
	var_id = ncid.get_var ( "regionCellMasks" );
	if(var_id == NULL){
		for(int i = 0; i < ncells; i++){
			for(int j = 0; j < nregions; j++){
				regioncellmasks[i*nregions + j] = 0;
			}
		}
	} else {
		(*var_id).get ( &regioncellmasks[0], ncells, nregions );
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_transectcellmasks ( string filename, int ncells, int ntransects, int transectcellmasks[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_TRANSECTCELLMASKS reads transectCellMasks
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    21 January 2016
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NCELLS, the number of nodes.
	//
	//    Input, int NTRANSECTS, the number of transects.
	//
	//    Output, int TRANSECTCELLMASKS[NCELLS] a mask on cells where each transect has a mask of 1 or 0 if the cell is part of the transect.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	NcError err(NcError::silent_nonfatal); // Don't error if the variable isn't found.
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading transectCellMasks" << endl;
#endif
	var_id = ncid.get_var ( "transectCellMasks" );
	if(var_id == NULL){
		for(int i = 0; i < ncells; i++){
			for(int j = 0; j < ntransects; j++){
				transectcellmasks[i*ntransects + j] = 0;
			}
		}
	} else {
		(*var_id).get ( &transectcellmasks[0], ncells, ntransects );
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_cellseedmask ( string filename, int ncells, int cellseedmask[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_CELLSEEDMASK reads cellSeedMask
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    21 January 2016
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NCELLS, the number of nodes.
	//
	//    Output, int CELLSEEDMASK[NCELLS] a mask on cells where a seeded flood fill has marked cells.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	NcError err(NcError::silent_nonfatal); // Don't error if the variable isn't found.
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading cellSeedMask" << endl;
#endif
	var_id = ncid.get_var ( "cellSeedMask" );
	if(var_id == NULL){
		for(int i = 0; i < ncells; i++){
			cellseedmask[i] = 0;
		}
	} else {
		(*var_id).get ( &cellseedmask[0], ncells );
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
/* }}} */

/* Vertex Reading Functions {{{*/
//****************************************************************************80
void netcdf_mpas_read_xyzvertex ( string filename, int nvertices, double xvertex[], double yvertex[], double zvertex[] ){/*{{{*/
	//****************************************************************************80

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_CELLS gets the cell center coordinates.
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
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NVERTICES, the number of vertices.
	//
	//    Output, double XVERTEX[NVERTICES], YVERTEXL[NVERTICES], 
	//    ZVERTEX[NVERTICES], the coordinates of the nodes.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading xVertex" << endl;
#endif
	var_id = ncid.get_var ( "xVertex" );
	(*var_id).get ( &xvertex[0], nvertices );
#ifdef _DEBUG
	cout << "   Reading yVertex" << endl;
#endif
	var_id = ncid.get_var ( "yVertex" );
	(*var_id).get ( &yvertex[0], nvertices );
#ifdef _DEBUG
	cout << "   Reading zVertex" << endl;
#endif
	var_id = ncid.get_var ( "zVertex" );
	(*var_id).get ( &zvertex[0], nvertices );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_latlonvertex ( string filename, int nvertices, double latvertex[], double lonvertex[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_LATLONVERTEX reads latVertex, lonVertex
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    11 February 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NVERTICES, the number of nodes.
	//
	//    Output, double LATVERTEX[NVERTICES], LONVERTEX[NVERTICES], the coordinates of the vertices.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading latVertex" << endl;
#endif
	var_id = ncid.get_var ( "latVertex" );
	(*var_id).get ( &latvertex[0], nvertices );
#ifdef _DEBUG
	cout << "   Reading lonVertex" << endl;
#endif
	var_id = ncid.get_var ( "lonVertex" );
	(*var_id).get ( &lonvertex[0], nvertices );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_areatriangle ( string filename, int nvertices, double areatriangle[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_AREATRIANGLE reads areaTriangle
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    02/12/2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NVERTICES, the number of nodes.
	//
	//    Output, double AREATRIANGLE[NVERTICES] the area of each dual cell
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading areaTriangle" << endl;
#endif
	var_id = ncid.get_var ( "areaTriangle" );
	(*var_id).get ( &areatriangle[0], nvertices );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_cellsonvertex ( string filename, int nvertices, int vertexdegree, int cellsonvertex[] ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_CELLSONVERTEX gets the cellsOnVertex information.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    13 February 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NEDGES, the number of edges.
	//
	//    Output, int CELLSONVERTEX[3*NVERTICES];
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading cellsOnVertex" << endl;
#endif
	var_id = ncid.get_var ( "cellsOnVertex" );
	(*var_id).get ( &cellsonvertex[0], nvertices, vertexdegree );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_kiteareasonvertex ( string filename, int nvertices, int vertexdegree, double kiteareasonvertex[] ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_KITEAREASONVERTEX gets the kiteAreasOnVertex information.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    02/13/2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NEDGES, the number of edges.
	//
	//    Input, int VERTEXDEGREE, the number of cells the border each vertex.
	//
	//    Output, double KITEAREASONVERTEX[VERTEXDEGREE*NVERTICES];
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading kiteAreasOnVertex" << endl;
#endif
	var_id = ncid.get_var ( "kiteAreasOnVertex" );
	(*var_id).get ( &kiteareasonvertex[0], nvertices, vertexdegree );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_edgesonvertex ( string filename, int nvertices, int vertexdegree, int edgesonvertex[] ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_EDGESONVERTEX gets the edgesOnVertex information.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    02/13/2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NVERTICES, the number of vertices.
	//
	//    Input, int VERTEXDEGREE, the number of cells that border each vertex.
	//
	//    Output, int EDGESONVERTEX[VERTEXDEGREE*NVERTICES], the indices of edges that radiate from a vertex.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading edgesOnVertex" << endl;
#endif
	var_id = ncid.get_var ( "edgesOnVertex" );
	(*var_id).get ( &edgesonvertex[0], nvertices, vertexdegree );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/

/* }}} */

/* Edge Reading Functions {{{*/
//****************************************************************************80
void netcdf_mpas_read_xyzedge ( string filename, int nedges, double xedge[], double yedge[], double zedge[] ){/*{{{*/
	//****************************************************************************80

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_XYZEDGES gets the edge coordinates.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    12 Feb. 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NEDGES, the number of edges.
	//
	//    Output, double XEDGE[NEDGES], YEDGE[NEDGES], 
	//    ZEDGE[NEDGES], the coordinates of the edges.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading xEdge" << endl;
#endif
	var_id = ncid.get_var ( "xEdge" );
	(*var_id).get ( &xedge[0], nedges );
#ifdef _DEBUG
	cout << "   Reading yEdge" << endl;
#endif
	var_id = ncid.get_var ( "yEdge" );
	(*var_id).get ( &yedge[0], nedges );
#ifdef _DEBUG
	cout << "   Reading zEdge" << endl;
#endif
	var_id = ncid.get_var ( "zEdge" );
	(*var_id).get ( &zedge[0], nedges );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_latlonedge ( string filename, int nedges, double latedge[], double lonedge[]){/*{{{*/
	//****************************************************************************80

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_LATLONEDGE gets the edge coordinates.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    12 Feb. 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NEDGES, the number of edges.
	//
	//    Output, double LATEDGE[NEDGES], LONEDGE[NEDGES], the coordinates of the edges.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading latEdge" << endl;
#endif
	var_id = ncid.get_var ( "latEdge" );
	(*var_id).get ( &latedge[0], nedges );
#ifdef _DEBUG
	cout << "   Reading lonEdge" << endl;
#endif
	var_id = ncid.get_var ( "lonEdge" );
	(*var_id).get ( &lonedge[0], nedges );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_verticesonedge ( string filename, int nedges, int verticesonedge[] ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_VERTICESONEDGE gets the verticesOnEdge information.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    13 February 2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NEDGES, the number of edges.
	//
	//    Output, int VERTICESONEDGE[2*NVERTICES];
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading verticesOnEdge" << endl;
#endif
	var_id = ncid.get_var ( "verticesOnEdge" );
	(*var_id).get ( &verticesonedge[0], nedges, 2 );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_cellsonedge ( string filename, int nedges, int cellsonedge[] ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_CELLSONEDGE gets the cellsOnEdge information.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    02/13/2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NEDGES, the number of edges.
	//
	//    Output, int CELLSONEDGE[2*NVERTICES];
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading cellsOnEdge" << endl;
#endif
	var_id = ncid.get_var ( "cellsOnEdge" );
	(*var_id).get ( &cellsonedge[0], nedges, 2 );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_dvedge ( string filename, int nedges, double dvedge[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_DVEDGE reads dvEdge
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    02/12/2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NEDGES, the number of nodes.
	//
	//    Output, double DVEDGE[NEDGES] the length between vertices of an edge.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading dvEdge" << endl;
#endif
	var_id = ncid.get_var ( "dvEdge" );
	(*var_id).get ( &dvedge[0], nedges );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_dcedge ( string filename, int nedges, double dcedge[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_DVEDGE reads dcEdge
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    02/12/2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NEDGES, the number of nodes.
	//
	//    Output, double DVEDGE[NEDGES] the length between the cells on an edge.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading dcEdge" << endl;
#endif
	var_id = ncid.get_var ( "dcEdge" );
	(*var_id).get ( &dcedge[0], nedges );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_angleedge ( string filename, int nedges, double angleedge[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_ANGLEEDGE reads angleEdge
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    02/12/2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NEDGES, the number of nodes.
	//
	//    Output, double ANGLEEDGE[NEDGES] the length between the cells on an edge.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading angleEdge" << endl;
#endif
	var_id = ncid.get_var ( "angleEdge" );
	(*var_id).get ( &angleedge[0], nedges );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_weightsonedge ( string filename, int nedges, int maxedges2, double weightsonedge[] ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_WEIGHTSONEDGE gets the weightsOnEdge information.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    02/13/2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NEDGES, the number of edges.
	//
	//    Input, int MAXEDGES2, two times the maximum number of edges around a cell.
	//
	//    Output, double WEIGHTSONEDGE[MAXEDGES2*NEDGES], the weights for each edge on edge.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading weightsOnEdge" << endl;
#endif
	var_id = ncid.get_var ( "weightsOnEdge" );
	(*var_id).get ( &weightsonedge[0], nedges, maxedges2 );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_edgesonedge ( string filename, int nedges, int maxedges2, int edgesonedge[] ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_EDGESONEDGE gets the edgesOnEdge information.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    02/13/2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NEDGES, the number of edges.
	//
	//    Input, int MAXEDGES2, two times the maximum number of edges around a cell.
	//
	//    Output, int EDGESONEDGE[MAXEDGES2*NEDGES], the weights for each edge on edge.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading edgesOnEdge" << endl;
#endif
	var_id = ncid.get_var ( "edgesOnEdge" );
	(*var_id).get ( &edgesonedge[0], nedges, maxedges2 );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_nedgesonedge ( string filename, int nedges, int nedgesonedge[] ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_NEDGESONEDGE gets the nEdgesOnEdge information.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    06/09/2014
	//
	//  Author:
	//
	//    Doug Jacobsen
	//
	//  Reference:
	//
	//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
	//    The NETCDF User's Guide,
	//    Unidata Program Center, March 2009.
	//
	//  Parameters:
	//
	//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
	//
	//    Input, int NEDGES, the number of edges.
	//
	//    Output, int NEDGESONEDGE[NEDGES], the number of edges on each edge.
	//
	NcVar *var_id;
	//
	//  Open the file.
	#ifdef _64BITPERIOD
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly, NULL, 0, NcFile::Period64Bits );
	#else
		NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	#endif
	//
	//
	//  Get the variable values.
	//
#ifdef _DEBUG
	cout << "   Reading nEdgesOnEdge" << endl;
#endif
	var_id = ncid.get_var ( "nEdgesOnEdge" );
	(*var_id).get ( &nedgesonedge[0], nedges );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/

/* }}} */

