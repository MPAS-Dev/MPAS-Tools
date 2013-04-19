#include <algorithm>

using namespace std;

double r8_huge ( ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    R8_HUGE returns a "huge" R8.
	//
	//  Discussion:
	//
	//    The value returned by this function is NOT required to be the
	//    maximum representable R8.  This value varies from machine to machine,
	//    from compiler to compiler, and may cause problems when being printed.
	//    We simply want a "very large" but non-infinite number.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    06 October 2007
	//
	//  Author:
	//
	//    John Burkardt
	//
	//  Parameters:
	//
	//    Output, double R8_HUGE, a "huge" R8 value.
	//
	double value;

	value = 1.0E+30;

	return value;
}/*}}}*/
//****************************************************************************80
double r8vec_max ( int n, double r8vec[] , double missing_value){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    R8VEC_MAX returns the value of the maximum element in an R8VEC.
	//
	//  Discussion:
	//
	//    An R8VEC is a vector of R8's.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    22 August 2010
	//
	//  Author:
	//
	//    John Burkardt
	//
	//  Parameters:
	//
	//    Input, int N, the number of entries in the array.
	//
	//    Input, double R8VEC[N], a pointer to the first entry of the array.
	//
	//    Output, double R8VEC_MAX, the value of the maximum element.  This
	//    is set to 0.0 if N <= 0.
	//
	int i;
	double value;

	if ( n <= 0 )
	{
		value = 0.0;
	}

	value = -1.0e34;

	for ( i = 0; i < n; i++ )
	{
		if(r8vec[i] != missing_value)
			value = max(value, r8vec[i]);
	}
	return value;
}/*}}}*/
//****************************************************************************80
double r8vec_min ( int n, double r8vec[] , double missing_value){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    R8VEC_MIN returns the value of the minimum element in an R8VEC.
	//
	//  Discussion:
	//
	//    An R8VEC is a vector of R8's.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    02 July 2005
	//
	//  Author:
	//
	//    John Burkardt
	//
	//  Parameters:
	//
	//    Input, int N, the number of entries in the array.
	//
	//    Input, double R8VEC[N], the array to be checked.
	//
	//    Output, double R8VEC_MIN, the value of the minimum element.
	//
	int i;
	double value;

	value = r8_huge ( );

	if ( n <= 0 )
	{
		return value;
	}

	for ( i = 0; i < n; i++ )
	{
		if(r8vec[i] != missing_value)
			value = min(value, r8vec[i]);
	}
	return value;
}/*}}}*/
//****************************************************************************80
void r8vec_min_max ( int n, double r8vec[], double &min, double&max, double missing_value ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    R8VEC_MIN_MAX returns the value of the minimum and maximum element in an R8VEC.
	//
	//  Discussion:
	//
	//    An R8VEC is a vector of R8's.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    02 July 2005
	//
	//  Author:
	//
	//    John Burkardt
	//
	//  Parameters:
	//
	//    Input, int N, the number of entries in the array.
	//
	//    Input, double R8VEC[N], the array to be checked.
	//
	//    Output, double R8VEC_MIN, the value of the minimum element.
	//
	//    Output, double R8VEC_MAX, the value of the maximum element.
	//
	int i;

	if(n <= 0){
		min = r8_huge();
		max = 0.0;
		return;
	}

	min = r8_huge();
	max = -1e34;

	for ( i = 0; i < n; i++ )
	{
		if(r8vec[i] != missing_value){
			min = std::min(r8vec[i], min);
			max = std::max(r8vec[i], max);
		}
	}

	return;
}/*}}}*/
