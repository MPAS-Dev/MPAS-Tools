#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

class pnt {/*{{{*/
	public:
		double x, y, z;
		double lat, lon;
		int idx;
		bool positiveLonRange;

		pnt(double x_, double y_, double z_, int idx_) {
			(*this).x = x_;
			(*this).y = y_;
			(*this).z = z_;
			(*this).idx = idx_;
			(*this).positiveLonRange = true;
			(*this).buildLat();
			(*this).buildLon();
		}

		pnt(double x_, double y_, double z_) {
			(*this).x = x_;
			(*this).y = y_;
			(*this).z = z_;
			(*this).idx = 0;
			(*this).positiveLonRange = true;
			(*this).buildLat();
			(*this).buildLon();
		}

		pnt() {
			(*this).x = 0.0;
			(*this).y = 0.0;
			(*this).z = 0.0;
			(*this).idx = 0;
			(*this).positiveLonRange = true;
			(*this).lat = 0.0;
			(*this).lon = 0.0;
		}

		friend pnt operator*(const double d, const pnt &p);
		friend std::ostream & operator<<(std::ostream &os, const pnt &p);
		friend std::istream & operator>>(std::istream &is, pnt &p);

		pnt& operator=(const pnt &p){/*{{{*/
			x = p.x;
			y = p.y;
			z = p.z;
			idx = p.idx;
			(*this).buildLat();
			(*this).buildLon();
			return *this;
		}/*}}}*/
		bool operator==(const pnt &p) const {/*{{{*/
			return (x == p.x) & (y == p.y) & (z == p.z);
		}/*}}}*/
		pnt operator-(const pnt &p) const {/*{{{*/
			double x_, y_, z_;

			x_ = x-p.x;
			y_ = y-p.y;
			z_ = z-p.z;

			return pnt(x_,y_,z_,0);
		}/*}}}*/
		pnt operator+(const pnt &p) const {/*{{{*/
			double x_, y_, z_;

			x_ = x+p.x;
			y_ = y+p.y;
			z_ = z+p.z;

			return pnt(x_,y_,z_,0);
		}/*}}}*/
		pnt operator*(double d) const {/*{{{*/
			double x_, y_, z_;
			x_ = x*d;
			y_ = y*d;
			z_ = z*d;
			return pnt(x_,y_,z_,0);
		}/*}}}*/
		pnt operator/(double d) const {/*{{{*/
			double x_, y_, z_;

			if(d == 0.0){
				std::cout << "pnt: operator/" << std::endl;
				std::cout << (*this) << std::endl;
			}

			assert(d != 0.0);
			x_ = x/d;
			y_ = y/d;
			z_ = z/d;
			return pnt(x_,y_,z_,0);
		}/*}}}*/
		pnt& operator/=(double d){/*{{{*/
			if(d == 0.0){
				std::cout << "pnt: operator /=" << std::endl << (*this) << std::endl;
			}
			assert(d != 0.0);
			x = x/d;
			y = y/d;
			z = z/d;
			return *this;
		}/*}}}*/
		pnt& operator+=(const pnt &p){/*{{{*/
			x += p.x;
			y += p.y;
			z += p.z;
			return *this;
		}/*}}}*/
		double operator[](int i) const {/*{{{*/
			if(i == 0){
				return x;
			} else if(i == 1){
				return y;
			} else {
				return z;
			}
		}/*}}}*/
		void normalize(){/*{{{*/
			double norm;

			norm = x*x + y*y + z*z;
			if(norm == 0){
				std::cout << "pnt: normalize" << std::endl;
				std::cout << x << " " << y << " " << z << " " << idx << std::endl;

				assert(norm != 0);
			}	
			norm = sqrt(norm);

			x = x/norm;
			y = y/norm;
			z = z/norm;
		}/*}}}*/
		double dot(const pnt &p) const {/*{{{*/
			double junk;
			junk = x*p.x+y*p.y+z*p.z;

			return junk;
		}/*}}}*/
		double dotForAngle(const pnt &p) const {/*{{{*/
			double junk;
			junk = x*p.x+y*p.y+z*p.z;
			if(junk > 1.0){
				junk = 1.0;
			}

			if(junk < -1.0){
				junk = -1.0;
			}
			return acos(junk);
		}/*}}}*/
		pnt cross(const pnt &p) const {/*{{{*/
			double x_, y_, z_;

			x_ = y*p.z - p.y*z;
			y_ = z*p.x - p.z*x;
			z_ = x*p.y - p.x*y;

			return pnt(x_,y_,z_,0);
		}/*}}}*/
		void fixPeriodicity(const pnt &p, const double xRef, const double yRef){/*{{{*/
			/* The fixPeriodicity function fixes the periodicity of the current point relative
			 *     to point p. It should only be used on a point in the x/y plane.
			 *     xRef and yRef should be the extents of the non-periodic plane.
			 */

			pnt dist_vec;

			dist_vec = (*this) - p;

			if(fabs(dist_vec.x) > xRef * 0.6){
#ifdef _DEBUG
				std::cout << "   Fixing x periodicity " << endl;
#endif
				(*this).x += -(dist_vec.x/fabs(dist_vec.x)) * xRef;
			}
			if(fabs(dist_vec.y) > yRef * 0.6){
#ifdef _DEBUG
				std::cout << "   Fixing y periodicity " << endl;
#endif
				(*this).y += -(dist_vec.y/fabs(dist_vec.y)) * yRef;
			}

		}/*}}}*/
		double magnitude() const{/*{{{*/
			return sqrt(x*x + y*y + z*z);
		}/*}}}*/
		double magnitude2() const {/*{{{*/
			return x*x + y*y + z*z;
		}/*}}}*/
		void rotate(const pnt &vec, const double angle){/*{{{*/
			double x_, y_, z_;
			double cos_t, sin_t;

			cos_t = cos(angle);
			sin_t = sin(angle);

			x_ = (cos_t + vec.x*vec.x *(1.0-cos_t)) * x
				+(vec.x*vec.y*(1.0-cos_t) - vec.z * sin_t) * y
				+(vec.x*vec.z*(1.0-cos_t) + vec.y * sin_t) * z;

			y_ = (vec.y*vec.x*(1.0-cos_t) + vec.z*sin_t) * x
				+(cos_t + vec.y*vec.x*(1.0 - cos_t)) * y
				+(vec.y*vec.z*(1.0-cos_t) - vec.x*sin_t) * z;

			z_ = (vec.z*vec.x*(1.0-cos_t)-vec.y*sin_t)*x
				+(vec.z*vec.y*(1.0-cos_t)-vec.x*sin_t)*y
				+(cos_t + vec.x*vec.x*(1.0-cos_t))*z;

			x = x_;
			y = y_;
			z = z_;
		}/*}}}*/
		void setPositiveLonRange(bool value){/*{{{*/
			(*this).positiveLonRange = value;
			(*this).buildLon();
		}/*}}}*/
		double getLat() const {/*{{{*/
			return (*this).lat;
		}/*}}}*/
		double getLon() const {/*{{{*/
			return (*this).lon;
		}/*}}}*/
		void buildLat() {/*{{{*/
			double dl;

			dl = sqrt((*this).x*(*this).x + (*this).y*(*this).y + (*this).z*(*this).z);
			(*this).lat = asin(z/dl);
		}/*}}}*/
		void buildLon() {/*{{{*/
			double lon;

			lon = atan2(y, x);

			// If the prime meridian is the minimum, this means the degree
			// range is 0-360, so we need to translate.
			if ( (*this).positiveLonRange ) {
				if(lon < 0.0) {
					lon = 2.0*M_PI + lon;
				}
			}

			(*this).lon = lon;
		}/*}}}*/
		double sphereDistance(const pnt &p) const {/*{{{*/
			double arg;

			arg = sqrt( powf(sin(0.5*((*this).getLat()-p.getLat())),2) +
					cos(p.getLat())*cos((*this).getLat())*powf(sin(0.5*((*this).getLon()-p.getLon())),2));
			return 2.0*asin(arg);

			/*
			pnt temp, cross;
			double dot;
			temp.x = x;
			temp.y = y;
			temp.z = z;

			cross = temp.cross(p);
			dot = temp.dot(p);



			return atan2(cross.magnitude(),dot);
			// */
		}/*}}}*/
	struct hasher {/*{{{*/
		size_t operator()(const pnt &p) const {
			uint32_t hash; 
			size_t i, key[3] = { (size_t)p.x, (size_t)p.y, (size_t)p.z };
			for(hash = i = 0; i < sizeof(key); ++i) {
				hash += ((uint8_t *)key)[i];
				hash += (hash << 10);
				hash ^= (hash >> 6);
			}
			hash += (hash << 3);
			hash ^= (hash >> 11);
			hash += (hash << 15);
			return hash;
		}
	};/*}}}*/
	struct idx_hasher {/*{{{*/
		size_t operator()(const pnt &p) const {
			return (size_t)p.idx;
		}
	};/*}}}*/
};/*}}}*/

inline pnt operator*(const double d, const pnt &p){/*{{{*/
	return pnt(d*p.x, d*p.y, d*p.z, 0);
}/*}}}*/

inline std::ostream & operator<<(std::ostream &os, const pnt &p){/*{{{*/
	os << '(' << p.x << ", " << p.y << ", " << p.z << ")  idx=" << p.idx << ' ';
	//os << p.x << " " << p.y << " " << p.z;
	return os;
}/*}}}*/
inline std::istream & operator>>(std::istream &is, pnt &p){/*{{{*/
	//is >> p.x >> p.y >> p.z >> p.idx;
	is >> p.x >> p.y >> p.z;
	return is;
}/*}}}*/

/* Geometric utility functions {{{ */
pnt gcIntersect(const pnt &c1, const pnt &c2, const pnt &v1, const pnt &v2){/*{{{*/
	/*
	 * gcIntersect is intended to compute the intersection of two great circles.
	 *   The great circles pass through the point sets c1-c2 and v1-v2.
	 */

	/*
	pnt n, m,c ;
	double dot;


	//n = c1 cross c2
	n = c1.cross(c2);

	//m = v1 cross v2
	m = v1.cross(v2);

	//c = n cross m
	c = n.cross(m);

	dot = c1.dot(c);

	dot = dot/fabs(dot);

	c = c*dot;

	c.normalize();

	return c;
	// */

	double n1, n2, n3;
	double m1, m2, m3;
	double xc, yc, zc;
	double dot;
	pnt c;

	n1 =  (c1.y * c2.z - c2.y * c1.z);
	n2 = -(c1.x * c2.z - c2.x * c1.z);
	n3 =  (c1.x * c2.y - c2.x * c1.y);

	m1 =  (v1.y * v2.z - v2.y * v1.z);
	m2 = -(v1.x * v2.z - v2.x * v1.z);
	m3 =  (v1.x * v2.y - v2.x * v1.y);

	xc =  (n2 * m3 - n3 * m2);
	yc = -(n1 * m3 - n3 * m1);
	zc =  (n1 * m2 - n2 * m1);

	dot = c1.x*xc + c1.y*yc + c1.z*zc;

	if (dot < 0.0) {
		xc = -xc;
		yc = -yc;
		zc = -zc;
	}

	c.x = xc;
	c.y = yc;
	c.z = zc;

	c.normalize();
	return c;
}/*}}}*/
pnt planarIntersect(const pnt &c1, const pnt &c2, const pnt &v1, const pnt &v2){/*{{{*/
	/*
	 * planarIntersect is intended to compute the point of intersection
	 *    of the lines c1-c2 and v1-v2 in a plane.
	 */
	double alpha_numerator, beta_numerator, denom;
	double alpha, beta;
	pnt c;

	alpha_numerator = (v2.x - v1.x) * (c1.y - v2.y) - (v2.y - v1.y) * (c1.x - v2.x);
//	beta_numerator = (c1.x - v2.x) * (c2.y - c1.y) - (c1.y - v2.y) * (c2.x - c1.x);
	denom = (c2.x - c1.x) * (v2.y - v1.y) - (c2.y - c1.y) * (v2.x - v1.x);

	alpha = alpha_numerator / denom;
//	beta = beta_numerator / demon;

	c = alpha * (v2 - v1) + v1;

	return c;
}/*}}}*/

double planarTriangleArea(const pnt &A, const pnt &B, const pnt &C){/*{{{*/
	/*
	 * planarTriangleArea uses Heron's formula to compute the area of a triangle in a plane.
	 */
	pnt ab, bc, ca;
	double s, a, b, c, e;

	ab = B - A;
	bc = C - B;
	ca = A - C;

	// Get side lengths
	a = ab.magnitude();
	b = bc.magnitude();
	c = ca.magnitude();

	// Semi-perimeter of TRI(ABC)
	s = (a + b + c) * 0.5;

	return sqrt(s * (s - a) * (s - b) * (s - c));
}/*}}}*/
double sphericalTriangleArea(const pnt &A, const pnt &B, const pnt &C){/*{{{*/
	/*
	 * sphericalTriangleArea uses the spherical analog of Heron's formula to compute
	 *    the area of a triangle on the surface of a sphere.
	 *
	 */
	double tanqe, s, a, b, c, e;

	a = A.sphereDistance(B);
	b = B.sphereDistance(C);
	c = C.sphereDistance(A);
	s = 0.5*(a+b+c);

	tanqe = sqrt(tan(0.5*s)*tan(0.5*(s-a))*tan(0.5*(s-b))*tan(0.5*(s-c)));
	e = 4.*atan(tanqe);

	return e;
}/*}}}*/

double planeAngle(const pnt &A, const pnt &B, const pnt &C, const pnt &n){/*{{{*/
	/*
	 * planeAngle computes the angles between the vectors AB and AC in a plane defined with
	 *    the normal vector n
	 */

	/*
	pnt ab;
	pnt ac;
	pnt cross;
	double cos_angle;

	ab = B - A;
	ac = C - A;
	cross = ab.cross(ac);
	cos_angle = ab.dot(ac)/(ab.magnitude() * ac.magnitude());

	cos_angle = std::min(std::max(cos_angle,-1.0),1.0);

	if(cross.x*n.x + cross.y*n.y + cross.z*n.z >= 0){
		return acos(cos_angle);
	} else {
		return -acos(cos_angle);
	}
	// */

	double ABx, ABy, ABz, mAB;
	double ACx, ACy, ACz, mAC;
	double Dx, Dy, Dz;
	double cos_angle;

	ABx = B.x - A.x;
	ABy = B.y - A.y;
	ABz = B.z - A.z;
	mAB = sqrt(ABx*ABx + ABy*ABy + ABz*ABz);

	ACx = C.x - A.x;
	ACy = C.y - A.y;
	ACz = C.z - A.z;
	mAC = sqrt(ACx*ACx + ACy*ACy + ACz*ACz);


	Dx =   (ABy * ACz) - (ABz * ACy);
	Dy = -((ABx * ACz) - (ABz * ACx));
	Dz =   (ABx * ACy) - (ABy * ACx);

	cos_angle = (ABx*ACx + ABy*ACy + ABz*ACz) / (mAB * mAC);

	if (cos_angle < -1.0) {
		cos_angle = -1.0;
	} else if (cos_angle > 1.0) {
		cos_angle = 1.0;
	}

	if ((Dx*n.x + Dy*n.y + Dz*n.z) >= 0.0) {
		return acos(cos_angle);
	} else {
		return -acos(cos_angle);
	}
}/*}}}*/
pnt pntFromLatLon(const double &lat, const double &lon){/*{{{*/
	/*
	 * pntFromLatLon constructs a point location in 3-Space
	 *    from an initial latitude and longitude location.
	 */
	pnt temp;
	temp.x = cos(lon) * cos(lat);
	temp.y = sin(lon) * cos(lat);
	temp.z = sin(lat);
	temp.normalize();
	temp.buildLat();
	temp.buildLon();
	return temp;
}/*}}}*/
/*}}}*/
