#include <vector>
#include <iostream>
#include <fstream>

class pnt {/*{{{*/
	public:
		double x, y, z;
		int isBdry;
		int idx;
		int vert_idx1, vert_idx2;

		pnt(double x_, double y_, double z_, int isBdry_, int idx_)
			:  x(x_), y(y_), z(z_), isBdry(isBdry_), idx(idx_) {	}

		pnt(double x_, double y_, double z_)
			: x(x_), y(y_), z(z_), isBdry(0), idx(0) { }

		pnt(double x_, double y_, double z_, int isBdry_)
			: x(x_), y(y_), z(z_), isBdry(isBdry_), idx(0) { }

		pnt()
			: x(0.0), y(0.0), z(0.0), isBdry(0), idx(0) { }

		friend pnt operator*(const double d, const pnt &p);
		friend std::ostream & operator<<(std::ostream &os, const pnt &p);
		friend std::istream & operator>>(std::istream &is, pnt &p);

		pnt& operator=(const pnt &p){/*{{{*/
			x = p.x;
			y = p.y;
			z = p.z;
			isBdry = p.isBdry;
			idx = p.idx;
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

			return pnt(x_,y_,z_,0,0);
		}/*}}}*/
		pnt operator+(const pnt &p) const {/*{{{*/
			double x_, y_, z_;

			x_ = x+p.x;
			y_ = y+p.y;
			z_ = z+p.z;

			return pnt(x_,y_,z_,0,0);
		}/*}}}*/
		pnt operator*(double d) const {/*{{{*/
			double x_, y_, z_;
			x_ = x*d;
			y_ = y*d;
			z_ = z*d;
			return pnt(x_,y_,z_,0,0);
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
			return pnt(x_,y_,z_,0,0);
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
				std::cout << x << " " << y << " " << z << " " << isBdry << " " << idx << std::endl;

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

			return pnt(x_,y_,z_,0,0);
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
		double getLat() const {/*{{{*/
			double dl;

			dl = sqrt((*this).x*(*this).x + (*this).y*(*this).y + (*this).z*(*this).z);
			return asin(z/dl);
		}/*}}}*/
		double getLon() const {/*{{{*/
			double lon;

			lon = atan2(y,x);

			if(lon < -M_PI){
				lon = lon + M_2_PI;
			} else if (lon > M_PI) {
				lon = lon - M_2_PI;
			}

			return lon+M_PI;

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
	struct edge_hasher {/*{{{*/
		size_t operator()(const pnt &p) const {
			uint32_t hash; 
			size_t i, key[2] = { (size_t)p.vert_idx1, (size_t)p.vert_idx2 };
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
};/*}}}*/
class tri {/*{{{*/
	public:
	int vi1, vi2, vi3;
	int ei1, ei2, ei3;
	int idx;

	tri() : vi1(0), vi2(0), vi3(0), idx(0) { }

	tri(int vi1_, int vi2_, int vi3_)
		: vi1(vi1_), vi2(vi2_), vi3(vi3_), idx(0) { }

	tri(int vi1_, int vi2_, int vi3_, int idx_)
		: vi1(vi1_), vi2(vi2_), vi3(vi3_), idx(idx_) { }


	friend std::ostream & operator<<(std::ostream &os, const tri &t);
	friend std::istream & operator>>(std::istream &is, tri &t);

	tri sortedTri(){/*{{{*/
		int v1, v2, v3, swp_v;
		v1 = vi1;
		v2 = vi2;
		v3 = vi3;

		if(v1 > v2){
			swp_v = v1;
			v1 = v2;
			v2 = swp_v;
		}

		if(v2 > v3){
			swp_v = v2;
			v2 = v3;
			v3 = swp_v;
		}

		if(v1 > v2){
			swp_v = v1;
			v1 = v2;
			v2 = swp_v;
		}

		return tri(v1,v2,v3);
	}/*}}}*/
	bool operator==(const tri &t) const {/*{{{*/
		return (vi1 == t.vi1) & (vi2 == t.vi2) & (vi3 == t.vi3);
	}/*}}}*/
	struct hasher {/*{{{*/
		size_t operator()(const tri &t) const {
			uint32_t hash; 
			size_t i, key[3] = { t.vi1, t.vi2, t.vi3 };
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
};/*}}}*/

inline pnt operator*(const double d, const pnt &p){/*{{{*/
	return pnt(d*p.x, d*p.y, d*p.z, 0, 0);
}/*}}}*/

inline std::ostream & operator<<(std::ostream &os, const pnt &p){/*{{{*/
	os << '(' << p.x << ", " << p.y << ", " << p.z << ") bdry=" << p.isBdry << " idx=" << p.idx << ' ';
	//os << p.x << " " << p.y << " " << p.z;
	return os;
}/*}}}*/
inline std::istream & operator>>(std::istream &is, pnt &p){/*{{{*/
	//is >> p.x >> p.y >> p.z >> p.isBdry >> p.idx;
	is >> p.x >> p.y >> p.z;
	return is;
}/*}}}*/

inline std::ostream & operator<<(std::ostream &os, const tri &t){/*{{{*/
	//return os << t.vi1 << " " << t.vi2 << " " << t.vi3;
	return os << t.vi1 << " " << t.vi2 << " " << t.vi3;
}/*}}}*/
inline std::istream & operator>>(std::istream &is, tri &t){/*{{{*/
	//return is >> t.vi1 >> t.vi2 >> t.vi3;
	return is >> t.vi1 >> t.vi2 >> t.vi3;
}/*}}}*/

void circumcenter(const pnt &A,const pnt &B,const pnt &C, pnt &cent){/*{{{*/
	double a, b, c;
	double pbc, apc, abp;

	a = (B-C).magnitude2();
	b = (C-A).magnitude2();
	c = (A-B).magnitude2();

	pbc = a*(-a + b + c);
	apc = b*( a - b + c);
	abp = c*( a + b - c);

	cent = (pbc*A + apc*B + abp*C)/(pbc + apc + abp);
}/*}}}*/
double circumradius(const pnt &A, const pnt &B, const pnt &C){/*{{{*/

	pnt ccenter;  
	pnt ac;

	circumcenter(A,B,C,ccenter);
	ac = A - ccenter;

	return ac.magnitude();

/*
	pnt a(0.0,0.0,0.0,0);
	pnt b(0.0,0.0,0.0,0);
	pnt sub(0.0,0.0,0.0,0);
	pnt cross(0.0,0.0,0.0,0);
	double top, bot;

	a = A-C;
	b = B-C;

	sub = a-b;
	cross = a.cross(b);

	top = a.magnitude() * b.magnitude() * sub.magnitude();
	bot = 2.0 * cross.magnitude();

	return top/bot; // */

/*	dx = a.x - b.x;
	dy = a.y - b.y;
	dz = a.z - b.z;

	sa = sqrt(dx*dx + dy*dy + dz*dz);

	dx = b.x - c.x;
	dy = b.y - c.y;
	dz = b.z - c.z;

	sb = sqrt(dx*dx + dy*dy + dz*dz);

	dx = c.x - a.x;
	dy = c.y - a.y;
	dz = c.z - a.z;

	sc = sqrt(dx*dx + dy*dy + dz*dz);

	dx = (sa+sb+sc)/2.0; // Semiperimeter
	dy = 4.0*sqrt(dx*(sa+sb-dx)*(sa+sc-dx)*(sb+sc-dx)); // Bottom of circumradius computation

	return (sa*sb*sc)/dy;*/
}/*}}}*/
double triAreaBAK(const pnt &A, const pnt &B, const pnt &C){/*{{{*/
	/**************************************************************************
	 * - This function calculates the area of the triangle A,B,C on the
	 *   surface of a sphere.
	 *
	 *   Input: A, B, C
	 *        A: vertex 1 of triangle
	 *        B: vertex 2 of triangle
	 *        C: vertex 3 of triangle
	 *   Output: (returned value area)
	 *   	area: surface area of triangle on sphere.
	 **************************************************************************/
	pnt u12, u23, u31;
	double a, b, c, s, tanqe, area;	
	double sign;

	area = 0.0;

	//Compute Surface normal for triangle to get "sign" of area
	u12 = B - A;
	u23 = C - A;

	u31 = u12.cross(u23);

	u31.normalize();
	sign = u31.magnitude();
	assert(sign != 0.0);

	//dot the surface norm with one of the points to get sign.
	sign = u31.dot(A);
	sign = sign/fabs(sign);

	a = A.dotForAngle(B);
	b = B.dotForAngle(C);
	c = C.dotForAngle(A);

	s = 0.5*(a+b+c);

	tanqe = sqrt(tan(0.5*s)*tan(0.5*(s-a))*tan(0.5*(s-b))*tan(0.5*(s-c)));

//	area = sign*4.0*atan(tanqe);
	area = 4.0*atan(tanqe);

	if(isnan(area))
		area = 0.0;

	return area;
}/*}}}*/
double triArea(const pnt &A, const pnt &B, const pnt &C){/*{{{*/
	double tanqe, s, a, b, c, e;

	a = A.sphereDistance(B);
	b = B.sphereDistance(C);
	c = C.sphereDistance(A);
	s = 0.5*(a+b+c);

	tanqe = sqrt(tan(0.5*s)*tan(0.5*(s-a))*tan(0.5*(s-b))*tan(0.5*(s-c)));
	e = 4.*atan(tanqe);

	return e;
}/*}}}*/
int isCcw(const pnt &p1, const pnt &p2, const pnt &p3){/*{{{*/
	// */
	double sign;
	pnt a, b, cross;

	a = p2-p1;
	b = p3-p1;
	cross = a.cross(b);

	sign = cross.dot(p1);

	if(sign > 0){
		return 1;
	} else {
		return 0;
	}// */

	/*
	double side_check;
	pnt a, cross;

	a = p3 - p1;
	a.normalize();
	cross = p1.cross(a);
	side_check = cross.dot(p2);	

	if(side_check < 0){
		return 1;
	} else {
		return 0;
	} // */
}/*}}}*/
pnt gcIntersect(const pnt &c1, const pnt &c2, const pnt &v1, const pnt &v2){/*{{{*/
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

double planeAngle(const pnt &A, const pnt &B, const pnt &C, const pnt &n){/*{{{*/
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
	pnt temp;
	temp.x = cos(lon) * cos(lat);
	temp.y = sin(lon) * cos(lat);
	temp.z = sin(lat);
	temp.normalize();
	return temp;
}/*}}}*/
bool flip_vertices(const pnt &c1, const pnt &c2, const pnt &v1, const pnt &v2){/*{{{*/
	pnt d_c, d_v;
	double ci, cj, ck;

	d_c = c2 - c1;
	d_v = v2 - v1;
	
	ci = d_c.y*d_v.z - d_c.z*d_v.y;
	cj = d_c.z*d_v.x - d_c.x*d_v.z;
	ck = d_c.x*d_v.y - d_c.y*d_v.x;

	if ((ci*c1.x + cj*c1.y + ck*c1.z) >= 0.0) {
		return false;
	} else {
		return true;
	}
}/*}}}*/
