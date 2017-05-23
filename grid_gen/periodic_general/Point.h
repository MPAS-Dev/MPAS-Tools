#ifndef _PointH
#define _PointH
#include <iostream>
#include <math.h>
using namespace std;

class Point
{
	private:
		double x, y;
		int boundary_point;
		int num;
	public:
		Point();
		Point(double x, double y, int boundary_point);
		~Point();
		void setX(double x);
		void setY(double y);
		void setXY(double x, double y);
		void setBoundaryPoint(int boundary_point);
		void setNum(int n);
		double getX() const;
		double getY() const;
		double distance(Point& p);
		int isBoundaryPoint() const;
		int getNum() const;
		Point operator+(Point p);
		Point operator-(Point p);
		Point operator*(double s);
		friend ostream& operator<<(ostream& output, const Point& p);
		friend bool operator<(Point const& lhs, Point const& rhs);
};
#endif
