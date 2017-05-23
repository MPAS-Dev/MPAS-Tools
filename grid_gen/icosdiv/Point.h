#ifndef _PointH
#define _PointH
#include <iostream>
#include <math.h>
using namespace std;

class Point
{
   private:
      double x, y, z;
      int num;
   public:
      Point();
      Point(double x, double y, double z);
      ~Point();
      void setX(double x);
      void setY(double y);
      void setZ(double z);
      void setXYZ(double x, double y, double z);
      void setNum(int n);
      double getX() const;
      double getY() const;
      double getZ() const;
      double distance(Point& p);
      int getNum() const;
      void normalize();
      Point operator+(Point p);
      Point operator-(Point p);
      Point operator*(double s);
      friend ostream& operator<<(ostream& output, const Point& p);
};

#endif
