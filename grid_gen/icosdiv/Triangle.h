#ifndef _TriangleH
#define _TriangleH
#include <iostream>
#include "Point.h"
using namespace std;

class Triangle
{
   private:
      Point * points[3];
   public:
      Triangle();
      Triangle(Point& a, Point& b, Point& c);
      ~Triangle();
      void setPoint(Point* p, int n);
      Point * getPoint(int n);
      friend ostream& operator<<(ostream& output, const Triangle& t);
};
#endif
