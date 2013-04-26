#include "Triangle.h"

/*
      Point * points[3];
   public:
      Triangle();
      Triangle(Point& a, Point& b, Point& c);
      ~Triangle();
      void setPoint(Point& p, int n);
      Point * getPoint(int n);
      friend ostream& operator<<(ostream& output, const Triangle& t)
*/

Triangle::Triangle()
{
    points[0] = NULL;
    points[1] = NULL;
    points[2] = NULL;
}


Triangle::Triangle(Point& a, Point& b, Point& c)
{
    points[0] = &a;
    points[1] = &b;
    points[2] = &c;
}


Triangle::~Triangle()
{
    // Nothing to do...
}


void Triangle::setPoint(Point* p, int n)
{
   // assert(n >= 0 && n <= 2);
   points[n] = p;
}


Point * Triangle::getPoint(int n)
{
   // assert(n >= 0 && n <= 2);
   return points[n];
}


ostream& operator<<(ostream& output, const Triangle& t)
{
   // output << "(" << p.x << ", " << p.y << ", " << p.z << ")";
   return output;
}
