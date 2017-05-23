#include "Point.h"

Point::Point()
{
   x = 0.0;
   y = 0.0;
   z = 0.0;
   num = 0;
}


Point::Point(double x, double y, double z)
{
   this->x = x;
   this->y = y;
   this->z = z;
}


Point::~Point()
{

}


void Point::setX(double x)
{
   this->x = x;
}


void Point::setY(double y)
{
   this->y = y;
}


void Point::setZ(double z)
{
   this->z = z;
}


void Point::setXYZ(double x, double y, double z)
{
   this->x = x;
   this->y = y;
   this->z = z;
}


void Point::setNum(int n)
{
   num = n;
}


double Point::getX() const
{
   return x;
}


double Point::getY() const
{
   return y;
}


double Point::getZ() const
{
   return z;
}


double Point::distance(Point& p)
{
   // Assume we're on the unit sphere
   return acos(p.getX()*x + p.getY()*y + p.getZ()*z);
}


int Point::getNum() const
{
   return num;
}


void Point::normalize()
{
   double mag;

   mag = sqrt(x*x + y*y + z*z);
   x = x / mag;
   y = y / mag;
   z = z / mag;
}


Point Point::operator+(Point p)
{
   Point retval;

   retval.x = x + p.x;
   retval.y = y + p.y;
   retval.z = z + p.z;
   retval.num = num;

   return retval;
}


Point Point::operator-(Point p)
{
   Point retval;

   retval.x = x - p.x;
   retval.y = y - p.y;
   retval.z = z - p.z;
   retval.num = num;

   return retval;
}


Point Point::operator*(double s)
{
   Point retval;

   retval.x = s * x;
   retval.y = s * y;
   retval.z = s * z;
   retval.num = num;

   return retval;
}


ostream& operator<<(ostream& output, const Point& p)
{
   output << p.num << " : " << p.x << " " << p.y << " " << p.z;
   // output << p.x << " " << p.y << " " << p.z;
   return output;
}
