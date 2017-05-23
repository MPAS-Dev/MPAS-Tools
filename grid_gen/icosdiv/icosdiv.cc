#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <set>
#include "Point.h"
#include "Triangle.h"

using namespace std;

void add_point(set<Point>& points, Point* newpt, int& np)
{
    set<Point>::iterator ip;

    // If the point doesn't exist, we assign it the next highest number
    //   and add it to the set
    ip = points.find(*newpt);

    if (ip == points.end()) {
        newpt->setNum(np++); 
        points.insert(*newpt);
    }

    // Otherwise, we want newpt to equal the existing point
    else {
        *newpt = *ip;
    }
}

inline bool operator<(Point const& lhs, Point const& rhs)
{
   if (lhs.getX() < rhs.getX())
      return true;
   else if (lhs.getX() > rhs.getX())
      return false;
   else
      if (lhs.getY() < rhs.getY())
         return true;
      else if (lhs.getY() > rhs.getY())
         return false;
      else
         if (lhs.getZ() < rhs.getZ())
            return true;
         else
            return false;
}


inline bool operator>(Point const& lhs, Point const& rhs)
{
   if (lhs.getX() > rhs.getX())
      return true;
   else if (lhs.getX() > rhs.getX())
      return false;
   else
      if (lhs.getY() > rhs.getY())
         return true;
      else if (lhs.getY() > rhs.getY())
         return false;
      else
         if (lhs.getZ() > rhs.getZ())
            return true;
         else
            return false;
}


inline bool operator==(Point const& lhs, Point const& rhs)
{
   if (lhs.getX() == rhs.getX() &&
       lhs.getY() == rhs.getY() &&
       lhs.getZ() == rhs.getZ())
      return true;
   else
      return false;
}


Point * great_circle_points(Point& p1, Point& p2, int n)
{
    double x1, x2, y1, y2, z1, z2;
    double x, y, z;
    double dtheta, dinc, dt, dx;
    Point * pl;
    int i;

    x1 = p1.getX(); y1 = p1.getY(); z1 = p1.getZ();
    x2 = p2.getX(); y2 = p2.getY(); z2 = p2.getZ();

    // For unit sphere, distance is the same as arc angle
    dtheta = p1.distance(p2); 
    dinc = dtheta / (double)(n-1);

    pl = new Point[n];

    pl[0].setXYZ(x1, y1, z1);
    pl[n-1].setXYZ(x2, y2, z2);

    // Fill in interior points
    for(i=1; i<n-1; i++) {
        dt = (double)(i)*dinc;

        if (dt <= 0.5*dtheta) {
            dx = 1.0-tan(0.5*dtheta-dt)/tan(0.5*dtheta);
            x = x1+0.5*dx*(x2-x1);
            y = y1+0.5*dx*(y2-y1);
            z = z1+0.5*dx*(z2-z1);
        }
        else {
            dt = dtheta-dt;
            dx = 1.0-tan(0.5*dtheta-dt)/tan(0.5*dtheta);
            x = x2+0.5*dx*(x1-x2);
            y = y2+0.5*dx*(y1-y2);
            z = z2+0.5*dx*(z1-z2);
        }
        pl[i].setXYZ(x, y, z);
        pl[i].normalize();
    }

    return pl;
}


Point ** divide_triangle(Point& p1, Point& p2, Point& p3, int glevel)
{
    int i, j;
    Point ** line;
    Point * p1p2, * p1p3;
    int n;
        
    line = new Point*[glevel];

    p1p2 = great_circle_points(p1, p2, glevel);
    p1p3 = great_circle_points(p1, p3, glevel);

    line[0] = new Point[1];
    line[1] = new Point[2];

    n = 0;

    line[0][0] = p1;
    line[0][0].setNum(n++); 

    line[1][0] = p1p2[1];
    line[1][0].setNum(n++); 

    line[1][1] = p1p3[1];
    line[1][1].setNum(n++); 

    for(i=2; i<glevel; i++) {
        line[i] = great_circle_points(p1p2[i], p1p3[i], i+1);
        for(j=0; j<i+1; j++)
            line[i][j].setNum(n++);
    }
        
    delete [] p1p2;
    delete [] p1p3;

    return line;
}


int main()
{
    Point icos[12];
    Point ** line;
    Point * p;
    Triangle * t;
    int tri[20][3];
    double x, y, z;
    ifstream fin;
    int i, j, k, np;
    
    set<Point> points;
    vector<Triangle> triangles;
    set<Point>::iterator ip;
    vector<Triangle>::iterator it;

    int div_factor = 76;

    // Read in 12 icosahedral vertices
    fin.open("locs.dat",ifstream::in);
    for(i=0; i<12; i++) {
        fin >> x >> y >> z;
        icos[i].setXYZ(x, y, z);
    }
    fin.close();

    // Read in triangulation of icosahedral points
    fin.open("tri.dat",ifstream::in);
    for(i=0; i<20; i++) {
        fin >> tri[i][0] >> tri[i][1] >> tri[i][2];
    }
    fin.close();
 
    np = 1;

    // In the code below, we actually know which points will be duplicated between
    //   the 20 large (icosahedral) triangles -- exactly those points along the perimeter
    //   of the triangle; so, we could take advantage of this information in the
    //   add_point() subroutine.


    // Subdivide each triangle
    for(k=0; k<20; k++) {
        line = divide_triangle(icos[tri[k][0]-1], icos[tri[k][1]-1], icos[tri[k][2]-1], div_factor);

        // Get triangulation
        for(i=1; i<div_factor; i++) {
            for(j=0; j<i-1; j++) {
//                cout << "Creating triangle from " << line[i][j].getNum() << " " << line[i-1][j].getNum() << " " << line[i][j+1].getNum() << endl;
                t = new Triangle;
                p = new Point; *p = line[i][j];   add_point(points, p, np); t->setPoint(p, 0);
                p = new Point; *p = line[i-1][j]; add_point(points, p, np); t->setPoint(p, 1);
                p = new Point; *p = line[i][j+1]; add_point(points, p, np); t->setPoint(p, 2);
                triangles.push_back(*t);

//                cout << "Creating triangle from " << line[i-1][j].getNum() << " " << line[i-1][j+1].getNum() << " " << line[i][j+1].getNum() << endl;
                t = new Triangle;
                p = new Point; *p = line[i-1][j];   add_point(points, p, np); t->setPoint(p, 0);
                p = new Point; *p = line[i-1][j+1]; add_point(points, p, np); t->setPoint(p, 1);
                p = new Point; *p = line[i][j+1];   add_point(points, p, np); t->setPoint(p, 2);
                triangles.push_back(*t);
            }
//            cout << "Creating triangle from " << line[i][i-1].getNum() <<  " " << line[i-1][i-1].getNum() << " " << line[i][i].getNum() <<endl;
            t = new Triangle;
            p = new Point; *p = line[i][i-1];   add_point(points, p, np); t->setPoint(p, 0);
            p = new Point; *p = line[i-1][i-1]; add_point(points, p, np); t->setPoint(p, 1);
            p = new Point; *p = line[i][i];     add_point(points, p, np); t->setPoint(p, 2);
            triangles.push_back(*t);
        }
        p = NULL; t = NULL;

        for(j=0; j<div_factor; j++)
            delete [] line[j];
        delete [] line;
    }

    for(ip = points.begin(); ip != points.end(); ip++) {
       cout << *ip << endl;
    }

/*
    for(it = triangles.begin(); it != triangles.end(); it++) {
       cout << *(it->getPoint(0)) << endl;
       cout << *(it->getPoint(1)) << endl;
       cout << *(it->getPoint(2)) << endl;
       cout << *(it->getPoint(0)) << endl;
       cout << endl;
       cout << endl;
    }
*/

    triangles.clear();
    points.clear();


    return 0;
}
