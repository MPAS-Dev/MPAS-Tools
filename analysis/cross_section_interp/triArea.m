function area=triArea(A,B,C,R)
% - This function calculates the area of the triangle A,B,C on the
%   surface of a sphere.
%
%   Input: A, B, C
%        A: vertex 1 of triangle
%        B: vertex 2 of triangle
%        C: vertex 3 of triangle
%        R: radius of sphere
%   Output: (returned value area)
%   	area: surface area of triangle on sphere.

R2inv = 1/R/R;

a = acos(dot(B,C)*R2inv);
b = acos(dot(C,A)*R2inv);
c = acos(dot(A,B)*R2inv);

s = 0.5*(a+b+c);

tanqe = sqrt(tan(0.5*s)*tan(0.5*(s-a))*tan(0.5*(s-b))*tan(0.5*(s-c)));

area = abs(4.0*atan(tanqe));

