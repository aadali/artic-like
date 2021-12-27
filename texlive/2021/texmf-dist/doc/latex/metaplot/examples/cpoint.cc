/* MetaContour, version 0.1 -------------------------------------- */
/* Copyright(C) 2004, Brooks Moses                                 */
/*                                                                 */
/* This version of MetaContour is made available under the Gnu     */
/* Public License; see metacontour_main.cc for details.            */
/*                                                                 */
/* This is a very pre-release version of MetaContour, distributed  */
/* primarily as an example of a use of MetaPlot.  It can be        */
/* compiled with gcc using the line:                               */
/*                                                                 */
/*   g++ metacontour.cc cpoint.cc metacontour_main.cc              */
/*                                                                 */
/*-----------------------------------------------------------------*/

/* Class for 2D points with one data value each */  
  
#include<sstream>
#include"cpoint.h"
using namespace std;


//-----cpoint--------------------------------------------------

cpoint::cpoint() {
  x = 0;
  y = 0;
  z = 0;
}

cpoint::cpoint(double xx, double yy, double zz) {
  x = xx;
  y = yy;
  z = zz;
}

cpoint cpoint::interpolateto(cpoint x2, double z3) {
  double w = (z3-z)/(x2.z-z);
  if (w > 1) w=1;
  if (w < 0) w=0;
  double x3 = x + (x2.x-x)*w;
  double y3 = y + (x2.y-y)*w;
  return cpoint(x3,y3,z3);
}

cpoint cpoint::partwayto(cpoint x2, double w) {
  if (w > 1) w=1;
  if (w < 0) w=0;
  double x3 = x + (x2.x-x)*w;
  double y3 = y + (x2.y-y)*w;
  double z3 = z + (x2.z-z)*w;
  return cpoint(x3,y3,z3);
}

// Derivation:
//   F = ax + by + cxy + d     x,y in [0,1]
//   Fx = a + by, Fy = b + cx
//   Thus: saddle point S is at (-b/c, -a/c)
//   Plugging this in, and solving, we get
//   F(S) = (f00*f11 - f01*f10) / (f00 + f11 - f01 - f10)
//
//   Note: mapping x1 = (0,0), x2 = (1,0), x3 = (1,1), x4 = (0,1)
//
double cpoint::saddlepoint(cpoint x1, cpoint x2, cpoint x3, cpoint x4) {
  return (x1.z * x3.z - x2.z * x4.z) / (x1.z + x3.z - x2.z - x4.z);
}

string cpoint::metapoint() {
  ostringstream pointstring;
  pointstring.setf(ios_base::fixed,ios_base::floatfield);
  pointstring.precision(5);
  pointstring << '(' << x << ',' << y << ')';  
  return pointstring.str();
}


//-----cparray-------------------------------------------------

cparray::cparray(int x, int y) {
  if (x<1) x=1;
  if (y<1) y=1;
  i = x;
  j = y;
  a = new valarray<cpoint>(x*y);
}

cparray::~cparray() {
  delete a;
}

cpoint& cparray::operator()(int x, int y)
{
  return (*a)[x*j+y];
}
