/* Simple capillary surface */  
  
#include<sstream>
#include<iostream>
#include<cmath>

using namespace std;

string mpoint(double x, double y) {
  ostringstream pointstring;
  pointstring.setf(ios_base::fixed,ios_base::floatfield);
  pointstring.precision(5);
  pointstring << '(' << x << ',' << y << ')';  
  return pointstring.str();
}

double arccosh(double x) {
  return log(x + sqrt(x*x-1));
}

double capillary(double z, double h, double d) {
  return arccosh(2.0*d/z) - arccosh(2.0*d/h)
         + sqrt(4.0 - h*h/(d*d)) - sqrt(4.0 - z*z/(d*d));
}

const double pi = 3.1415926535897932384626433832795;

int main() {
  double theta = pi/6.0;
  double d = 1.0;
  double h = sqrt(2.0 * d*d * (1.0 - sin(theta)));
  double y, z;

  double zmax = h; double zmin = -0.5;
  double ymax = capillary(0.03*h,h,d);
  
  cout << "picture capillary.fplot; capillary.fplot := nullpicture;\n";
  cout << "addto capillary.fplot contour " << mpoint(0.0,1.0);
  for(int i=99;i>2;i--) {
    z = (i/100.0) * h;
    y = capillary(z,h,d);
    cout << " .. " << mpoint(y/ymax,(z-zmin)/(zmax-zmin));
  }
  cout << " -- " << mpoint(1.0,0.0);
  cout << " -- " << mpoint(0.0,0.0);  
  cout << " -- cycle;\n\n";

  cout << "picture capillary.lplot; capillary.lplot := nullpicture;\n";
  cout << "addto capillary.lplot doublepath " << mpoint(0.0,1.0);
  for(int i=99;i>2;i--) {
    z = (i/100.0) * h;
    y = capillary(z,h,d);
    cout << " .. " << mpoint(y/ymax,(z-zmin)/(zmax-zmin));
  }
  cout << ";\n\n";

  cout.setf(ios_base::fixed,ios_base::floatfield);
  cout.precision(5);
  cout << "numeric capillary.xleft; ";
  cout << "capillary.xleft = 0.0;\n";
  cout << "numeric capillary.xright; ";
  cout << "capillary.xright = " << ymax << ";\n";
  cout << "numeric capillary.ybot; ";
  cout << "capillary.ybot = " << zmin << ";\n";
  cout << "numeric capillary.ytop; ";
  cout << "capillary.ytop = " << zmax << ";\n\n";
  
  cout << "pair capillary.contactpoint; ";
  cout << "capillary.contactpoint = (0.0, 1.0);\n";
  cout << "numeric capillary.contactangle; ";
  cout << "capillary.contactangle = " << theta * 180/pi << ";\n";
}
