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
  double theta = pi/4.0;
  double d = 1.0;
  double h = sqrt(2.0 * d*d * (1.0 - sin(theta)));
  double y, z;
  
  cout << "picture capillary; capillary := nullpicture;\n";
  cout << "addto capillary contour " << mpoint(0.0,h);
  
  for(int i=99;i>2;i--) {
    z = (i/100.0) * h;
    y = capillary(z,h,d);
    cout << " .. " << mpoint(y,z);
  }
  
  cout << " -- " << mpoint(y,-0.5);
  cout << " -- " << mpoint(0.0,-0.5);  
  cout << " -- cycle;\n";
}
