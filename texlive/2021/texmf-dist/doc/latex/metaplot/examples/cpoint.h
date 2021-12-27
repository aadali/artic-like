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

#include<string>
#include<valarray>
using namespace std;
  
struct cpoint {
  double x;
  double y;
  double z;

  cpoint();
  cpoint(double, double, double);
  
  cpoint interpolateto(cpoint, double);
  cpoint partwayto(cpoint, double);
  static double saddlepoint(cpoint, cpoint, cpoint, cpoint);

  string metapoint();
};


class cparray {
  valarray<cpoint>* a;
  int i;
  int j;
public:
  cparray(int x, int y);        // note: no default constructor
  cparray(const cparray&);
  cparray& operator=(const cparray&);
  ~cparray();

  int size() const { return i*j; }
  int i_dim() const { return i; }
  int j_dim() const { return j; }

  cpoint& operator()(int x, int y);         // Fortran-style subscripts
  cpoint operator()(int x, int y) const;

  valarray<cpoint>& array() { return *a; }
};


