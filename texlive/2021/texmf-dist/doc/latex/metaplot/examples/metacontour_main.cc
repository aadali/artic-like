/* MetaContour, version 0.1 -------------------------------------- */
/* Copyright(C) 2004, Brooks Moses                                 */
/*                                                                 */
/* This version of MetaContour is made available under the Gnu     */
/* Public License; see below for details.                          */
/*                                                                 */
/* This is a very pre-release version of MetaContour, distributed  */
/* primarily as an example of a use of MetaPlot.  It can be        */
/* compiled with gcc using the line:                               */
/*                                                                 */
/*   g++ metacontour.cc cpoint.cc metacontour_main.cc              */
/*                                                                 */
/*-----------------------------------------------------------------*/

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


#include<string>
#include<vector>
#include<iostream>
#include<cmath>
#include"metacontour.h"
using namespace std;

int main() {
  int il = 21;
  int jl = 11;
  double pi = atan(1.0)*4.0;
  
  cparray GridData(il,jl);

  for (int i=0; i<il; i++) {
    double theta = 2.0 * pi * double(i) / (il-1.0);
    for (int j=0; j<jl; j++) {
      double r = 1.5 * (double(j)/(jl-1.0))*(double(j)/(jl-1.0)) + 0.03;
      
      double x = - r * cos(theta);
      double y = r * sin(theta);
      
      double z = (-0.5 * log(r)) + x;
      
      GridData(i,j) = cpoint(0.5*(x+1.0),0.5*(y+1.0),z);
    }
  }

  vector<double> ContourList(31);

  cout << "color contourcolor[];\n";
  for (int k=0; k<31; k++) {
    ContourList[k] = 0.1 * k - 1.0;
    cout << "contourcolor[" << k << "] = " << (1.0 - k/50.0) << "white;\n";
  }

  string PlotName = "contplotA";
  
  MetapostContours testplot(GridData, ContourList, PlotName);
  
  cout << '\n' << testplot.Grid();
  cout << testplot.LinePlot();
  cout << testplot.FillPlot();
  
}

