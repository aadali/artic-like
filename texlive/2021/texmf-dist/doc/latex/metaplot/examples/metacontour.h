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

/* Class containing Metapost contour plot streams */

#include<string>
#include<sstream>
#include<vector>
#include"cpoint.h"

class MetapostContours {

  ostringstream GridStream;
  ostringstream FillPlotStream;
  ostringstream LinePlotStream;

  vector<double> ContourList;

  void PlotSquare (cpoint, cpoint, cpoint, cpoint);
  void PlotCorner (cpoint, cpoint, cpoint, cpoint, int);
  void PlotRect   (cpoint, cpoint, cpoint, cpoint, int);
  void PlotPent   (cpoint, cpoint, cpoint, cpoint, int);
  void PlotHex    (cpoint, cpoint, cpoint, cpoint, int, double);

  void DrawFilled (int, cpoint, cpoint, cpoint);
  void DrawFilled (int, cpoint, cpoint, cpoint, cpoint);
  void DrawFilled (int, cpoint, cpoint, cpoint, cpoint, cpoint);
  void DrawFilled (int, cpoint, cpoint, cpoint, cpoint, cpoint, cpoint);
  void DrawLine   (int, cpoint, cpoint);
  void DrawGrid   (cpoint, cpoint);
  
  string PlotName;
  string GridCommand;
  string DrawCommand;
  string FillCommand;
  
public:
  MetapostContours(cparray&, vector<double>&, string&);

  string Grid() const { return GridStream.str(); }
  string FillPlot() const { return FillPlotStream.str(); }
  string LinePlot() const { return LinePlotStream.str(); }

  vector<double> ContourListOut() const { return ContourList; }
  
};
