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

#include"metacontour.h"
using namespace std;

MetapostContours::MetapostContours(cparray& PlotData,
                                   vector<double>& ContourListIn,
                                   string& PlotNameIn)
{
  ContourList = ContourListIn;
  PlotName = PlotNameIn;

  GridCommand = "addto " + PlotName + ".Grid doublepath ";
  DrawCommand = "addto " + PlotName + ".LinePlot doublepath ";
  FillCommand = "addto " + PlotName + ".FillPlot contour ";

  GridStream << "picture " << PlotName << ".Grid; "
             << PlotName << ".Grid := nullpicture;\n";
  LinePlotStream << "picture " << PlotName << ".LinePlot; "
                 << PlotName << ".LinePlot := nullpicture;\n";
  FillPlotStream << "picture " << PlotName << ".FillPlot; "
                 << PlotName << ".FillPlot := nullpicture;\n";


  /* Construct contours */

  cpoint x00; cpoint x01; cpoint x10; cpoint x11;
  for (int j=1; j<PlotData.j_dim(); j++) {
    x00 = PlotData(0,j-1);
    x01 = PlotData(0,j);
    DrawGrid(x00, x01);
  }
  for (int i=1; i<PlotData.i_dim(); i++) {
    x00 = PlotData(i-1,0);
    x10 = PlotData(i,0);
    DrawGrid(x00, x10);
    for (int j=1; j<PlotData.j_dim(); j++) {
      x00 = PlotData(i-1, j-1);
      x01 = PlotData(i-1,  j );
      x10 = PlotData( i , j-1);
      x11 = PlotData( i ,  j );
      DrawGrid(x10, x11);
      DrawGrid(x01, x11);
      PlotSquare(x00, x10, x11, x01);
    }
  }
}


/* Start plotting a given grid square */
//
// Search Region: x3--------x2
//                |oooooooooo|
//                |oooooooooo|
//                |oooooooooo|
//                |oooooooooo|
//                x4--------x1
//
void MetapostContours::PlotSquare(cpoint x1, cpoint x2,
                             cpoint x3, cpoint x4)
{
  /* Organize things so x4 is smallest */

  if (x1.z < x4.z) {   // 1 is lowest of (1,4)
    if (x2.z < x1.z) {   // 2 is lowest of (1,2,3,4)
      if (x3.z < x2.z)  {   // 3 is lowest of all
        cpoint xtemp = x4; x4 = x3; x3 = x2; x2 = x1; x1 = xtemp;
      } else {             // 2 is lowest of all
        cpoint xtemp = x4; x4 = x2; x2 = xtemp; xtemp = x3; x3 = x1; x1 = xtemp;
      }
    } else {             // 1 is lowest of (1,2,4)
      if (x3.z < x1.z)  {   // 3 is lowest of all
        cpoint xtemp = x4; x4 = x3; x3 = x2; x2 = x1; x1 = xtemp;
      } else {             // 1 is lowest of all
        cpoint xtemp = x4; x4 = x1; x1 = x2; x2 = x3; x3 = xtemp;
      }
    }
  } else {             // 4 is lowest of (1,4)
    if (x2.z < x4.z) {   // 2 is lowest of (1,2,4)
      if (x3.z < x2.z)  {   // 3 is lowest of all
        cpoint xtemp = x4; x4 = x3; x3 = x2; x2 = x1; x1 = xtemp;
      } else {             // 2 is lowest of all
        cpoint xtemp = x4; x4 = x2; x2 = xtemp; xtemp = x3; x3 = x1; x1 = xtemp;
      }
    } else {             // 4 is lowest of (1,2,4)
      if (x3.z < x4.z)  {   // 3 is lowest of all
        cpoint xtemp = x4; x4 = x3; x3 = x2; x2 = x1; x1 = xtemp;
      }                    // 4 is lowest of all, do nothing
    }
  }  
    
  /* Find initial level, draw square, start rest if appropriate */
  if (ContourList[ContourList.size()-1] < x4.z) {
    DrawFilled(ContourList.size()-1, x1, x4, x3, x2);
  } else {
    int contourLevel = 0;
    while (ContourList[contourLevel] < x4.z &&
           contourLevel < ContourList.size()-1) contourLevel++;
    DrawFilled(contourLevel-1, x1, x4, x3, x2);
    PlotPent(x1, x2, x3, x4, contourLevel);
  }
}


/* Continue plotting for a corner of a grid square */
//
// Search Region: x3--------x2
//                |..........|
//                |........oo|
//                |......oooo|
//                |....oooooo|
//                x4--------x1
//
void MetapostContours::PlotCorner(cpoint x1, cpoint x2,
                             cpoint x3, cpoint x4, int contourLevel)
{
  if (contourLevel >= ContourList.size()) return;
  double zc2 = ContourList[contourLevel];
  if (x1.z > zc2) {
    /* another corner */
    cpoint y1 = x4.interpolateto(x1, zc2);
    cpoint y2 = x2.interpolateto(x1, zc2);
    DrawFilled(contourLevel, y1, y2, x1);
    DrawLine(contourLevel, y1, y2);
    PlotCorner(x1, x2, x3, x4, contourLevel+1);
  }
}


/* Continue plotting for a two-corner side of a grid square */
//
// Search Region: x3--------x2
//                |.....ooooo|
//                |.....ooooo|
//                |.....ooooo|
//                |.....ooooo|
//                x4--------x1
//
void MetapostContours::PlotRect(cpoint x1, cpoint x2,
                           cpoint x3, cpoint x4, int contourLevel)
{
  if (contourLevel >= ContourList.size()) return;
  double zc2 = ContourList[contourLevel];
  if (x1.z > zc2) {
    if (x2.z > zc2) {
      /* another rectangle */
      cpoint y1 = x4.interpolateto(x1, zc2);
      cpoint y2 = x3.interpolateto(x2, zc2);
      DrawFilled(contourLevel, y1, y2, x2, x1);
      DrawLine(contourLevel, y1, y2);
      PlotRect(x1, x2, x3, x4, contourLevel+1);
    } else {
      /* corner around x1 */
      cpoint y1 = x4.interpolateto(x1, zc2);
      cpoint y2 = x2.interpolateto(x1, zc2);
      DrawFilled(contourLevel, y1, y2, x1);
      DrawLine(contourLevel, y1, y2);
      PlotCorner(x1, x2, x3, x4, contourLevel+1);
    }
  } else {
    if (x2.z > zc2) {
      /* corner around x2 */
      cpoint y1 = x1.interpolateto(x2, zc2);
      cpoint y2 = x3.interpolateto(x2, zc2);
      DrawFilled(contourLevel, y1, y2, x2);
      DrawLine(contourLevel, y1, y2);
      PlotCorner(x2, x3, x4, x1, contourLevel+1);
    }
  }
}

  
/* Continue plotting for a three-corner pentagon of a grid square */
//
// Search Region: x3--------x2
//                |oooooooooo|
//                |..oooooooo|
//                |....oooooo|
//                |......oooo|
//                x4--------x1
//
void MetapostContours::PlotPent(cpoint x1, cpoint x2,
                           cpoint x3, cpoint x4, int contourLevel)
{
  if (contourLevel >= ContourList.size()) return;
  double zc2 = ContourList[contourLevel];
  if (x1.z > zc2) {
    if (x3.z > zc2) {
      if (x2.z < zc2) {
        /* a saddle-point hexagon, or two corners */
        cpoint y1 = x4.interpolateto(x1, zc2);
        cpoint y2 = x2.interpolateto(x1, zc2);
        cpoint y3 = x2.interpolateto(x3, zc2);
        cpoint y4 = x4.interpolateto(x3, zc2);
        double zcenter = cpoint::saddlepoint(x1, x2, x3, x4);
        if (zcenter < zc2) {
          /* two corners around x1 and x3 */
          DrawFilled(contourLevel, y1, y2, x1);
          DrawLine(contourLevel, y1, y2);
          PlotCorner(x1, x2, x3, x4, contourLevel+1);
          DrawFilled(contourLevel, y3, y4, x3);
          DrawLine(contourLevel, y3, y4);
          PlotCorner(x3, x4, x1, x2, contourLevel+1);
        } else {
          /* a saddle-point hexagon */
          DrawFilled(contourLevel, y1, y4, x3, y3, y2, x1);
          DrawLine(contourLevel, y1, y4);
          DrawLine(contourLevel, y3, y2);
          PlotHex(x1, x2, x3, x4, contourLevel+1, zcenter);
        }          
      } else {
        /* another pentagon */
        cpoint y1 = x4.interpolateto(x1, zc2);
        cpoint y2 = x4.interpolateto(x3, zc2);
        DrawFilled(contourLevel, y1, y2, x3, x2, x1);
        DrawLine(contourLevel, y1, y2);
        PlotPent(x1, x2, x3, x4, contourLevel+1);
      }
    } else {                // (x3.z < zc2)
      if (x2.z > zc2) {
        /* rectangle around x1, x2 */
        cpoint y1 = x4.interpolateto(x1, zc2);
        cpoint y2 = x3.interpolateto(x2, zc2);
        DrawFilled(contourLevel, y1, y2, x2, x1);
        DrawLine(contourLevel, y1, y2);
        PlotRect(x1, x2, x3, x4, contourLevel+1);
      } else {
        /* corner around x1 */
        cpoint y1 = x4.interpolateto(x1, zc2);
        cpoint y2 = x2.interpolateto(x1, zc2);
        DrawFilled(contourLevel, y1, y2, x1);
        DrawLine(contourLevel, y1, y2);
        PlotCorner(x1, x2, x3, x4, contourLevel+1);
      }
    }
  } else {                  // (x1.z < zc2)
    if (x3.z > zc2) {
      if (x2.z > zc2) {
        /* rectangle around x2, x3 */
        cpoint y1 = x1.interpolateto(x2, zc2);
        cpoint y2 = x4.interpolateto(x3, zc2);
        DrawFilled(contourLevel, y1, y2, x3, x2);
        DrawLine(contourLevel, y1, y2);
        PlotRect(x2, x3, x4, y1, contourLevel+1);
      } else {
        /* corner around x3 */
        cpoint y1 = x2.interpolateto(x3, zc2);
        cpoint y2 = x4.interpolateto(x3, zc2);
        DrawFilled(contourLevel, y1, y2, x3);
        DrawLine(contourLevel, y1, y2);
        PlotCorner(x3, x4, x1, x2, contourLevel+1);
      }
    } else {
      if (x2.z > zc2) {
        /* corner around x2 */
        cpoint y1 = x1.interpolateto(x2, zc2);
        cpoint y2 = x3.interpolateto(x2, zc2);
        DrawFilled(contourLevel, y1, y2, x2);
        DrawLine(contourLevel, y1, y2);
        PlotCorner(x2, x3, x4, x1, contourLevel+1);
      }
    }
  }
}


/* Continue plotting for a two-corner hex in a grid square */
//
// Search Region: x3--------x2
//                |oooooo....|
//                |oooooooo..|
//                |..oooooooo|
//                |....oooooo|
//                x4--------x1
//
void MetapostContours::PlotHex(cpoint x1, cpoint x2,
                          cpoint x3, cpoint x4,
                          int contourLevel, double zcenter)
{
  if (contourLevel >= ContourList.size()) return;
  double zc2 = ContourList[contourLevel];
  if (x1.z > zc2) {
    if (x3.z > zc2) {
      /* another saddle-point hexagon, or two corners */
      cpoint y1 = x4.interpolateto(x1, zc2);
      cpoint y2 = x2.interpolateto(x1, zc2);
      cpoint y3 = x2.interpolateto(x3, zc2);
      cpoint y4 = x4.interpolateto(x3, zc2);
      if (zcenter < zc2) {
        /* two corners around x1 and x3 */
        DrawFilled(contourLevel, y1, y2, x1);
        DrawLine(contourLevel, y1, y2);
        PlotCorner(x1, x2, x3, x4, contourLevel+1);
        DrawFilled(contourLevel, y3, y4, x3);
        DrawLine(contourLevel, y3, y4);
        PlotCorner(x3, x4, x1, x2, contourLevel+1);
      } else {
        /* another saddle-point hexagon */
        DrawFilled(contourLevel, y1, y4, x3, y3, y2, x1);
        DrawLine(contourLevel, y1, y4);
        DrawLine(contourLevel, y3, y2);
        PlotHex(x1, x2, x3, x4, contourLevel+1, zcenter);
      }          
    } else {
      /* corner around x1 */
      cpoint y1 = x4.interpolateto(x1, zc2);
      cpoint y2 = x2.interpolateto(x1, zc2);
      DrawFilled(contourLevel, y1, y2, x1);
      DrawLine(contourLevel, y1, y2);
      PlotCorner(x1, x2, x3, x4, contourLevel+1);
    }
  } else {
    if (x3.z > zc2) {
      /* corner around x3 */
      cpoint y1 = x2.interpolateto(x3, zc2);
      cpoint y2 = x4.interpolateto(x3, zc2);
      DrawFilled(contourLevel, y1, y2, x2);
      DrawLine(contourLevel, y1, y2);
      PlotCorner(x3, x4, x1, x2, contourLevel+1);
    }
  }
}
  
void MetapostContours::DrawGrid (cpoint x1, cpoint x2) {
  GridStream << GridCommand << x1.metapoint() << "--" << x2.metapoint()
             << ";\n";
}

void MetapostContours::DrawLine   (int contourLevel, cpoint x1, cpoint x2) {
  LinePlotStream << DrawCommand << x1.metapoint() << "--" << x2.metapoint()
                 << " withcolor contourcolor"
                 << contourLevel+1 << ";\n";
}

void MetapostContours::DrawFilled (int contourLevel, cpoint x1, cpoint x2, cpoint x3) {
  FillPlotStream << FillCommand << x1.metapoint() << "--" << x2.metapoint()
                 << "--" << x3.metapoint()
                 << "--cycle withcolor contourcolor" << contourLevel+1 << ";\n";
}

void MetapostContours::DrawFilled (int contourLevel, cpoint x1, cpoint x2, cpoint x3,
                                   cpoint x4) {
  FillPlotStream << FillCommand << x1.metapoint() << "--" << x2.metapoint()
                 << "--" << x3.metapoint() << "--" << x4.metapoint()
                 << "--cycle withcolor contourcolor" << contourLevel+1 << ";\n";
}

void MetapostContours::DrawFilled (int contourLevel, cpoint x1, cpoint x2, cpoint x3,
                                   cpoint x4, cpoint x5) {
  FillPlotStream << FillCommand << x1.metapoint() << "--" << x2.metapoint()
                 << "--" << x3.metapoint() << "--" << x4.metapoint()
                 << "--" << x5.metapoint()
                 << "--cycle withcolor contourcolor" << contourLevel+1 << ";\n";
}

void MetapostContours::DrawFilled (int contourLevel, cpoint x1, cpoint x2, cpoint x3,
                                   cpoint x4, cpoint x5, cpoint x6) {
  FillPlotStream << FillCommand << x1.metapoint() << "--" << x2.metapoint()
                 << "--" << x3.metapoint() << "--" << x4.metapoint()
                 << "--" << x5.metapoint() << "--" << x6.metapoint()
                 << "--cycle withcolor contourcolor" << contourLevel+1 << ";\n";
}


