% Figures 8-9 from _MetaPlot, MetaContour, and Other Collaborations
% with METAPOST_, preprint version, given at the Practical TeX 2004
% conference.

prologues:=0;
input metaplot    % MetaPlot macros


% mc2.mp is produced by the very pre-release version of metacontour
% that should be included with this file; the program consists of
% the files metacontour.cc, metacontour.h, cpoint.cc, cpoint.h,
% and metacontour_main.cc.  There is no input file; the data to
% be plotted is generated within the program.
%
% Note that the resulting file is not in true "plot object" form;
% we need to manually clip it and define the plot scale here.
% This will, of course, be fixed in later versions of metacontour.

input mc2
clip contplotA.LinePlot to (0,0)--(0,1)--(1,1)--(1,0)--cycle;
clip contplotA.FillPlot to (0,0)--(0,1)--(1,1)--(1,0)--cycle;
clip contplotA.Grid to (0,0)--(0,1)--(1,1)--(1,0)--cycle;

contplotA.xleft = -1.0;
contplotA.xright = 1.0;
contplotA.ybot = -1.0;
contplotA.ytop = 1.0;


% Define pens, for convenience.

pen thickline; thickline := pencircle scaled 2pt;
pen thinline; thinline := pencircle scaled 1pt;
pen hairline; hairline := pencircle scaled 0.25pt;


% Figure 8: Sample graph created by MetaContour and MetaPlot,
% showing potential lines for a combination of a linear gradient
% and a point source, plotted on a polar grid.

plot_instantiate(cplotA, contplotA);

cplotA.pagewidth = 2.25in;
plot_setequalaxes(cplotA);
cplotA.llft = (0,0);

beginfig(8);
  draw cplotA.plot(LinePlot) withpen thinline;
  draw cplotA.plot(Grid) withpen hairline;
  linejoin:=mitered;
  draw cplotA.llft--cplotA.lrt--cplotA.urt--cplotA.ulft--cycle withpen thickline;
  draw plot_xtickscale(cplotA)(
    cplotA.llft, cplotA.lrt,
    0.08in, 0.06in, down, 0.0, 0.5, "%3f")
    withpen thinline;
  draw plot_xtickscale(cplotA)(
    cplotA.llft, cplotA.lrt,
    0.05in, 0.06in, down, 0.0, 0.1, "")
    withpen thinline;
  draw plot_ytickscale(cplotA)(
    cplotA.llft, cplotA.ulft,
    0.08in, 0.06in, left, 0.0, 0.5, "%3f")
    withpen thinline;
  draw plot_ytickscale(cplotA)(
    cplotA.llft, cplotA.ulft,
    0.05in, 0.06in, left, 0.0, 0.1, "")
    withpen thinline;
endfig;  
  

% Figure 9: Another sample graph created by MetaContour and MetaPlot,
% illustrating a filled contour-plot style rather than using contour
% lines.

beginfig(9);
  draw cplotA.plot(FillPlot);
  draw cplotA.plot(Grid) withpen hairline;
  linejoin:=mitered;
  draw cplotA.llft--cplotA.lrt--cplotA.urt--cplotA.ulft--cycle withpen thickline;
  draw plot_xtickscale(cplotA)(
    cplotA.llft, cplotA.lrt,
    0.08in, 0.06in, down, 0.0, 0.5, "%3f")
    withpen thinline;
  draw plot_xtickscale(cplotA)(
    cplotA.llft, cplotA.lrt,
    0.05in, 0.06in, down, 0.0, 0.1, "")
    withpen thinline;
  draw plot_ytickscale(cplotA)(
    cplotA.llft, cplotA.ulft,
    0.08in, 0.06in, left, 0.0, 0.5, "%3f")
    withpen thinline;
  draw plot_ytickscale(cplotA)(
    cplotA.llft, cplotA.ulft,
    0.05in, 0.06in, left, 0.0, 0.1, "")
    withpen thinline;
endfig;

end
