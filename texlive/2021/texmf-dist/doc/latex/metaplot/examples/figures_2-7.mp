% Figures 2-7 from _MetaPlot, MetaContour, and Other Collaborations
% with METAPOST_, preprint version, given at the Practical TeX 2004
% conference.

prologues:=0;
input metaplot    % MetaPlot macros

% capillary.mp and capillary2.mp are both produced as output of the 
% program capillary.cc.  The relevant parameters in the program
% (lines 29-30) are:
%
% For capillary.mp:
%    double theta = pi/4.0;
%    double d = 1.0;
%
% For capillary2.mp:
%    double theta = pi/6.0;
%    double d = 1.0;

input capillary   % capillary plot object
input capillary2  % capillaryb plot object


% Figure 2: The capillary surface, in its unadorned form as plot
% object elements scaled to 2.0in by 0.75in.

plot_instantiate(plotA, capillary)
plotA.pageleft = 0.0;
plotA.pagebottom = 0.0;
plotA.pagewidth = 2.0in;
plotA.pageheight = 0.75in;

beginfig(2)
  draw plotA.plot(fplot) withcolor 0.85white;
  draw plotA.plot(lplot) withpen pencircle scaled 1pt;
endfig;


% Figure 3: The capillary surface, with equal y and z scales, a page
% height of 0.75in, and appropriate annotations.

plot_instantiate(plotB, capillary)
plotB.pageleft = 0.0;
plotB.pagebottom = 0.0;
plotB.pageheight = 0.75in;
plot_setequalaxes(plotB);

beginfig(3)
  draw plotB.plot(fplot) withcolor 0.85white;
  linecap := butt;
  pickup pencircle scaled 1pt;
  % z-axis (vertical)
  z1 = (plotB.pageleft, plotB.pagebottom);
  z2 = (plotB.pageleft, plotB.pagetop + 0.1in);
  % y-axis (horizontal)
  z3 = (plotB.pageleft, plot_ypageloc(plotB,0.0));
  z4 = (plotB.pageright + 0.1in,
          plot_ypageloc(plotB,0.0));
  drawarrow z1 -- z2;
  label.top(btex $z$ etex, z2);
  drawarrow z3 -- z4;
  label.rt(btex $y$ etex, z4);
  pickup pencircle scaled 0.25pt;
  % Label for contact angle
  z5 = plotB.plot(contactpoint);
  z6 = z5 + 0.24in
        * dir(-90 + capillary.contactangle);
  z7 = z5 + 0.18in
        * dir(-90 + 0.5*capillary.contactangle);
  draw z5 -- z6;
  label(btex $\theta$ etex, z7);
endfig;


% Figure 4: The capillary surface with parameters and page height as
% in Figure 3, but with $\theta = \pi/6$.

plot_instantiate(plotC, capillaryb)
plotC.pageleft = 0.0;
plotC.pagebottom = 0.0;
plotC.pageheight = 0.75in;
plot_setequalaxes(plotC);

beginfig(4)
  draw plotC.plot(fplot) withcolor 0.85white;
  linecap := butt;
  pickup pencircle scaled 1pt;
  % z-axis (vertical)
  z1 = (plotC.pageleft, plotC.pagebottom);
  z2 = (plotC.pageleft, plotC.pagetop + 0.1in);
  % y-axis (horizontal)
  z3 = (plotC.pageleft, plot_ypageloc(plotC,0.0));
  z4 = (plotC.pageright + 0.1in,
          plot_ypageloc(plotC,0.0));
  drawarrow z1 -- z2;
  label.top(btex $z$ etex, z2);
  drawarrow z3 -- z4;
  label.rt(btex $y$ etex, z4);
  pickup pencircle scaled 0.25pt;
  % Label for contact angle
  z5 = plotC.plot(contactpoint);
  z6 = z5 + 0.24in
        * dir(-90 + capillaryb.contactangle);
  z7 = z5 + 0.18in
        * dir(-90 + 0.5*capillaryb.contactangle);
  draw z5 -- z6;
  label(btex $\theta$ etex, z7);
endfig;


% Figure 5: Two capillary surfaces, as in Figure 3 and Figure 4, showing
% the difference in the curves as a result of varying $\theta$.

plot_instantiate(plotD, capillary)
plotD.pageleft = 0.0;
plotD.pagebottom = 0.0;
plotD.pageheight = 0.75in;
plot_setequalaxes(plotD);

plot_instantiate(plotE, capillaryb)
plot_zpageloc(plotD, 0.0, 0.0)
 = plot_zpageloc(plotE, 0.0, 0.0);
plot_zpageloc(plotD, 1.0, 1.0)
 = plot_zpageloc(plotE, 1.0, 1.0);

beginfig(5)
  linecap := butt;
  pickup pencircle scaled 1pt;
  draw plotD.plot(fplot) withcolor 0.85white;
  draw plotE.plot(lplot) dashed evenly
    withpen pencircle scaled 0.5pt;
  % z-axis (vertical)
  z1 = (plotD.pageleft, plotD.pagebottom);
  x2 = plotD.pageleft;
  y2 = max(plotD.pagetop, plotE.pagetop) + 0.1in;
  % y-axis (horizontal)
  z3 = (plotD.pageleft, plot_ypageloc(plotD,0.0));
  x4 = max(plotD.pageright, plotE.pageright) + 0.1in;
  y4 = plot_ypageloc(plotD,0.0);
  drawarrow z1 -- z2;
  label.top(btex $z$ etex, z2);
  drawarrow z3 -- z4;
  label.rt(btex $y$ etex, z4);
endfig;
end


% Figure 6: A repeat of Figure 5, with simple grid labels added.

beginfig(6)
  linecap := butt;
  pickup pencircle scaled 1pt;
  draw plotD.plot(fplot) withcolor 0.85white;
  draw plotE.plot(lplot) dashed evenly
    withpen pencircle scaled 0.5pt;
  % z-axis (vertical)
  z1 = (plotD.pageleft, plotD.pagebottom);
  x2 = plotD.pageleft;
  y2 = max(plotD.pagetop, plotE.pagetop) + 0.1in;
  % y-axis (horizontal)
  z3 = (plotD.pageleft, plot_ypageloc(plotD,0.0));
  x4 = max(plotD.pageright, plotE.pageright) + 0.1in;
  y4 = plot_ypageloc(plotD,0.0);
  drawarrow z1 -- z2;
  label.top(btex $z$ etex, z2);
  drawarrow z3 -- z4;
  label.rt(btex $y$ etex, z4);
  x5 = plotD.pageleft;
  x6 = x4;
  y5 = y6 = plotD.pagebottom;
  draw plot_xtickscale(plotD)(
    z5, z6, 0.08in, 0.06in, down,
    0.0, 1.0, "%3f")
    withpen pencircle scaled 0.5pt;
  y7 = plotD.pagebottom;
  y8 = y2;
  x7 = x8 = plotD.pageleft;
  draw plot_ytickscale(plotD)(
    z7, z8, 0.08in, 0.06in, left,
    0.0, 0.5, "%3f")
    withpen pencircle scaled 0.5pt;
endfig;


% Figure 7: A repeat of Figure 5 again, with more advanced grid labels.

beginfig(7)
  linecap := butt;
  pickup pencircle scaled 1pt;
  draw plotD.plot(fplot) withcolor 0.85white;
  draw plotE.plot(lplot) dashed evenly
    withpen pencircle scaled 0.5pt;
  % z-axis (vertical)
  z1 = (plotD.pageleft, plotD.pagebottom);
  x2 = plotD.pageleft;
  y2 = max(plotD.pagetop, plotE.pagetop) + 0.1in;
  % y-axis (horizontal)
  z3 = (plotD.pageleft, plot_ypageloc(plotD,0.0));
  x4 = max(plotD.pageright, plotE.pageright) + 0.1in;
  y4 = plot_ypageloc(plotD,0.0);
  drawarrow z1 -- z2;
  label.top(btex $z$ etex, z2);
  drawarrow z3 -- z4;
  label.rt(btex $y$ etex, z4);
  x5 = plotD.pageleft;
  x6 = x4 - 0.1in;
  y5 = y6 = plotD.pagebottom - 0.06in;
  draw plot_xtickscale(plotD)(
    z5, z6, 0.08in, 0.06in, down,
    0.5, 1.0, "%3f")
    withpen pencircle scaled 0.5pt;
  draw plot_xtickscale(plotD)(
    z5, z6, 0.08in, 0.06in, down,
    0.0, 1.0, "")
    withpen pencircle scaled 0.5pt;
  draw plot_xtickscale(plotD)(
    z5, z6, 0.04in, 0.06in, down,
    0.0, 0.1, "")
    withpen pencircle scaled 0.5pt;
  y7 = y4;
  y8 = y2 - 0.1in;
  x7 = x8 = plotD.pageleft - 0.06in;
  draw plot_ytickscale(plotD)(
    z7, z8, 0.08in, 0.06in, left,
    0.0, 0.5, "%3f")
    withpen pencircle scaled 0.5pt;
  draw plot_ytickscale(plotD)(
    z7, z8, 0.04in, 0.06in, left,
    0.0, 0.1, "")
    withpen pencircle scaled 0.5pt;
endfig;

end
