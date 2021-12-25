% Figure 1 from _MetaPlot, MetaContour, and Other Collaborations
% with METAPOST_, preprint version, given at the Practical TeX 2004
% conference.

prologues :=0;

% cap1.mp is the output from the cap1.cc program.
input cap1

% Define pens, for convenience.
pen thickline; thickline := pencircle scaled 2pt;
pen thinline; thinline := pencircle scaled 1pt;
pen hairline; hairline := pencircle scaled 0.25pt;

beginfig(1)
  pair contactpoint; contactpoint := ulcorner capillary scaled 0.5in;
  pair rightend; rightend := (xpart(lrcorner capillary) * 0.5in, 0);
  linecap := butt;
  draw (capillary scaled 0.5in) withcolor (0.85,0.85,0.85);
  drawarrow (0,-0.25in) -- (0,0.5in) withpen thinline;
  label.top(btex $z$ etex,(0,0.5in));
  drawarrow (0,0) -- (rightend + (0.1in,0)) withpen thinline;
  label.rt(btex $y$ etex, (rightend + (0.1in,0)));
  draw ((0,0) -- (0.24in, -0.24in)) shifted contactpoint withpen hairline;
  label(btex $\theta$ etex, (contactpoint + (0.07in, -0.18in)));
endfig;

end
  