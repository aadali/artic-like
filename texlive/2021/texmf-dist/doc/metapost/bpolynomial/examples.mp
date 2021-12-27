input bpolynomial;
input metafun
input graph
input hatching
prologues := 3;

transform T;
T := identity xscaled 10mm yscaled 1mm;

beginfig(1);
  newBPolynomial.f(2, 0, -3, -1);
  draw f.getPath(-2, 2) transformed T;
endfig;


beginfig(2);
  newBPolynomial.f(2, 0, -3, -1);
  draw f.getPath(-2, 2) transformed T;
  draw f.getTangent(-0.5)(-1, 1) transformed T;
endfig;


beginfig(3);
  dotlabeldiam := 2bp;
  labeloffset := 10bp;
  newBPolynomial.f(2, 0, -3, -1);
  draw f.getPath(-2, 2) transformed T;
  x := -0.5;
  show (x, f.eval(x));
  draw f.getTangent(x)(-1, 1) transformed T;
  dotlabel.top(btex $(-0.5, 0.25)$ etex, (x, f.eval(x)) transformed T);
endfig;

labeloffset := 3bp;

beginfig(4);
  newBPolynomial.f(2, 0, -3, -1);
  draw f'.getPath(-2, 2) transformed T;
  draw f'.getTangent(-0.25)(-1, 1) transformed T;
endfig;


beginfig(5);
  T := identity scaled 10mm;
  newBSqrRoot.s(1,0,0);
  newBCubRoot.c(1,0,0);
  draw s.getPath(0,6) transformed T;
  draw c.getPath(0,6) transformed T;
  draw s.getTangent(3)(-2, 2) transformed T;
endfig;


T := identity xscaled 10mm yscaled 1mm;
beginfig(6);
  draw getBezierFromPolynomial(2, 0, -3, -1)(-2, 2) transformed T;
endfig;


beginfig(11);
path f, g;
  xmin := -7; xmax := 7;
  ymin := -7; ymax := 7;
  %%% Define polynomials f and g.
  newBPolynomial.f(0.3, 0, -3, -1);
  f := f.getPath(xmin, xmax);
  newBPolynomial.g(0, 0.5, -2, 0);
  g := g.getPath(xmin, xmax);
  %%% Draw graph.
  draw begingraph(10cm, 6cm);
    setrange(xmin,ymin, xmax,ymax);
    autogrid(grid.bot, grid.lft) dashed evenly withcolor .9white;
    drawoptions(withpen pencircle scaled 1bp);
    gdraw f;
    gdraw g dashed evenly scaled 2;
    drawoptions();
  endgraph;
  show f;
  %%% Write table with some points of f to log file.
  show "Polynomial: " & decimal f.a & "x^3 + " &
    decimal f.b & "x^2 + " & decimal f.c & "x + " & decimal f.d;
  for x=-5 upto 5:
    show (x, f.eval(x));
  endfor
endfig;


beginfig(12);
  xmin := -6; xmax := 6;
  ymin := -6; ymax := 6;
  newBPolynomial.f(0.3, -0.5, -0.5, -1);
  draw begingraph(10cm, 6cm);
    setrange(xmin,ymin, xmax,ymax);
    autogrid(grid.bot, grid.lft) dashed evenly withcolor .9white;
    drawoptions(withpen pencircle scaled 1bp);
    %%% Draw f and its derivatives f', f'', f'''.
    gdraw f.getPath(xmin, xmax);
    gdraw f'.getPath(xmin, xmax) dashed evenly scaled 2;
    gdraw f''.getPath(xmin, xmax) dashed withdots
      withpen pencircle scaled 2bp;
    gdraw f'''.getPath(-5, 5) withcolor .6white;
    %%% Draw tangents and mark points.
    x := 2;
    drawoptions(withcolor (1, 0.6, 0.6));
    gdraw f.getTangent(x)(-2, 2);
    gdraw f'.getTangent(x)(-1, 1);
    gdraw f''.getTangent(x)(-2, 2);
    gdraw f'''.getTangent(x)(-2, 2);
    drawoptions(withcolor (0.6, 0.6, 1));
    dotlabeldiam := 2.5bp;
    gdotlabel("", (x, f.eval(x)));
    gdotlabel("", (x, f'.eval(x)));
    gdotlabel("", (x, f''.eval(x)));
    gdotlabel("", (x, f'''.eval(x)));
    drawoptions();
  endgraph;
endfig;


beginfig(13);
path f, g, A;
  xmin := -3; xmax := 6;
  ymin := -3; ymax := 6;
  newBPolynomial.f(-0.25, 0.5, 2, -1);
  newBPolynomial.g(0, 0.5, -2, 0);
  f := f.getPath(-2.5, 5.5);
  g := g.getPath(-1.5, 5.5);
  %%% Find area between f and g.
  A := buildcycle(g, reverse f);
  draw begingraph(10cm, 6cm);
    setrange(xmin,ymin, xmax,ymax);
    autogrid(grid.bot, grid.lft) dashed evenly withcolor .9white;
    %%% Fill area with transparent colour.
    gfill A withcolor transparent (1, .3, (1, 0.5, 0));
    drawoptions(withpen pencircle scaled 1bp);
    gdraw f;
    gdraw g dashed evenly scaled 2;
    drawoptions();
  endgraph;
endfig;


beginfig(14);
path f, g, A;
  T := identity xscaled 10mm yscaled 6mm;
  %%% Draw coordinate system.
  xmin := -3; xmax := 6;
  ymin := -3; ymax := 6;
  drawoptions(withpen pencircle scaled 1bp withcolor 0.8white);
  drawarrow ((xmin,0)--(xmax,0)) transformed T;
  drawarrow ((0,ymin)--(0,ymax)) transformed T;
  newBPolynomial.f(-0.25, 0.5, 2, -1);
  newBPolynomial.g(0, 0.5, -2, 0);
  f := f.getPath(-2.5, 4.2);
  g := g.getPath(-1, 5);
  A := buildcycle(g, reverse f);
  %%% Fill area with pattern.
  drawoptions();
  hatchoptions(withcolor (0.6, 0.3, 0.3));
  hatchfill A transformed T
    withcolor (-45, 2mm, -0.5bp) withcolor (45, 2mm, -0.5bp);
  drawoptions(withpen pencircle scaled 1bp);
  draw f transformed T;
  draw g transformed T dashed evenly scaled 2;
endfig;


%%%  The following figures work around a bug in metafun's
%%%  mp-form.mp package for the original figures 11 to 13.
%%%  The bug shows up when rendering negative numbers
%%%  on corrdinate axes using macro 'format'.
beginfig(21);
path f, g;
  xmin := -7; xmax := 7;
  ymin := -7; ymax := 7;
  %%% Define polynomials f and g.
  newBPolynomial.f(0.3, 0, -3, -1);
  f := f.getPath(xmin, xmax);
  newBPolynomial.g(0, 0.5, -2, 0);
  g := g.getPath(xmin, xmax);
  %%% Draw graph.
  draw begingraph(10cm, 6cm);
    setrange(xmin,ymin, xmax,ymax);
%     autogrid(grid.bot, grid.lft) dashed evenly withcolor .9white;
    for i=xmin+1 step 2 until xmax-1:
      grid.bot(decimal i, i) dashed evenly withcolor .9white;
    endfor
    for i=ymin+1 step 2 until ymax-1:
      grid.lft(decimal i, i) dashed evenly withcolor .9white;
    endfor
    drawoptions(withpen pencircle scaled 1bp);
    gdraw f;
    gdraw g dashed evenly scaled 2;
    drawoptions();
  endgraph;
  show f;
  %%% Write table with some points of f to log file.
  show "Polynomial: " & decimal f.a & "x^3 + " &
    decimal f.b & "x^2 + " & decimal f.c & "x + " & decimal f.d;
  for x=-5 upto 5:
    show (x, f.eval(x));
  endfor
endfig;


beginfig(22);
  xmin := -6; xmax := 6;
  ymin := -6; ymax := 6;
  newBPolynomial.f(0.3, -0.5, -0.5, -1);
  draw begingraph(10cm, 6cm);
    setrange(xmin,ymin, xmax,ymax);
%     autogrid(grid.bot, grid.lft) dashed evenly withcolor .9white;
    for i=xmin step 2 until xmax:
      grid.bot(decimal i, i) dashed evenly withcolor .9white;
    endfor
    for i=ymin step 2 until ymax:
      grid.lft(decimal i, i) dashed evenly withcolor .9white;
    endfor
    drawoptions(withpen pencircle scaled 1bp);
    %%% Draw f and its derivatives f', f'', f'''.
    gdraw f.getPath(xmin, xmax);
    gdraw f'.getPath(xmin, xmax) dashed evenly scaled 2;
    gdraw f''.getPath(xmin, xmax) dashed withdots
      withpen pencircle scaled 2bp;
    gdraw f'''.getPath(-5, 5) withcolor .6white;
    %%% Draw tangents and mark points.
    x := 2;
    drawoptions(withcolor (1, 0.6, 0.6));
    gdraw f.getTangent(x)(-2, 2);
    gdraw f'.getTangent(x)(-1, 1);
    gdraw f''.getTangent(x)(-2, 2);
    gdraw f'''.getTangent(x)(-2, 2);
    drawoptions(withcolor (0.6, 0.6, 1));
    dotlabeldiam := 2.5bp;
    gdotlabel("", (x, f.eval(x)));
    gdotlabel("", (x, f'.eval(x)));
    gdotlabel("", (x, f''.eval(x)));
    gdotlabel("", (x, f'''.eval(x)));
    drawoptions();
  endgraph;
endfig;


beginfig(23);
path f, g, A;
  xmin := -3; xmax := 6;
  ymin := -3; ymax := 6;
  newBPolynomial.f(-0.25, 0.5, 2, -1);
  newBPolynomial.g(0, 0.5, -2, 0);
  f := f.getPath(-2.5, 5.5);
  g := g.getPath(-1.5, 5.5);
  %%% Find area between f and g.
  A := buildcycle(g, reverse f);
  draw begingraph(10cm, 6cm);
    setrange(xmin,ymin, xmax,ymax);
%     autogrid(grid.bot, grid.lft) dashed evenly withcolor .9white;
    for i=xmin+1 step 2 until xmax:
      grid.bot(decimal i, i) dashed evenly withcolor .9white;
    endfor
    for i=ymin+1 step 2 until ymax:
      grid.lft(decimal i, i) dashed evenly withcolor .9white;
    endfor
    %%% Fill area with transparent colour.
    gfill A withcolor transparent (1, .3, (1, 0.5, 0));
    drawoptions(withpen pencircle scaled 1bp);
    gdraw f;
    gdraw g dashed evenly scaled 2;
    drawoptions();
  endgraph;
endfig;

end
