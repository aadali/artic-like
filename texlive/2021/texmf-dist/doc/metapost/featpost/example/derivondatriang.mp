
verbatimtex
\documentclass{article}
\usepackage{palatino}
\usepackage{concmath}
\usepackage{ccfonts}
\begin{document}
etex

def marktic( expr Where, IsVertical, Size ) =
  begingroup
    pair low, hig;
    if IsVertical:
      low = ( xpart Where, (ypart Where) - 0.5*Size );
      hig = ( xpart Where, (ypart Where) + 0.5*Size );
    else:
      low = ( (xpart Where) - 0.5*Size, ypart Where );
      hig = ( (xpart Where) + 0.5*Size, ypart Where );
    fi;
    draw low--hig
  endgroup
enddef;
  
beginfig(1);
  numeric u, xsize, ysize, ticlen, xstep, yamp, yshift;
  path o.dtriang, xax, yax;
  picture previouspic;
  u = 8mm;
  xsize = 9u;
  ysize = 2.7u;
  ticlen = 0.2u;
  yamp = 1.8u;
  xstep = 1.6u;
  yshift = 6.0u;
  z0 = origin;
  z1 = ( xsize, 0 );
  xax = z0--z1;
  z2 = ( 0, -ysize );
  z3 = ( 0, ysize );
  yax = z2--z3;
  z4 = ( 0, -yamp );
  z40 = ( 0, yamp );
  z5 = ( xstep, yamp );
  z50 = ( xstep, 0 );
  z6 = ( 2*xstep, -yamp );
  z60 = ( 2*xstep, 0 );
  z7 = ( 3*xstep, yamp );
  z70 = ( 3*xstep, 0 );
  z8 = ( 4*xstep, -yamp );
  z80 = ( 4*xstep, 0 );
  z9 = ( 5*xstep, yamp );
  z90 = ( 5*xstep, 0 );
  o.dtriang = z4--z5--z6--z7--z8--z9;
  drawarrow xax;
  drawarrow yax;
  marktic( z4, false, ticlen );
  marktic( z40, false, ticlen );
  marktic( z50, true, ticlen );
  marktic( z60, true, ticlen );
  marktic( z70, true, ticlen );
  marktic( z80, true, ticlen );
  marktic( z90, true, ticlen );
  label.bot( btex $x$ etex, z1 );
  label.bot( btex 1 etex, z50 );
  label.bot( btex 2 etex, z60 );
  label.bot( btex 3 etex, z70 );
  label.bot( btex 4 etex, z80 );
  label.bot( btex 5 etex, z90 );
  previouspic = currentpicture;
  draw previouspic shifted (0,-yshift);
  label.lft( btex $\frac{d\, f}{d\, x}$ etex scaled 1.35, (0,ysize-yshift) );
  label.lft( btex $f$ etex, z3 );
  draw o.dtriang withpen pencircle scaled 2pt;
  label.lft( btex 2 etex, z40 );
  label.lft( btex -2 etex, z4 );
endfig;

verbatimtex
\end{document}
etex

end.
