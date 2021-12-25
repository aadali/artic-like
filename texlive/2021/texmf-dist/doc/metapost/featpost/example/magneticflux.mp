% magneticflux.mp
% L. Nobre G. 
% IYP (2005)

input featpost3Dplus2D;

verbatimtex
\documentclass{article}
\usepackage{palatino,mathpazo}
\begin{document}
etex

TDAtiplen := 0.1;              % tdarrow
TDAhalftipbase := 0.04;         % Three-Dimensional
TDAhalfthick := 0.02;           % Arrow
ForePen := pencircle scaled 10pt;
BackPen := pencircle scaled 6pt;

numeric ranga;
ranga = 30+uniformdeviate(120);

def ellip( expr TheVal ) =
    ( (cosd(360*TheVal),sind(360*TheVal)/sind(ranga),0) )
enddef;

beginfig(2); 
  f := ( 4, 1, 2 );
  Spread := 90;
  pen blackpen, whitepen;
  numeric majoraxe, i, j, gridstepn, gridstep, gridstart, halen, rang;
  numeric ellfocus, totdist, marginl;
  path espira, partes, covere;
  color avec, midvec;
  halen = 1;
  marginl = 0.14;
  rang = 45+30*floor(uniformdeviate(3.9));
  show rang;
  avec = (0,cosd(rang),sind(rang))*halen;
  blackpen = pencircle scaled 5pt;
  whitepen = pencircle scaled 3pt;
  majoraxe = 1.2 + uniformdeviate(0.8);
  ellfocus = majoraxe +-+ 1;
  gridstepn= 15;
  gridstep = 2*majoraxe/(gridstepn-1);
  gridstart= 0.5*gridstepn*gridstep;
  espira = ellipticpath( black, -red, -green*majoraxe );
  draw espira withpen blackpen;
  undraw espira withpen whitepen;
  for j=-1 step gridstep until 1:
    for i=-gridstart step gridstep until gridstart:
      midvec := (j,i,0);
      totdist := conorm(midvec-(0,-ellfocus,0))+conorm(midvec-(0,ellfocus,0));
      if 2*majoraxe-marginl > totdist: 
	tdarrow( midvec, midvec+avec );
%	draw rp(midvec) withpen whitepen;
      fi;
    endfor;
  endfor;
  partes = subpath (10.3,25.6) of espira;
  covere = subpath (10,26) of espira;
  draw partes withpen blackpen;
  undraw covere withpen whitepen;
endfig;

beginfig(1);
  f := ( 4, 1.9, 1.22 );
  Spread := 90;
  numeric nodc, lsc;
  nodc = 10;
  lsc = 2.4;
  closedline( true, 200, 0.5, 1.75 )( ellip );
  numeric i, j, circlesize, polar, phi;
  color currpos, ini, fin, avec;
  string comm;
  avec = (0,cosd(ranga),sind(ranga));
  circlesize = 1.86/nodc;
  for i=0.93 step -circlesize until 0.4*circlesize:
    polar := floor(6.28318*i/circlesize);
    for j=1 upto polar:
      phi := 360*j/polar;
      currpos := i*(cosd(phi),sind(phi)/sind(ranga),0);
      comm:="tdarrow("
	& cstr( currpos )
	& ","
	& cstr( currpos + avec )
	& ");";
      getready( comm, currpos );
    endfor;
  endfor;
  doitnow;
  pickup pencircle scaled 2pt;
  drawarrow rp((0,1.07/sind(ranga),0))--rp((0,2.47,0));
  drawarrow rp((1.07,0,0))--rp((1.67,0,0));
  angline((0,2.37,0),(0,0.93/sind(ranga),0)+avec,(0,0.93/sind(ranga),0),
    0.3, btex $\alpha$ etex scaled lsc, urt);
  label.bot( btex $x$ etex scaled lsc, rp((1.67,0,0)) );
  label.bot( btex $y$ etex scaled lsc, rp((0,2.47,0)) );  
endfig;
  
  
verbatimtex
\end{document}
etex

end.
