% quartertorus.mp
% L. Nobre G.
% 2012

input featpost3Dplus2D;

prologues := 3;

verbatimtex
\documentclass{article} 
\usepackage{mathpazo} 
\begin{document} 
etex

f := (3,5,2.5);

beginfig(1);
  color a, b, c, g, h, i, j, l, m, n;
  numeric fac, siz, s;
  fac = 0.05;
  siz = 0.3;
  s = 2;
  a = blue-red;
  b = g-fac*a;
  c = -a-fac*g;
  g = -green;
  h = 2*g;
  i = red+g;
  j= -red;
  l = blue;
  m = 0.6red-g;
  n = -0.7*blue;
  drawoptions( withpen pencircle scaled 2.5pt );
  quartertorus( black, h, i, siz );
  quartertorus( black, g, a, siz );
  quartertorus( black, b, c, siz );
  quartertorus( black, -g, j, siz );
  quartertorus( n, l, m, siz );
  label( btex $a$ etex scaled s, rp( a ) );
  label.llft( btex $b$ etex scaled s, rp( b ) );
  label( btex $c$ etex scaled s, rp( c ) );
  label.urt( btex $g$ etex scaled s, rp( g ) );
  label( btex $-g$ etex scaled s, rp( -g ) );
  label( btex $h$ etex scaled s, rp( h ) );
  %label( btex $i$ etex scaled s, rp( i ) );
  label( btex $j$ etex scaled s, rp( j ) );
  label( btex $l$ etex scaled s, rp( l+n ) );
  label( btex $m$ etex scaled s, rp( m+n ) );
endfig;

beginfig(2);
  color col[];
  numeric ray, ang, len, bigray;
  ray = 0.15;
  bigray = 0.4;
  ang = 60;
  len = 0.7;
  col1 = (-(ray+bigray)/cosd(ang),0,0);
  col2 = (bigray*cosd(ang),0,bigray*sind(ang));
  col3 = (bigray*cosd(ang),0,-bigray*sind(ang));
  col4 = ((ray+bigray)*cosd(ang),0,(ray+bigray)*sind(ang));
  col5 = ((ray+bigray)*cosd(ang),0,-(ray+bigray)*sind(ang));
  col6 = (-len*sind(ang),0,len*cosd(ang));
  col7 = (-len*sind(ang),0,-len*cosd(ang));
  col8 = col1+col4+col6;
  col9 = col1+col5+col7;
  rigorousdisc( 0, true, col1+col3, ray, col7 );
  rigorousdisc( 0, true, col1+col2, ray, col6 );
  quartertorus( col1, col3, col2, ray );
  draw rp(col9)--rp(black)--rp(col8);
endfig;

verbatimtex 
\end{document} 
etex

end.