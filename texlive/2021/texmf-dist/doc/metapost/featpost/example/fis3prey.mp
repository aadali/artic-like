% fis3prey.mp
% L. Nobre G. 
% 2006

input featpost3Dplus2D;

verbatimtex
\documentclass{article}
\usepackage{palatino,mathpazo}
\begin{document}
etex

beginfig(1);
  Spread := 50;
  f := (4,6,3.5);
  path malhacub, rama, ramb, ramc, malhatet;
  pen grossa, nuclear, negat;
  color p[];
  numeric lengax, lscal, arrowmarg, arrowshift;
  pair tx, ty, tz;
  lengax = 0.7;
  lscal = 2;
  arrowshift = 0.85;
  arrowmarg = 0.13;
  grossa = pencircle scaled 2.7pt;
  nuclear = pencircle scaled 2.1pt;
  negat = pencircle scaled 1.8pt;
  p1 = (1,1,-1);
  p2 = (-1,1,-1);
  p3 = (-1,-1,-1);
  p4 = (1,-1,-1);
  p5 = (1,1,1);
  p6 = (-1,1,1);
  p7 = (-1,-1,1);
  p8 = (1,-1,1);

  p0 = p2 + arrowshift*N((1,0,1)) + arrowmarg*N((-1,0,1));
  p9 = p5 - arrowshift*N((1,0,1)) + arrowmarg*N((-1,0,1));

  z1 = rp(p1);
  z2 = rp(p2);
  z3 = rp(p3);
  z4 = rp(p4);
  z5 = rp(p5);
  z6 = rp(p6);
  z7 = rp(p7);
  z8 = rp(p8);

  z0 = rp(p0);
  z9 = rp(p9);
  
  tx = rp( p1 - (lengax,0,0) ); 
  ty = rp( p1 - (0,lengax,0) ); 
  tz = rp( p1 + (0,0,lengax) ); 

  malhacub = z1--z2--z3--z4--z1--z5--z6--z7--z8--z5;
  malhatet = z2--z7--z4--z5--cycle;
  rama = z2--z6;
  ramb = z3--z7;
  ramc = z4--z8;
  drawoptions( dashed evenly );
  draw malhacub;
  draw rama;
  draw ramb;
  draw ramc;
  drawoptions( withpen grossa );
  draw malhatet;

  drawarrow z1--tx withpen negat;
  drawarrow z1--ty withpen negat;
  drawarrow z1--tz withpen negat;
  %show "alert for fis3prey.mp";
  label.lrt(btex $x$ etex scaled lscal, tx);
  label.llft(btex $y$ etex scaled lscal, ty);
  label.lft(btex $z$ etex scaled lscal, tz);
  drawarrow z0--z9 withpen negat;
endfig;

beginfig(2);
  Spread := 50;
  f := 2.3*(6,2,3.5);
  numeric lang, cray, lpos, llen, smar;
  path circ, rigcirc, lone, ltwo, carrpath;
  pen pena, penb, penc;
  pena = pencircle scaled 1.3mm;
  penb = pencircle scaled 0.95mm;
  penc = pencircle scaled 0.4mm;
  color outcol, incol;
  outcol = black;
  incol = 0.876white;
  lang = 30;
  cray = 1;
  lpos = 0.57;
  llen = 2.2;
  smar = 0.15;
  circ = goodcirclepath( black, green, cray );
  rigcirc = rigorouscircle( black, green, cray+smar );
  carrpath = subpath (3.2,4.5) of rigcirc;
  lone = rp((-lpos,llen*cosd(lang),llen*sind(lang)))
       --rp(-(lpos,llen*cosd(lang),llen*sind(lang)));
  ltwo = rp((lpos,-llen*cosd(lang),llen*sind(lang)))
       --rp((lpos,llen*cosd(lang),-llen*sind(lang)));
  drawarrow carrpath withpen penc withcolor outcol;
  draw circ withpen pena withcolor outcol;
  draw circ withpen penb withcolor incol;
  draw rp((-cray,0,0))--rp((-lpos,0,0)) dashed evenly;
  draw lone withpen pena withcolor outcol;
  draw lone withpen penb withcolor incol;
  draw rp((lpos,0,0))--rp((-lpos,0,0)) dashed evenly;
  draw rp((0,0,-cray))--rp((0,0,cray)) dashed evenly;
  draw ltwo withpen pena withcolor outcol;
  draw ltwo withpen penb withcolor incol;
  draw rp((cray,0,0))--rp((lpos,0,0)) dashed evenly;
  draw subpath (1,15) of circ withpen pena withcolor outcol;
  draw subpath (0,16) of circ withpen penb withcolor incol;
%  label.top( btex $+$ etex, point 3 of rigcirc );
endfig;

beginfig(3);
  TDAtiplen := 0.26;              % tdarrow
  TDAhalftipbase := 0.12;         % Three-Dimensional
  TDAhalfthick := 0.05;           % Arrow
  Spread := 23;
  numeric len, rad, tle, tcr, tof;
  len = 10;
  rad = 0.15;
  tle = 1.5;
  tcr = 4;
  tof = 0.35;
  f := 10*(1.7,0.8,1.05);
  tdarrow( (tcr,-1-tof,-1), (tcr-tle,-1-tof,-1) );
  rigorousdisc( 0, true, (-2*len,-1,-1), rad, 3*len*red );
  tdarrow( (tcr,1-tof,-1), (tcr-tle,1-tof,-1) );
  draw rp((-1,-1+rad,-1))--rp((-1,1-rad,-1)) dashed evenly;
  draw rp((1,-1+rad,-1))--rp((1,1-rad,-1)) dashed evenly;
  rigorousdisc( 0, true, (-2*len,1,-1), rad, 3*len*red ); 
  tdarrow( (-1,-tcr,1+tof), (-1,-tcr+tle,1+tof) );
  draw rp((-1,-1,-1+rad))--rp((-1,-1,1-rad)) dashed evenly;
  draw rp((-1,1,-1+rad))--rp((-1,1,1-rad)) dashed evenly;
  draw rp((1,-1,-1+rad))--rp((1,-1,1-rad)) dashed evenly;
  draw rp((1,1,-1+rad))--rp((1,1,1-rad)) dashed evenly;
  rigorousdisc( 0, true, (-1,-len,1), rad, 2*len*green ); 
  draw rp((-1+rad,-1,1))--rp((1-rad,-1,1)) dashed evenly;
  draw rp((-1+rad,1,1))--rp((1-rad,1,1)) dashed evenly;
  rigorousdisc( 0, true, (1,-len,1), rad, 2*len*green ); 
  tdarrow( (1,-tcr,1-tof), (1,-tcr+tle,1-tof) );
  label( btex $I_A$ etex, rp((tcr-0.5*tle,-1-2*tof,-1)) );
  label( btex $I_B$ etex, rp((tcr-0.5*tle,1-2*tof,-1)) );
  label( btex $I_C$ etex, rp((-1,-tcr+0.5*tle,1+2*tof)) );
  label( btex $I_D$ etex, rp((1,-tcr+0.5*tle,1-2*tof)) );
  produce_vga_border;
endfig;

verbatimtex
\end{document}
etex

end;
