% twoDcolision.mp
% L. Nobre G.
% 2013

verbatimtex \documentclass{article}\begin{document} etex;

beginfig(1);
  ahangle := 30;
  numeric angthe, angalf, xpoin;
  numeric radthe, vini;
  numeric ysize, xsize, margs, units, i, dim;
  pair vone, vinic, vtwo;
  path xaxis, yaxis, vonep, vinicp, vtwop, auxil;
  pen thick, usual, vthin;
  angthe = 65;
  dim = 50;
  units = 0.8mm;
  ahlength := 2.5units;
  margs = 8units;
  xsize = 75units;
  ysize = 40units;
  radthe = 23units;
  vini = xsize-2.5*margs;
  thick = pencircle scaled 0.5units;
  usual = pencircle scaled 0.2units;
  vthin = pencircle scaled 0.1units;
  yaxis = ( 0, -vini/2 )--( 0, ysize-margs );
  xaxis = ( -margs, 0 )--( xsize-margs, 0 );
  pickup usual;
  drawarrow xaxis;
  label.bot( btex $x$ etex, point 1 of xaxis );
  drawarrow yaxis;
  label.lft( btex $y$ etex, point 1 of yaxis );
  pickup thick;
  vinic = vini*right;
  vinicp = origin--vinic;
  drawarrow vinicp;
  label.urt( btex $\vec{v}_i$ etex, vinic );
  vone = vini*cosd(angthe)*dir(angthe);
  vonep = origin--vone;
  drawarrow vonep; 
  label.ulft( btex $\vec{v}_1$ etex, vone );
  vtwo = vini*sind(angthe)*dir(angthe-90);
  vtwop = origin--vtwo;
  drawarrow vtwop; 
  label.lrt( btex $\vec{v}_2$ etex, vtwo );
  pickup vthin;
  draw fullcircle scaled radthe cutafter vonep;
  draw fullcircle scaled vini shifted (vinic/2);
  label.urt( btex $\Theta$ etex, ( radthe/2, 0 ) );
  z1 = dir(angthe)*margs/2;
  z2 = dir(angthe-45)*margs/sqrt(2);
  z3 = dir(angthe-90)*margs/2;
  draw z1--z2--z3;
endfig;

verbatimtex \end{document} etex

end;
