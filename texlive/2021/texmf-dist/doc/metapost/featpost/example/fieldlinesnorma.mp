% fieldlinesnorma.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

Spread := 25;
f := (10,-2,7);

color posa, posb, perpd;
numeric kfac;
posa = -5*green;
posb = 5*green;
perpd = (1,0,1);
kfac = -10;

def vecfunc( expr posit )=
  begingroup
    save modula, modulb, field;
    color field;
    numeric modula, modulb;
    modula = conorm( posit-posa );
    modulb = conorm( posit-posb );
    field = ( posit-posa )/(modula**3)+kfac*( posit-posb )/(modulb**3);
    ( N(field) )
  endgroup
enddef;

def potfunc( expr posit )=
  begingroup
    save field;
    color field;
    field = vecfunc( posit );
    ( ncrossprod( field, perpd ) )
  endgroup
enddef;
  
beginfig(1);
  numeric i, lena, lenb, numa, numb, sa, sb, diffstep, ray;
  color stp;
  path oneline;
  pen grossa, fina;
  grossa = pencircle scaled 3pt; 
  fina = pencircle scaled 2pt; 
  lena = 85;
  lenb = 20;
  numa = 30;
  numb = 18;
  diffstep = 0.1;
  ray = 0.5;
  sa = 75;
  sb = 90;
  for i=sa step numa until (360+sa):
    stp := posa+( ray*cosd(i), ray*sind(i), 0 );
    oneline := fieldlinepath( lena, stp, diffstep, vecfunc );
    draw oneline withpen grossa;
%    draw oneline withpen fina withcolor background;
  endfor;
  for i=sb step numb until (360+sb):
    stp := posb+( ray*cosd(i), ray*sind(i), 0 );
    oneline := fieldlinepath( lenb, stp, -diffstep, vecfunc );
    draw oneline withpen grossa;
%    draw oneline withpen fina withcolor background;
  endfor;
  draw rp(posa) withpen grossa;
  draw rp(posb) withpen grossa;
endfig;

Spread := 17;
f:= (10,-6,7);

beginfig(2);
  numeric i, lena, lenb, numa, numb, sa, sb, diffstep, ang;
  color stp;
  path oneline;
  pen grossa, fina;
  grossa = pencircle scaled 3pt; 
  fina = pencircle scaled 2pt; 
  lena = 85;
  lenb = 125;
  numa = 0.03;
  numb = 0.30;
  diffstep = 0.2;
  ang = -90;
  sa = 1.5;
  sb = 4.5;
  for i=sa step numa until (1.5+sa):
    stp := posa+i*cosd(ang)*green+i*sind(ang)*ncrossprod(perpd,green);
    oneline := fieldlinepath( lena, stp, diffstep, potfunc );
    draw oneline withpen grossa;
%    draw oneline withpen fina withcolor background;
  endfor;
  ang := 180;
  for i=sb step numb until (3.7+sb):
    stp := posb+i*cosd(ang)*green+i*sind(ang)*ncrossprod(perpd,green);
    oneline := fieldlinepath( lenb, stp, diffstep, potfunc );
    draw oneline withpen grossa;
%    draw oneline withpen fina withcolor background;
  endfor;
  draw rp(posa) withpen grossa;
  draw rp(posb) withpen grossa;
endfig;

end.
