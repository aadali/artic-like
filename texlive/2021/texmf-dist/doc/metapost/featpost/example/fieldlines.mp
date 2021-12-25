% fieldlines.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

Spread := 30;
f := (10,0,0);

color posa, posb;
numeric kfac;
posa = -5*green;
posb = 5*green;
kfac = 2.5;

def vecfunc( expr posit )=
  begingroup
    save modula, modulb, field;
    color field;
    numeric modula, modulb;
    modula = conorm( posit-posa );
    modulb = conorm( posit-posb );
    field = ( posit-posa )/(modula**3)+kfac*( posit-posb )/(modulb**3);
    ( field )
  endgroup
enddef;
  
beginfig(1);
  numeric i, lena, lenb, numa, numb, sa, sb, diffstep, ray;
  color stp;
  path oneline;
  pen grossa, fina;
  grossa = pencircle scaled 3pt; 
  fina = pencircle scaled 2pt; 
  lena = 105;
  lenb = 145;
  numa = 72;
  numb = 30;
  diffstep = 1;
  ray = 0.5;
  sa = 36;
  sb = 75;
  for i=sa step numa until (360+sa):
    stp := posa+( 0, ray*cosd(i), ray*sind(i) );
    oneline := fieldlinepath( lena, stp, diffstep, vecfunc );
    draw oneline withpen grossa;
%    draw oneline withpen fina withcolor background;
  endfor;
  for i=sb step numb until (360+sb):
    stp := posb+( 0, ray*cosd(i), ray*sind(i) );
    oneline := fieldlinepath( lenb, stp, diffstep, vecfunc );
    draw oneline withpen grossa;
%    draw oneline withpen fina withcolor background;
  endfor;
  draw rp(posa) withpen grossa;
  draw rp(posb) withpen grossa;
endfig;

end.
