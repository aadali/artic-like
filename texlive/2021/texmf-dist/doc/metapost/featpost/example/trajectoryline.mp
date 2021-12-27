% trajectoryline.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

Spread := 30;
f := (10,0,0);

color posa, posb, axer;
numeric kfac, refrad;
posa = -5*green;
posb = 5*green;
axer = blue;
kfac = 2.5;
refrad = 1;

def vecfunc( expr posit )=
  begingroup
    save modula, modulb, field;
    color field;
    numeric modula, modulb;
    modula = conorm( posit-posa );
    modulb = conorm( posit-posb );
    field = ( posit-posa )/(modula**3)+kfac*( posit-posb )/(modulb**3);
    ( -field )
  endgroup
enddef;

def magfunc( expr posit )=
  begingroup
    save perpd, refvec, dist;
    color perpd, refvec;
    numeric dist;
    perpd = ccrossprod( axer, posit );
    refvec = posit - refrad*ncrossprod( perpd, axer );
    dist = conorm( refvec );
    ( ncrossprod( refvec, perpd )/(dist**2) )
  endgroup
enddef;

beginfig(1);
  numeric i, lena, numa, sa, diffstep, ray;
  color stp, svl;
  path oneline;
  pen grossa, fina;
  grossa = pencircle scaled 3pt; 
  fina = pencircle scaled 2pt; 
  lena = 400;
  numa = 58;
  diffstep = 0.12;
  ray = 3.5;
  sa = 168;
  for i=sa step numa until (180):
    stp := ( 0, ray*cosd(i), ray*sind(i) );
    svl := N( -stp );
    oneline := trajectorypath( lena, stp, svl, diffstep, vecfunc );
    draw oneline withpen grossa;
    draw oneline withpen fina withcolor background;
  endfor;
  draw rp(posa) withpen grossa;
  draw rp(posb) withpen grossa;
endfig;

beginfig(2);
  numeric i, lena, numa, sa, diffstep, ray;
  color stp, svl;
  path oneline;
  pen grossa, fina;
  grossa = pencircle scaled 3pt; 
  fina = pencircle scaled 2pt; 
  lena = 100;
  numa = 4;
  diffstep = 1;
  ray = 5.5;
  sa = 97;
  for i=sa step numa until (120):
    stp := ( 0, ray*cosd(i), ray*sind(i) );
    svl := 0.8*N( -stp );
    oneline := magnetictrajectorypath( lena, stp, svl, diffstep, magfunc );
    draw oneline withpen grossa;
    draw oneline withpen fina withcolor background;
  endfor;
  draw rp(black) withpen grossa;
endfig;

end.

