% fallinthewind.mp
% L. Nobre G.
% 2014

input featpost3Dplus2D;

prologues := 1;

Spread := 9;
f := 4.5*(2.5,-1.3,1.7);
viewcentr := (0,0,2);
background := 0.9white;

color referencewind;
numeric dragcoef, dragcoeg, horangle, quadcoef;
horangle = 90;
referencewind = red*cosd( horangle ) + green*sind( horangle );
dragcoef = 0.065;
dragcoeg = 0.009;
quadcoef = 0.99;

def vecfunc( expr posit, veloc )=
  begingroup
    save heig, wind, field, wint, blow;
    color wind, field, blow;
    numeric heig, wint;
    heig = Z( posit );
    wind = referencewind*quadcoef*sqrt( heig+0.05 );
    wind := wind-veloc;
    wint = conorm( wind );
    blow = N( wind );
    field = (dragcoef+dragcoeg*wint)*wind-blue;
    ( field )
  endgroup
enddef;

beginfig(1);
  numeric i, j, k, numa, sa, ea, diffstep, ray;
  color stp, svl, arrbe, arren;
  path oneline[];
  pen grossa, fina;
  arrbe = (0,-8,5);
  arren = (0,-5.5,5);
  setthestage( 18, 20 );
  TDAtiplen := 0.75;      
  TDAhalftipbase := 0.3; 
  TDAhalfthick := 0.15;        
  grossa = pencircle scaled 1.8pt; 
  fina = pencircle scaled 1.2pt; 
  diffstep = 0.015;
  ray = 3.95;
  numa = 3.45;
  sa = 38;
  ea = 142;
  j = 0;
  for i=sa step numa until ea:
    j := j+1;
    stp := black;
    svl := ray*(green*cosd(i)+blue*sind(i));
    oneline[j] := dragtrajectorypath( stp, svl, diffstep, vecfunc );
    draw oneline[j] withpen grossa withcolor white;
  endfor;
  for k=1 upto j:
    draw oneline[k] withpen fina withcolor red;
  endfor;
  tdarrow( arrbe, arren );
  draw rp(black)--rp((0,0,4.5)) dashed withdots withpen fina withcolor blue;
  produce_vga_border;
endfig;

end.

ray := 3.25;
numeric min, rst;
min = 0.5;
rst = 0.1;
for i=min step rst until ray:
    %show i;
  numeric xa, xb, xc, fa, fb, fc;
  xa = 85;
  svl := i*(green*cosd(xa)+blue*sind(xa));
  oneline := dragtrajectorypath( black, svl, diffstep, vecfunc );
  fa = PrintStep; 
  xb = 89;
  svl := i*(green*cosd(xb)+blue*sind(xb));  
  oneline := dragtrajectorypath( black, svl, diffstep, vecfunc );
  fb = PrintStep;
  xc = 0;
  forever:
    exitif ( abs( fb - fa ) < 0.005 );
    xc := xb - fb*(xb-xa)/(fb-fa);
    svl := i*(green*cosd(xc)+blue*sind(xc));  
    oneline := dragtrajectorypath( black, svl, diffstep, vecfunc );
    fc := PrintStep;
    exitif ( abs( fc - fa ) < 0.0005 );
    xa := xb;
    fa := fb;
    xb := xc;
    fb := fc;
  endfor;
  message decimal( i ) & "   " & decimal( xc );
endfor;
message " ";



