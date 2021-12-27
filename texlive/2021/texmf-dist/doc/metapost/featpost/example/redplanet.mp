% redplanet.mp
% L. Nobre G.
% 2006

input featpost3Dplus2D;

f := (4,0,-1);
def cf( expr Pair ) =
  begingroup
    save xc, yc, xa, ya;
    numeric xc, yc, xa, ya;
    xc = xpart Pair;
    yc = ypart Pair;
    xa = 1.25*(xc+1);
    ya = 1.25*(yc+1);
    ( rp( ( cosd(ya)*cosd(xa), cosd(ya)*sind(xa), sind(ya) ) ) )
  endgroup
enddef;

beginfig(1);
  numeric u;
  color yellow, thered, theblue;
  path martelo, foice, monitor, ecran, mar, foi, rotfoi, rotmar, pa, pb, pc;
  u = 1pt;
  yellow = 0.87*(0.95*red+0.7*green);
  thered = 0.66red;
  theblue= 0.87*(0.95*blue+0.27*red);

  numeric cabosal, cabesal, cabolar, cabelar;
  numeric cabolen, cabeclen, cabecontrol, martelrot;
  pair martelpos;
  cabosal = 2u;
  cabesal = 6u;
  cabolar = 8u;
  cabelar = 10u;
  cabolen = 70u;
  cabeclen = 7u;
  cabecontrol = 4u;
  martelrot = 50;
  martelpos = dir(180)*0.03*cabolen;
  z100 = origin;
  z101 = ( cabesal, 0 );
  z102 = ( cabesal, cabosal+cabelar-cabolen );
  z103 = ( cabesal+cabolar, cabosal+cabelar-cabolen );
  z104 = ( cabesal+cabolar, 0 );
  z105 = ( cabesal+cabolar+cabeclen, 0 );
  z106 = ( cabesal+cabolar+cabecontrol, cabelar );
  z107 = ( cabesal+cabolar, cabelar );
  z108 = ( cabesal+cabolar, cabelar+cabosal );
  z109 = ( cabesal, cabelar+cabosal );
  z110 = ( cabesal, cabelar );
  z111 = ( 0, cabelar );
  martelo = z100--
  z101--
  z102--
  z103--
  z104--
  z105..controls z106 and 0.5[z106,z107]..
  z107--
  z108--
  z109--
  z110--
  z111--cycle;
  mar = martelo rotated martelrot shifted martelpos;
		 
  numeric fcabolen, fcabolar, foicecontrol, laminang;
  numeric foicerot, foicaber, flamilar, fcabomar, refang;
  numeric fmaomarg, fmaorais, fmaospac, fprotrud, basangul;
  pair foicepos;
  laminang = 9;
  refang = 10;
  basangul = 23;
  fcabomar = 3u;
  fcabolar = 8u - fcabomar;
  flamilar = 4u;
  fcabolen = 22u;
  foicecontrol = 2u;
  fmaomarg = 6u;
  fmaorais = 2u;
  fmaospac = 14u;
  fprotrud = 25u;
  foicaber = 42u;
  foicerot = -50;
  foicepos = dir(270)*0.325*foicaber;
  z200 = origin;
  z201 = ( -fcabomar, 0 );
  z202 = ( -fcabomar, -fcabolen );
  z203 = ( fcabolar, -fcabolen );
  z204 = ( fcabolar+foicecontrol, -fcabolen );
  z205 = ( fcabolar+fmaomarg, fmaorais-fmaospac );
  z206 = ( fcabolar, -fmaospac );
  z207 = ( fcabolar, 0 );
  z208 = ( -fprotrud, foicaber );
  z209 = ( 0, flamilar );
  foice = z200--
  z201--
  z202--
  z203..controls z204 and 0.5[z204,z205]..
  z205{dir(-160)}..{left}
  z206--
  z207{dir(basangul)}..tension 0.8..{dir(180+refang+laminang)}
  z208{dir(refang-laminang)}..tension 0.8..{dir(180+basangul)}
  z209--cycle;
%  draw foice;
%  dotlabels.top(203,204,205);
  foi = foice rotated foicerot shifted foicepos;

   numeric l, k;
%   k = length foi;
%   l = length mar;
%   show (k,l);
  rotfoi = startahead( foi, 6 );
  rotmar = startahead( mar, 1 );
% drawarrow subpath (0,1) of rotfoi withpen pencircle scaled 1pt withcolor red;
% drawarrow subpath (0,1) of rotmar withpen pencircle scaled 1pt withcolor red;

  fill rigorousfearpath( black, 1 ) withcolor thered;
  pa = crossingline(rotfoi,rotmar,0.001);
  pb = closedpathinspace( pa, 16, cf );
  pc = closedpathinspace( lasermachine(pa,3pt,0.3), 16, cf );
  fill pc withcolor theblue;
  fill pb withcolor yellow;
  
%  fill rotfoi withcolor yellow+thered;
%  fill rotmar withcolor 0.7blue;

  produce_auto_scale;

%   for l=0 step 15 until 360:
%     show angle(dir(l));
%   endfor;
endfig;


end;


