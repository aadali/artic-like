% twoholes.mp
% L. Nobre G. 
% 2004

input featpost3Dplus2D;

f := 1.5*(0,1.5,1);
Spread := 100;

beginfig(1);
  numeric stagelevel, holediam, outerdiam, stagemargin;
  numeric numberofstagesidesquares, stagesidesize;
  numeric centercoord;
  color centerposdisc, centerpostorus;
  stagelevel = 0.2;
  holediam = 0.6;
  stagemargin = 0.2;
  numberofstagesidesquares = 20;
  outerdiam = holediam+4*stagelevel;
  stagesidesize = 2*( outerdiam + stagemargin );
  centercoord = 0.3535*outerdiam;
  HoriZon := -stagelevel;
  centerpostorus = ( centercoord, -centercoord, 0 );
  centerposdisc = ( -centercoord, centercoord, HoriZon );
  setthestage( numberofstagesidesquares, stagesidesize );
  smoothtorus( centerpostorus, blue,
    0.5*holediam+stagelevel, stagelevel );
  rigorousdisc( 0.5*holediam, true, centerposdisc,
    0.5*outerdiam, 1.5*stagelevel*blue );
endfig;
end.
