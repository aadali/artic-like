% antimattermeteor.mp
% L. Nobre G.
% 2014

input featpost3Dplus2D;

prologues := 1;

f := (30,0,0);
Spread := 100;

numeric initdist, timestep, Gfactor, vprobe;
inidist = 10;
timestep = 0.1;
Gfactor = 2;
vprobe = 2.2;
numeric numsteps, numtries, diststep, i, sigma, sigmang, angl;
color inipos, inivel, dstp, veloc;
path traject;
pen bluepen, redpen;
bluepen = pencircle scaled 0.9mm;
redpen = pencircle scaled 0.55mm;

def gravitfield( expr relatposprobe ) =
  ( Gfactor*relatposprobe/((conorm( relatposprobe ))**3) )
enddef;

vardef parabol( expr index ) =
  traject:=trajectorypath(
    numsteps,
    inipos+index*dstp,
    inivel,
    timestep)(gravitfield);
  draw traject;
enddef;

vardef hiperbo( expr veloc ) =
  traject:=trajectorypath(
    numsteps,
    inipos,
    veloc,
    timestep)(gravitfield);
  draw traject;
enddef;

beginfig(1);
  numtries = 36;
  numsteps = 95;
  sigma = 3.15;
  sigmang = 16;
  diststep = sigma/numtries;
  dstp = (0,0,diststep);
  inipos = (0,-inidist,0);
  inivel = (0,vprobe,0);
  drawoptions( withpen bluepen withcolor blue );
  for i=-numtries upto numtries:
    parabol( i );
  endfor;
  drawoptions( withpen redpen withcolor red );
  for i=-numtries upto numtries:
    angl := i*sigmang/numtries;
    veloc := vprobe*(0,cosd(angl), sind(angl));
    hiperbo( veloc );
  endfor;
  drawoptions();
  fill fullcircle scaled (4cm) shifted (rp(black));
  produce_auto_scale;
endfig;

end.
