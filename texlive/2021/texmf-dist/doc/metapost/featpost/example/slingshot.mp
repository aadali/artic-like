% slingshot.mp
% L. Nobre G.
% 2008

input featpost3Dplus2D;

% bash: export TEX=latex
% tcsh: setenv TEX latex

%prologues := 0;

% verbatimtex
% \documentclass{article} 
% \usepackage{beton,concmath,ccfonts} 
% \begin{document} 
% etex

f := (30,0,0);
Spread := 100;
MaxFearLimit := 50;

numeric initdist, timestep, Gfactor, vplanet, vprobe;
numeric targetangle, cutlinangle;
inidist = 10;
timestep = 0.1;
Gfactor = 7;
vplanet = 1.3;
vprobe = 1.2;
numeric numsteps, cutlinelen, distfraction, numtries, diststep, i, ip;
numeric a, aa, ab, ac;
color inipos, inivel, dstp, iter;
path cutlinepath, traject;

def gravitfield( expr relatposprobe ) =
  ( -Gfactor*relatposprobe/((conorm( relatposprobe ))**3) )
enddef;

def costfu( expr actang ) =
  (((cutlinangle-actang-targetangle)/targetangle)**2)
enddef;

vardef scalarfu( expr index ) =
  traject:=trajectorypath(
    numsteps,
    inipos+index*dstp,
    inivel,
    timestep)(gravitfield);
  draw traject;
  ip := xpart( traject intersectiontimes cutlinepath );
  a := costfu( angle(direction ip of traject) );
  ( a )
enddef;

beginfig(1);
  targetangle = 50;
  cutlinelen = 15;
  distfraction = 0.74;
  numtries = 36;
  numsteps = 135;
  diststep = distfraction*inidist/numtries;
  dstp = (0,0,diststep);
  inipos = (0,0,-inidist);
  inivel = (0,-vplanet,vprobe);
  cutlinangle = angle( (-vplanet,vprobe) );
  %cartaxes( 1, 1, 1 );
  cutlinepath = rp(black)--rp(cutlinelen*N(inivel));
  draw cutlinepath withcolor red;
  drawoptions( withcolor blue );
  for i=0 upto numtries:
    show (i,scalarfu( i ));
  endfor;
  drawoptions();
  iter = (0,0.5*numtries,numtries);
  message "--------------------------------------------------------";
  forever:
    iter := minimizestep( iter )( scalarfu );
    exitif (Z(iter)-X(iter)<0.5) or (Y(iter)<0.005);
  endfor;
  draw traject withcolor green withpen (pencircle scaled (2.5mm)) ;
  produce_auto_scale;
endfig;

% verbatimtex 
% \end{document} 
% etex

end.
