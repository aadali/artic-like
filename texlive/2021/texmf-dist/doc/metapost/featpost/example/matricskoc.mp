% matricskoc.mp
% Pedro Mendonca
% 2006

input featpost3Dplus2D;

prologues:=2;

beginfig(1);
  ParallelProj := false;
  SphericalDistortion := false;
  Spread:=256;
  f:=(5,-2,2);
  cartaxes(-.5,.5,.5);
  for $=-8 step 1 until 0:
    kindofcube(true,true,($,0,0),90,0,0,.3,.05,0.1);
  endfor;
endfig;

end;


