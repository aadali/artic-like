% simplecar.mp
% L. Nobre G. 
% 2004

input featpost3Dplus2D;

%f := 0.25*(2,2.15,15.3);
f := (2,0.15,1.1);
Spread := 320;
HoriZon := -0.2;

beginfig(1);
  color refpos, angcol, lencol, froncol, rearcol;
  refpos = black;
  angcol = black;
  lencol = (0.6, 0.35, 0.28);
  froncol = (0.026, 0.12, 0.06);
  rearcol = froncol;
  setthearena( 12, 1.3 );
  simplecar( refpos, angcol, lencol, froncol, rearcol );
endfig;

end;


