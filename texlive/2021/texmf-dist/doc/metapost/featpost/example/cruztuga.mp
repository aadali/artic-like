% cruztuga.mp
% L. Nobre G.
% 2007

prologues := 1;

beginfig(1);
  numeric u, da, db;
  u = 1cm;
  da = 1-1/sqrt(2);
  db = 3/sqrt(2)-2;
  z1 = u*(1.5,0.5);
  z2 = u*(1.5,2.5-da);
  z3 = u*(3.5+da,4.5);
  z4 = u*(-3.5-da,4.5);
  z5 = u*(db-1.5,3.5);
  z6 = u*(1.5-db,3.5);
  z7 = u*(0.5,2.5+db);
  z8 = 0.5u*(1,1);
  z10= 6u*(1,1);
  z11= 6u*(-1,1);
  z12= 6u*(-1,-1);
  z13= 6u*(1,-1);
  path octave;
  picture tmppic;
  octave = z1--z2--z3--z4--z5--z6--z7--z8--cycle;
  fill octave withcolor red;
  draw currentpicture xscaled -1;
  draw currentpicture rotated 90;
  draw currentpicture rotated 180;
  tmppic = currentpicture;
  currentpicture := nullpicture;
  fill z10--z11--z12--z13--cycle;
  draw tmppic;
endfig;

end.