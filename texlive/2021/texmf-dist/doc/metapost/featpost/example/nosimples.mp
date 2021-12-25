% nosimples.mp
% L. Nobre G.
% 2012

prologues := 1;

beginfig(1);
  numeric u, ray, side, i, pensize, ang, fac, mar;
  u = 2mm;
  ray = 2.9u;
  size = 10u;
  pensize = 0.45u;
  fac = 2.5;
  mar = 1.0;
  pickup pencircle scaled pensize;
  for i=1 upto 5:
    ang := i*72+18;
    z[i] = size*dir(ang);
  endfor;
  z6  = z1+ray*(dir(-36)+dir(54));
  z61 = z1+ray*(dir(-36)+dir(-18));
  z11 = z1+ray*(dir(216)+dir(126));
  z7  = z2+ray*(dir(216)+dir(126));
  z8  = z3+ray*(dir(216)+dir(198));
  z12 = z2+ray*(dir(216)+dir(198));
  z21 = z1+ray*(dir(216)+dir(198));
  z24 = z4+ray*(dir(216)+dir(198));
  z13 = z3+2*ray*down;
  z9  = z4+2*ray*down;
  z15 = z5+ray*(dir(-36)+dir(54));
  z10 = z5+ray*(dir(-36)+dir(-18));
  z14 = z4+ray*(dir(-36)+dir(-18));
  z63 = z3+ray*(dir(-36)+dir(-18));
  z16 = fac[z1,z2];
  z17 = (fac-mar)[z5,z3];
  z18 = (fac-mar)[z2,z4];
  z19 = fac[z1,z5];
  z35 = (z10--z14) intersectionpoint (z5--z19);
  z34 = (z13---z9..z14---z10) intersectionpoint (z4--z18);
  z33 = (z9---z13..z8---z12) intersectionpoint (z3--z17);
  z32 = (z2--z16) intersectionpoint (z8--z12);
  z31 = (z1--z2) intersectionpoint (z21--z24);
  z50 = (z21--z24) intersectionpoint (z3--z5);
  z55 = (z3--z5) intersectionpoint (z2--z4);
  z70 = (z61--z63) intersectionpoint (z2--z4);
  z75 = (z61--z63) intersectionpoint (z21--z24);
  draw z9---z13..z8---z12..z7;
  draw z55--z18--z19--z5;
  draw (z1{dir(144)}..z11..z21---z24..z9..{dir(54)}z14) cutafter (z2--z18);
  draw z35---z10..z15..z5---z50;
  draw z33..z17--z16--z32;
  draw z31---z1..z6--z15;
  draw z7--z11;
%  draw z1..z6..z61---z63..z13..z8;
  draw z70--z75;
%  dotlabels.top(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,55);
%  dotlabels.lft(21,24,31,32,33,34,35,50);
%  dotlabels.bot(61,63,70,75);
endfig;

end.