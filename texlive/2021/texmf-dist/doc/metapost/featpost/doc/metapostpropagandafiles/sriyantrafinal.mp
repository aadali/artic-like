% sriyantrafinal.mp
% L. Nobre G.
% 2012

outputformat := "svg";
outputtemplate := "%j-%c.svg";

numeric u;
u = 5cm;

def ap( expr poia, poib ) =
  begingroup
    pair dirv;
    dirv = 2*u*unitvector(poib-poia);
    ( (poia-dirv)--(poia+dirv) )
  endgroup
enddef;

vardef makesriyantrafinal  =
  path outercircle;
  outercircle = fullcircle scaled (2*u);
  draw outercircle;
  numeric pA, pB, pC, pD, pE;
  pA = u*0.80556;
  pB = u*0.30064;
  pC = -u*0.17032;
  pD = u*0.31772;
  pE = u*0.51138;
  z1 = (pD,pA);   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  z2 = (-pD,pA);
  z3 = (0,pC);    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  z4 = ((-u,-pB)--(u,-pB)) intersectionpoint reverse outercircle; %%%%%%%%
  z5 = (-x4,-pB);
  z6 = (0,u);
  draw z4--z5--z6--cycle;
  z7 = (z1--z3) intersectionpoint (z4--z6);
  z8 = (-x7,y7);
  z9 = 0.5[z7,z8];
  z10 = 0.5[z4,z5];
  z11 = (0,-pA);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  y12 = y7;       
  x12 = pE;       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  z13 = whatever[z12,z11];
  z13 = whatever[z4,z6];
  z14 = ((-u,y13)--(u,y13)) intersectionpoint outercircle;
  z15 = (-x14,y13);
  z16 = (0,-u);
  draw z14--z15--z16--cycle;
  z17 = (-x12,y12);
  draw z12--z11--z17--cycle;
  z18 = whatever[z1,z3];
  z18 = whatever[z14,z15];
  z19 = whatever[z11,z12];
  z19 = whatever[z4,z5];
  z20 = 0.5[z14,z15];
  y21 = y11;
  z21 = whatever[z20,z19];
  z22 = (-x21,y21);
  draw z21--z22--z20--cycle;  
  z23 = whatever[z14,z16];
  z23 = whatever[z4,z5];
  z24 = whatever[z21,z20];
  z24 = whatever[z14,z16];
  x25 = 0;
  z25 = whatever[z23,z18];
  y26 = y25;
  z26 = whatever[z3,z7];
  z27 = (-x26,y26);
  draw z27--z26--z3--cycle;
  y28 = y24;
  z28 = whatever[z23,z25];
  z29 = (-x28,y28);
  draw z25--z28--z29--cycle;
  z30 = (z3--(u,pC)) intersectionpoint (z13--z11);
  z31 = (-x30,y30);
  draw z30--z31--z9--cycle;
  z32 = 0.5[z28,z29];
  z33 = (z3--z1) intersectionpoint (z30--z9);
  z34 = (-x33,y33);
  z35 = whatever[z25,z28];
  z35 = whatever[z33,z34];
  z36 = (-x35,y35);
  draw z32--z35--z36--cycle;
  z37 = (z21--z20) intersectionpoint (z1--z3);
  z38 = (-x37,y37);
  z39 = whatever[z37,z38];
  z39 = whatever[z30,z9];
  z40 = (-x39,y39);
  draw z10--z39--z40--cycle;
%  path pe, pf, pg;
%  pe = ap(z21,z20);
%  pf = ap(z35,z32);
%  pg = ap(z31,(u,y30));
%  filldraw buildcycle( pe, pf, pg ) withcolor red;
enddef;

beginfig(1);
  drawoptions( withpen pencircle scaled 3pt );
  makesriyantrafinal;
endfig;

end.