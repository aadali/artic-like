% profaux.mp
% L. Nobre G.
% 2005

beginfig(1);
  numeric u, penbig, pensma, aangle, barhei, cirray, ularg, umarg, xhei, xcen;
  numeric atray;
  u = 1cm;
  cirray = 5u;
  atray  = 6.7u;
  aangle = 33;
  penbig = 1.85u;
  pensma = 0.75u;
  barhei = -1u;
  ularg = penbig-pensma;
  umarg = 0.5*pensma;
  color cux, coa;
  cux = 0.78red;
  coa = 0.78green;
  pen bigpen, smapen;
  bigpen = pencircle scaled penbig;
  smapen = pencircle scaled pensma;
  path theu;
  z11 = (-0.5*ularg,0.5*ularg+umarg);
  z12 = ( 0.5*ularg,0.5*ularg+umarg);
  theu = (z11--(halfcircle scaled ularg rotated 180)--z12);
  draw theu;
  xhei = ularg+umarg;
  xcen = 0.5*umarg;
  z13 = (ularg+1.5*pensma)*right + xcen*up;
  z14 = z13 + 0.5*(ularg,xhei);
  z15 = z13 - 0.5*(ularg,xhei);
  z16 = z13 + 0.5*(-ularg,xhei);
  z17 = z13 - 0.5*(-ularg,xhei);
  draw (z14--z15);
  draw (z16--z17);
  picture tmpic;
  tmpic = currentpicture;
  currentpicture := nullpicture;
  path horba, theci, lsid, rsid;
  theci = fullcircle scaled (2*cirray);
  z2 = atray*up;
  z1 = z2 - 2*atray*dir(90-0.5*aangle);
  z3 = z2 + 2*atray*dir(0.5*aangle-90);
  z5 = barhei*up;
  z4 = z5 + cirray*left;
  z6 = z5 + cirray*right;
  z18= (z2--z1) intersectionpoint theci;
  z19= (z2--z3) intersectionpoint theci;
  lsid = z18--z1;
  rsid = z19--z3;
  z7 = (z4--z5) intersectionpoint lsid;
  z8 = (z5--z6) intersectionpoint rsid;
  horba = z7--z8;
  z9 = (z7--z1) intersectionpoint theci;
  z10= (z8--z3) intersectionpoint theci;
  drawoptions( withpen bigpen withcolor coa );
  draw theci;
  draw z18--z9;
  draw z19--z10;
  draw horba;
  draw tmpic shifted (z10+umarg*(right+down)) withpen smapen withcolor cux;
%  drawoptions( withcolor blue );
%  dotlabels.top(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19);
endfig;
end.



  