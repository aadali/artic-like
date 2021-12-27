% rothexagrid.mp
% L. Nobre G.
% IYP (2005)

vardef locator( expr Indi, Indj ) =
   ((Indi+0.5*Indj)*5mm,0.5*sqrt(3)*Indj*5mm)
enddef;

vardef gridmark( expr Loc ) =
  draw Loc--(Loc+sqrt(3)*down*5mm/3);
  draw Loc--(Loc+( 0.5,sqrt(3)/6)*5mm);
  draw Loc--(Loc+(-0.5,sqrt(3)/6)*5mm);
enddef;

beginfig(1);
  numeric i, j, fdist;
  picture tmppic;
  path clipper;
  color forec, backc;
  transform T;
  fdist = 10cm;
  forec = (0.8,1,0.8);
  backc = (0.8,0.8,0.8);
  z1 = locator( 0, 11 );
  z2 = locator( 4, 22 );
  z3 = locator( 19, 11 );
  z4 = locator( 15, 0 );
  origin = z1 transformed T;
  (0,fdist) = z2 transformed T;
  (fdist,fdist) = z3 transformed T;
  clipper = (z1--z2--z3--z4--cycle) transformed T ;
  drawoptions( withpen pencircle scaled 1pt );
  for i = 0 upto 19:
    for j = 0 upto 22:
      gridmark( locator( i, j ) );
    endfor;
  endfor;
  tmppic = currentpicture;
  currentpicture := nullpicture;
  fill clipper withcolor backc;
  draw tmppic transformed T withcolor forec;
  clip currentpicture to clipper;
%  dotlabels.top( 1,2,3,4 );
endfig;
end.