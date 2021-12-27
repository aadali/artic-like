
beginfig(1);
  pair po, pa, pb, pc, pd, pe, pf, pg;
  numeric u, ang, bigray, smaray, i;
  path smacirc, heigray;
  u = 0.2mm;
  ang = 40;
  bigray = 250u;
  smaray = 205u;
  po = (105mm,145mm);
  pa = po - bigray*dir( ang );
  pb = po + (2*smaray-bigray)*dir( ang );
  pc = po + bigray*dir( ang );
  pe = po + bigray*dir( ang + 90 );
  pf = po + (smaray-bigray)*dir( ang );
  pg = po - bigray*dir( ang + 90 );
  heigray = po--pe;
  smacirc = (fullcircle scaled (2*smaray)) shifted pf;
  pd = heigray intersectionpoint smacirc;
  draw pg--pe;
  draw smacirc;
  draw (fullcircle scaled (2*bigray)) shifted po;
  draw pc--pa--pd--pb;
  draw pd--pe withcolor red;
  draw pb--pc withcolor red;
  dotlabel.bot( "O", po );
  dotlabel.llft( "A", pa );
  dotlabel.bot( "B", pb );
  dotlabel.urt( "C", pc );
  dotlabel.top( "D", pd );
  dotlabel.ulft( "E", pe );
  dotlabel.top( "F", pf );
endfig;

end.
  