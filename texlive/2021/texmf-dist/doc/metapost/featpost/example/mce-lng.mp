% mce-lng.mp
% L. Nobre G.
% 2009

prologues := 1;

beginfig(1);
  numeric u, anga, angb, posa, posb, thic, marg, angm, stwo, rest;
  numeric plac, n;
  path pta, ptb, ptc, ptd, pte, ptf, ptfa, ptfb, pa, pb;
  picture trib, half, unib;
  transform ta, tb;
  u = 10cm;
  anga = 120;
  angb = 135;
  posa = 0.8;
  posb = 0.23;
  thic = 0.15;
  marg = 0.07;
  angm = 10;
  plac = 0.5;
  stwo = u*sqrt(2);
  rest = u/sqrt(2);
  z1  = u*(0,posa+marg);
  z2  = u*(0,posa);
  z3  = u*(0,posa-thic);
  z4  = u*(0,posa-thic-marg);
  z6  = u*(0,-posb+0.5*thic);
  z7  = u*(0,-posb-0.5*thic);
  z8  = (posb-0.5*thic)[(0,u),(u,0)];
  z9  = (posb+0.5*thic)[(0,u),(u,0)];
  z12 = posb[(0,u),(u,0)];
  z11 = z12+rest*(marg,marg);
  z13 = z12-rest*(thic,thic);
  z14 = z12-rest*(marg+thic,marg+thic);
  z22 = posa[(0,-u),(u,0)];
  z21 = z22+rest*(marg,-marg);
  z23 = z22-rest*(thic,-thic);
  z24 = z22-rest*(marg+thic,-marg-thic);
  z33 = (plac*x43,y43);
  z34 = (plac*x43+stwo*marg,y43);
  z32 = (plac*x43-stwo*thic,y43);
  z31 = (plac*x43-stwo*(thic+marg),y43);
  z43 = (posa-thic)[(0,-u),(u,0)];
  pta = z6{dir(angb+0.5*angm-90)}..{dir(angb+0.5*angm-45)}z8;
  ptb = z7{dir(angb-0.5*angm-90)}..{dir(angb-0.5*angm-45)}z9;
  ptc = z4{dir(anga-angm-90)}..{dir(-45)}z14..{dir(-135)}z24..z34;
  ptd = z3{dir(anga-angm-90)}..{dir(-45)}z13..{dir(-135)}z23..{up}z33
  ..{dir(anga-angm-135)}z43;
  pte = z2{dir(anga-90)}..{dir(-45)}z12..{dir(-135)}z22..{up}z32
  ..{dir(anga-135)}z22;
  ptf = z1{dir(anga-90)}..{dir(-45)}z11..{dir(-135)}z21..{up}z31;
  ta = identity scaled (1/sqrt(2)) shifted (0,rest) rotated 45; 
  tb = identity scaled (1/sqrt(2)) shifted (0,-rest) rotated -45; 
  ptfa = ptf transformed ta; 
  ptfb = ptf transformed tb;
  pa = buildcycle( pte, ptc, reverse ptd, reverse ptfa );
%   draw pta;
%   draw ptb;
%   draw ptc;
%   draw ptd;
%   draw pte;
%   draw ptf;
%   draw ptfb;
  pb = buildcycle( pta, ptc, reverse ptb, reverse ptfb );
  fill pa;
  fill pb;
%  draw (0,-u)--(u,0)--(0,u)--cycle;
  trib = currentpicture;
  currentpicture := nullpicture;
  drawoptions( shifted (u,0) );
  draw trib;
  draw trib transformed ta;
  draw trib transformed tb;
  drawoptions();
  draw trib scaled 0.5 shifted (0.5u,0) rotated 90;
  draw trib scaled 0.5 shifted (0.5u,0) rotated (-90);
  half = currentpicture;
  draw half rotated 180;
  unib := currentpicture;
  draw currentpicture scaled 0.5;
  draw unib scaled 0.25;
%  draw currentpicture scaled 0.25;
%  draw unib scaled 0.0625 withcolor blue;
%  dotlabels.urt( 1,2,3,4,6,7,8,9,11,12,13,14,21,22,23,24,31,32,33,34,43);
%  dotlabels.urt(31,32,33,34,43);
endfig;

end.