% sniperlogo.mp
% L. Nobre G.
% 2012

prologues := 1;

beginfig(1);
  numeric u, ray, siz, xdis, ang, mar, i;
  u = 1cm;
  siz = 5u;
  ray = siz/6;
  fill fullcircle scaled (2*(siz+ray));
  xdis = (siz-2*ray) +-+ (2*ray);
  ang = angle(xdis,2*ray);
  z1 = (siz-ray)*dir(ang);
  z2 = (siz-ray)*dir(90-ang);
  z3 = (ray,(siz-2*ray)*cosd(ang));
  z4 = z3 xscaled -1;
  z5 = z2 xscaled -1;
  z6 = z1 xscaled -1;
  for i=7 upto 16:
    z[i] = z[i-4] rotated 90;
  endfor;
  for i=0 upto 3:
    z[21+i] = (siz-ray)*dir(45+90*i);
  endfor;
  path sniperpath;
  sniperpath = z1{dir(ang+90)}..z21..{dir(180-ang)}z2..z3{down}..{up}z4..
    z5{dir(180+ang)}..z22..{dir(270-ang)}z6..z7{right}..{left}z8..
    z9{dir(270+ang)}..z23..{dir(-ang)}z10..z11{up}..{down}z12..
    z13{dir(ang)}..z24..{dir(90-ang)}z14..z15{left}..{right}z16..cycle;
  unfill sniperpath rotated 45;
%  draw sniperpath;
%  drawoptions( withcolor red );
%  dotlabels(1,2,3,4,5,6,7,8,9,21,22)
endfig;

end.