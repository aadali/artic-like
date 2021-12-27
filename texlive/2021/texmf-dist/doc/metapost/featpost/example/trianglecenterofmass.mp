% trianglecenterofmass.mp
% L. Nobre G.
% 2010

%

prologues := 1;

beginfig(1);
  numeric u;
  u = 1cm;
  z1 = origin;
  z2 = (1u,0);
  z3 = (1+uniformdeviate(4))*dir(5+uniformdeviate(80))*u;
  draw z1--z2--z3--cycle withpen pencircle scaled 2pt;
  numeric a, b, c, d;
  a = (ypart z3)/3;
  z4 = (0,a);
  if (xpart z3) < (xpart z2):
    d = xpart z2;
  else:
    d = xpart z3;
  fi;
  z5 = (d,a);
  draw z4--z5 withcolor red;
  z12 = (xpart z3,0);
  z23 = whatever[z2,z3];
  z230 = unitvector((-ypart (z3-z2),xpart (z3-z2)));
  z23 = whatever[z1,z230];
  z13 = whatever[z1,z3];
  z130 = unitvector((ypart (z3-z1),-xpart (z3-z1)));
  z13 = whatever[z2,z2+z130];
  b = abs(z23)/3;
  c = abs(z13-z2)/3;
  z6 = z2+b*z230;
  z7 = z3+b*z230;
  draw z6--z7 withcolor green;
  z8 = z3+c*z130;
  z9 = z1+c*z130;
  draw z8--z9 withcolor blue;
endfig;

end.