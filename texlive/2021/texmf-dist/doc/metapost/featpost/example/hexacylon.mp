% hexacylon.mp
% L. Nobre G.
% 2006

beginfig(1);
  numeric i, actang, hexdiam, curvray, myang;
  path outerp, innerp;
  hexdiam = 15cm;
  curvray = 8.5mm;
  myang = 36;
  for i=1 upto 6:
    actang := (i-1)*60-90;
    z[2*i-1] = 0.5*hexdiam*dir( actang ) + curvray*dir( actang - 120 );
    z[2*i]   = 0.5*hexdiam*dir( actang ) + curvray*dir( actang + 120 );
    z[20+i] = 0.5[z[2*i-1],z[2*i]] - curvray*dir( actang );
  endfor;
  outerp = for i=1 upto 6: z[2*i-1]...z[2*i]--- endfor cycle;
  innerp = z22---z23{dir(myang+180)}..{dir(myang+180)}z25---
           z26{dir(myang)}..{dir(myang)}cycle;
  fill outerp withcolor (0.87,0.12,0.08);
  fill innerp withcolor (0.08,0.12,0.87);
endfig;
end.