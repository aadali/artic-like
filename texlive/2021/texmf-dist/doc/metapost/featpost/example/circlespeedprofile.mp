% circlespeedprofile.mp
% L. Nobre G.
% 2012

prologues := 1;

beginfig(1);
  numeric u, s, a, ai, as, af, dist, b, ray;
  pen one;
  pair poi, ref, tip, cen, contactpoi, toppoi;
  u = 10cm;
  s = 0.5;
  ai = 270;
  as = -8;
  af = -90-as;
  one = pencircle scaled 0.015u;
  draw (-0.5u*(1+s,1))--(u,-0.5u) withpen one withcolor 0.5white;
  draw fullcircle scaled u withpen one;
  contactpoi = (0,-0.5u);
  toppoi = (0,0.5u);
  cen = 0.5*(contactpoi+toppoi+(s*u,0));
  ray = abs( cen-toppoi );
  draw fullcircle scaled 2ray shifted cen withpen one withcolor red;
  ref = (0,-0.5u);
  ahangle := 30;
  ahlength := 0.03u;
  for a=ai step as until af:
    poi := 0.5u*dir(a);
    dist := abs( poi-ref );
    tip := poi+s*dist*dir(0.5*a-45);
    drawarrow poi--tip;
  endfor;
  % path prof;
  % prof = for a=ai step as until af:
  %     hide( poi := 0.5u*dir(a);
  %     dist := abs( poi-ref);
  %     tip := (poi+s*dist*dir(0.5*a-45));
  %     drawarrow poi--tip;)
  %     tip.. endfor cycle;
  %  draw prof withpen one withcolor red;
endfig;
 
end.