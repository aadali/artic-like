% electrospiral.mp
% L. Nobre G. 
% 2004

numeric turnstep;
turnstep = 1cm;

def spiral( expr Ang ) =
  begingroup
    numeric k, r;
    k = turnstep/360.0;
    r = k*Ang;
    ( r*dir(Ang) )
  endgroup
enddef;
    
beginfig(1);
  numeric i, angstep, actx, ray, corner, sider;
  path spir, cutter;
  pair thispos, sw, nw, ne, se, actpos;
  angstep = 30;
  ray = 0.7mm;
  corner = 4.2;
  sider = 10.5;
  spir = origin..
  for i=angstep step angstep until 4000-1:
    spiral(i)..
  endfor
  spiral(4000);
  draw spir;
  actx = 1;
  forever:
    thispos := (turnstep*actx,0);
    for i=1 upto actx:
      actpos := thispos+ray*normaldeviate*dir(uniformdeviate(359));
      draw actpos withpen pencircle scaled 1.5mm;
    endfor;
    actx := 2*actx;
    exitif actx > 9;
  endfor;
  sw = turnstep*(-corner,-corner);
  nw = turnstep*(-corner,corner);
  ne = turnstep*(sider,corner);
  se = turnstep*(sider,-corner);
  cutter = sw--nw--ne--se--cycle;
  clip currentpicture to cutter;
endfig;

end;

