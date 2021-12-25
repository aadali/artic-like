% rosadosventos.mp
% L. Nobre G.
% 2013

prologues := 1;

def getcrossingscale( expr Ang, Lev ) =
  begingroup
    numeric angoflevel;
    pair tmpair, northpoi, levelpoi;
    angoflevel = 90-180/(2**Lev);
    northpoi = (0,1);
    levelpoi = dir(angoflevel);
    tmpair = northpoi+whatever*dir(Ang-90);
    tmpair = levelpoi+whatever*dir(angoflevel-Ang);
    ( cosd(Ang)*abs(levelpoi-tmpair) )
  endgroup
enddef;

beginfig(1);
  numeric level, i, j, ang, u, scal, penscal, startang;
  path basetriang, worktriang;
  color fillcol;
  pen pathpen;
  pair shiftvec;
  level = 4;
  ang = 18;
  u = 7cm;
  fillcol = (0.77,0.88,0.99);
  penscal = 50;
  pathpen = pencircle scaled (u/penscal);
  basetriang = (0,1)--(-sind(ang)/cosd(ang),0)--(sind(ang)/cosd(ang),0)--cycle;
  basetriang := basetriang scaled u;
  fill basetriang withcolor fillcol;
  fill basetriang rotated 180 withcolor fillcol;
  draw basetriang withpen pathpen;
  draw basetriang rotated 180 withpen pathpen;
  for i=1 upto level:
    scal := getcrossingscale( ang, i );
    shiftvec := (0,u*(1-scal));
    worktriang := basetriang scaled scal shifted shiftvec rotated -90;
    startang := -90+360/(2**(i+1));
    for j=startang step 360/(2**i) until 270:
      fill worktriang rotated j withcolor (fillcol/(2**i));
      draw worktriang rotated j withpen pathpen;
    endfor;
  endfor;
endfig;

end.

  