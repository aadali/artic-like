% crossingline.mp
% L. Nobre G.
% 2007

input featpost3Dplus2D;

prologues := 1;

beginfig(3);
path foi, mar, foimar;
foi=
(-2.09207,0)..controls (-0.8113,1.52635) and (0.46948,3.05269)
 ..(1.75026,4.57904)..controls (16.50494,-7.80182) and (31.25961,-20.18265)
 ..(46.01431,-32.5635)..controls (47.72202,-30.52838) and (49.42972,-28.49324)
 ..(51.13742,-26.45813)..controls (36.38272,-14.07727) and (21.62805,-1.69644)
 ..(6.87337,10.68442)..controls (8.36761,12.46516) and (9.86185,14.2459)
 ..(11.3561,16.02664)..controls (1.80319,20.141) and (1.80319,20.141)
 ..(-0.75836,17.0883)..controls (-1.26715,17.51523) and (-1.77594,17.94215)
 ..(-2.28471,18.36908)..controls (-3.99242,16.33397) and (-5.70012,14.29883)
 ..(-7.40782,12.2637)..controls (-6.89905,11.83678) and (-6.39026,11.40985)
 ..(-5.88147,10.98293)..controls (-7.16225,9.45657) and (-8.44302,7.93024)
 ..(-9.7238,6.40388)..controls (-7.17989,4.26926) and (-4.63599,2.13463)
 ..cycle;
mar = 
(0,-12.55302)..controls (-0.6404,-11.78986) and (-1.28078,-11.02667)
 ..(-1.92117,-10.2635)..controls (-7.51778,-14.95969) and (-13.11438,-19.65587)
 ..(-18.71098,-24.35205)..controls (-17.00328,-26.38718) and (-15.29556,-28.42232)
 ..(-13.58786,-30.45744)..controls (-12.30708,-31.98378) and (-12.30708,-31.98378)
 ..(-2.1138,-28.6326)..controls (-3.95111,-27.57187) and (-6.11067,-26.96915)
 ..(-7.48247,-25.33434)..controls (-3.921,-22.34586) and (-0.35953,-19.35738)
 ..(3.20195,-16.3689)..controls (28.64638,-44.3316) and (76.49767,8.8977)
 ..(12.84161,37.23851)..controls (64.17233,7.60352) and (23.4439,-32.40083)
 ..(3.05269,-9.99147)..controls (2.03513,-10.84532) and (1.01756,-11.69917)
 ..cycle;
foimar = crossingline( mar, foi, 0.001 );
draw foimar;
show foimar;
endfig;

beginfig(2);
  numeric u, elhei, elwid, elpos, bawid, bahei, thstp, thdel, i;
  numeric refang, reft, tfin, tini, spikang;
  path el, ba, aux, cutel, horns;
  pair reftpair;
  color yellow, thered;
  pen penpoint, pendot;
  path togetp;
  yellow = 0.87*(0.95*red+0.7*green);
  thered = 0.66red;
  u = 2mm;
  thdel = 25;
  thstp = 5;
  elhei = 30u;
  elwid = 44u;
  elpos = 26u;
  bawid = 27u;
  bahei = 10u;
  penpoint = pencircle scaled 0.5u;
  pendot = pencircle scaled u;
  el = (elwid,elpos)..for i=thstp step thstp until 360-thstp:
    (elwid*cosd(i),elpos+elhei*sind(i)).. endfor cycle;
  ba = (bawid,-bahei){right}..{left}(bawid,0)--
    (elwid*cosd(90-thdel),elpos+elhei*sind(90-thdel))--
    (bawid/3.0,0)--(0,elpos+elhei)--(-bawid/3.0,0)--
    (elwid*cosd(90+thdel),elpos+elhei*sind(90+thdel))--
    (-bawid,0){left}..{right}(-bawid,-bahei)--cycle;
  aux = (0,elpos)--
    1.1[(0,elpos),(elwid*cosd(270+2thdel),elpos+elhei*sind(270+2thdel))];
  reftpair = el intersectiontimes aux;
  reft = xpart reftpair;
  refang = angle( direction reft of el );
  tfin = directiontime dir(180-refang) of el;
  tini = directiontime dir(180+refang) of el - length el;
  cutel = subpath (tini,tfin) of el;
  spikang = (90-angle((0,elpos+elhei)-(-bawid,0)));
  horns = cutel{dir(-refang-spikang)}..
    {left}(0,elpos-elhei+bawid/3.0)..
    {dir(refang+spikang)}cycle;
  draw ba;
  draw horns;
  drawarrow subpath (0,1) of ba withpen pendot withcolor red;
  drawarrow subpath (0,1) of horns withpen pendot withcolor red;
  togetp = crossingline( ba, horns, 0.001 );
  draw togetp withpen penpoint;
  show togetp;
endfig;

beginfig(1);
  numeric u, oray, odev, nl, circrot, i;
  pair vect;
  u = 1cm;
  oray = 1.1u;
  odev = 0.4*oray;
  nl = 5;
  vect = (oray,0);
  circrot = -90;
%  circrot = -180/nl;
  path outerp, innerp, togetp;
  innerp = fullcircle scaled 2u rotated circrot shifted vect;
  outerp = for i=0 upto nl-1:
    (oray-odev)*dir(i*360/nl){dir(i*360/nl+90)}...
    (oray+odev)*dir((i+0.5)*360/nl){dir((i+0.5)*360/nl+90)}...
    endfor cycle;
  togetp = crossingline( outerp, innerp, 0.001 );
  fill togetp;
  draw outerp withcolor red withpen pencircle scaled 1.2pt;
  drawarrow subpath (0,1) of outerp
            withcolor green withpen pencircle scaled 0.7pt;
  draw innerp withcolor blue;
  draw point 0 of innerp withcolor green withpen pencircle scaled 1pt;
endfig;
end.