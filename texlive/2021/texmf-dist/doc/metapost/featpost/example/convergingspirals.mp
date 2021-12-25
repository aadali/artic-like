% convergingspirals.mp
% L. Nobre G. - a twisted mind...
% 2007

prologues := 1;

vardef oper( expr Valu, Slop ) =
  (Valu*Slop)
enddef;

vardef hiperspiral( expr Val, Phi, Pitch ) =
  dir(57.29578*Val*(Pitch**2)+Phi)/oper(Val,Pitch)
enddef;

beginfig(5); %%%%%%%%%%%%%% This figure has one very small bug... %%%%%%%
  boolean adddepth;
  adddepth = true;  
  pi = 3.14159;
  co = 180/pi;
  mar = 0.20cm;
  prmax = 6cm;
  prmin = 0.85cm;
  hrmax = 7cm;
  ptmax = 1.699*pi;
  pstpnum = 100;
  numarquid = 5;
  multipl = 5; %%%%%%%%%%%%%%% must be odd for checkerboard %%%%%%%%%%%%%%
  checkerdepth = multipl+3;
  fac = hrmax/prmax;
  alfa = 2*pi/numarquid;
  pitcha = prmax/ptmax;
  pitchb = prmax/(ptmax-alfa);
  numhiperb = numarquid*multipl;
  ptstp = ptmax/pstpnum;
  ptstb = (ptmax-alfa)/pstpnum;
  umin = prmax/10;
  vmin = prmax/20;
  ustp = (hrmax-umin)/pstpnum;
  vstp = (hrmax-vmin)/pstpnum;

  path n[]m[], c[], cutter;
  cutter = fullcircle rotated 180 scaled (2prmin);
  for i=1 upto 8:
    n0m[i] = origin;
  endfor;
  for i=1 upto pstpnum:
    tp := i*ptstp;
    tb := i*ptstb;
    pu := umin+i*ustp;
    vu := vmin+i*vstp;
    n0m1 := n0m1--dir(co*tp)*(pitcha*tp+mar);
    n0m2 := n0m2--dir(co*tp)*pitcha*tp;
    n0m3 := n0m3--dir(co*(tb+alfa))*pitchb*tb;
    n0m4 := n0m4--dir(co*(tb+alfa))*(pitchb*tb-mar);
    n0m5 := n0m5--dir(co*pitcha/pu)*(pu-fac*mar);
    n0m6 := n0m6--dir(co*pitcha/pu)*pu;
    n0m7 := n0m7--dir(co*pitchb/vu)*vu;
    n0m8 := n0m8--dir(co*pitchb/vu)*(vu+fac*mar);
  endfor;
  pair flast, hlast;
  flast = point pstpnum of n0m6;
  hlast = point pstpnum of n0m8;
  adif = angle( hlast ) - angle( flast );
  n0m5 := n0m5 rotated adif;
  n0m6 := n0m6 rotated adif;
  for i=1 upto 8:
    n0m[i] := (reverse n0m[i]) cutafter cutter ;
  endfor;  
  for j=1 upto numarquid:
    for i=1 upto 4:
      n[j]m[i] := n0m[i] rotated (j*alfa*co);
%      draw n[j]m[i];
    endfor;
  endfor;
  for j=1 upto numhiperb:  
    for i=5 upto 8:
      n[j]m[i] := n0m[i] rotated (j*360/numhiperb);
%      draw n[j]m[i];
    endfor;
  endfor;
  distance = pstpnum;
  actualspik = 0;
  boolean intersectexists;
  pair itstm;
  for i=1 upto numhiperb:
    itstm := n1m2 intersectiontimes n[i]m6;
    intersectexists := ( itstm dotprod (-1,1) ) <> 0;
    if intersectexists:
      if ( (xpart itstm) < distance ):
	distance := xpart itstm;
	actualspik := i;
      fi;
    fi;
  endfor;
  numeric seq[];
  seq0 = actualspik;
  for i=1 upto 2checkerdepth:
    if (seq[i-1]-1<1):
      seq[i] := numhiperb;
    else:
      seq[i]=seq[i-1]-1;
    fi;
  endfor;
  path ca, cb;
  ca = (reverse n1m3)--n1m2;
  cb = (reverse n[seq1]m7)--n[seq1]m6;
  c0 = buildcycle( cb, reverse n1m1 ); %%%%% may reverse
  c1 = buildcycle( ca, n[seq0]m8 );
  for i=2 step 2 until 2checkerdepth:
    if i < 2*multipl-1 :
      ca := (reverse n[seq[i-2]]m6)--n[seq[i-2]]m7;
      c[i] := buildcycle( ca, n2m1 );
    else:
      c[i]:=buildcycle(n[seq[i-2]]m7,n2m1,reverse n[seq[i-2]]m6,reverse n0m4);
    fi;
    c[i+1]:=buildcycle(n[seq[i-2]]m5,n1m3,reverse n[seq[i]]m8,reverse n1m2);
  endfor;
  if adddepth:
    depthadded = 1;
  else:
    depthadded = 0;
  fi;
  for i=1 step 2 until 2*checkerdepth+depthadded:
    for j=1 upto numarquid:
      fill c[i-1] rotated (j*alfa*co) withcolor 0.85*(red+green);
      fill c[i] rotated (j*alfa*co) withcolor 0.7*green;
    endfor;
  endfor;
  if not (multipl <> 5):
    cb := (reverse n[seq3]m7)--n[seq3]m6;
    ca := buildcycle( cb, reverse n1m1 ); %%%%% may reverse
    for j=1 upto numarquid:
      fill ca rotated (j*alfa*co) withcolor 0.85*(red+green);
    endfor;
  fi;
endfig;

beginfig(4);
  numeric pstpnum, prmax, i, j, pi, ptmax, ptstp, th, co, mar;
  numeric numarquid, numhiperb, multipl, pitcha, pitchb, alfa;
  path a, b, c, d, e, f, g, h;
  numeric ptstb, hrmax, fac, sn, sm, rn, rm;
  pi = 3.14159;
  co = 180/pi;
  mar = 0.16mm;
  prmax = 5mm;
  hrmax = 7mm;
  ptmax = 1.05*pi;
  pstpnum = 100;
  numarquid = 5;
  multipl = 1;
  rn = 4;
  rm = 2;
  sn = 8;
  sm = 6;
  fac = hrmax/prmax;
  alfa = 2*pi/numarquid;
  pitcha = prmax/ptmax;
  pitchb = prmax/(ptmax-alfa);
  numhiperb = numarquid*multipl;
  ptstp = ptmax/pstpnum;
  ptstb = (ptmax-alfa)/pstpnum;
  a = origin;
  b = origin;
  c = origin;
  d = origin;
  for i=1 upto pstpnum:
    th := i*ptstp;
    a := a--dir(co*th)*pitcha*th;
    b := b--dir(co*th)*(pitcha*th+mar);
  endfor;
  for i=1 upto pstpnum:
    th := i*ptstb;
    c := c--dir(co*(th+alfa))*pitchb*th;
    d := d--dir(co*(th+alfa))*(pitchb*th-mar);
  endfor;
  for i=1 upto numarquid:
%    unfill buildcycle( b , reverse d ) rotated (i*alfa*co);
    fill a--reverse c--cycle rotated (i*alfa*co);
%     draw a rotated (i*alfa*co) withcolor red+green;
%     draw b rotated (i*alfa*co) withcolor red;
%     draw c rotated (i*alfa*co) withcolor green;
%     draw d rotated (i*alfa*co) withcolor blue;
  endfor;
  numeric umin, vmin, ustp, vstp, pu, pv;
  umin = prmax/20;
  vmin = prmax/20;
  ustp = (hrmax-umin)/pstpnum;
  vstp = (hrmax-vmin)/pstpnum;
  e = origin;
  f = origin;
  g = origin;
  h = origin;
  for i=1 upto pstpnum:
    pu := umin+i*ustp;
    e := e--dir(co*pitcha/pu)*(pu-((fac*mar)**2));
    f := f--dir(co*pitcha/pu)*pu;
    vu := vmin+i*vstp;
    g := g--dir(co*pitchb/vu)*(vu+sqrt(fac*mar));
    h := h--dir(co*pitchb/vu)*vu;
  endfor;
  pair flast, hlast;
  numeric adif;
  path hiperbolio, hiperbolic;
  flast = point pstpnum of f;
  hlast = point pstpnum of h;
  adif = angle( hlast ) - angle( flast );
%   draw e withcolor red+blue;
%   draw f withcolor red+blue;
%   draw g rotated (-adif) withcolor green+blue;  
%   draw h rotated (-adif) withcolor green+blue;
  hiperbolio = (subpath (rn,pstpnum) of g rotated (-adif))...
               (subpath (pstpnum,rm) of e)..cycle;
  hiperbolic = (subpath (sn,pstpnum) of h rotated (-adif))...
               (subpath (pstpnum,sm) of f)..cycle;
  for i=1 upto numhiperb:
    unfill hiperbolio rotated (i*360/numhiperb);
  endfor;
  for i=1 upto numhiperb:
    fill hiperbolic rotated (i*360/numhiperb) withcolor 0.4*red+0.2*blue;
  endfor;
  picture storepic;
  storepic = currentpicture;
  currentpicture := nullpicture;
  draw storepic scaled 12;
endfig;

beginfig(3);
  numeric pstpnum, prmax, i, j, pi, ptmax, ptstp, th, co, mar;
  numeric numarquid, numhiperb, multipl, pitcha, pitchb, alfa;
  path a, b, c, d;
  pi = 3.14159;
  co = 180/pi;
  mar = 0.32;
  prmax = 5;
  ptmax = 4*pi;
  pstpnum = 100;
  multipl = 3;
  numarquid = 3;
  alfa = 2*pi/numarquid;
  pitcha = prmax/ptmax;
  pitchb = prmax/(ptmax-alfa);
  numhiperb = numarquid*multipl;
  ptstp = ptmax/pstpnum;
  a = origin;
  b = origin;
  c = origin;
  d = origin;
  for i=1 upto pstpnum:
    th := i*ptstp;
    a := a...dir(co*th)*pitcha*th;
    b := b...dir(co*th)*(pitcha*th-mar);
    c := c...dir(co*(th+alfa))*pitchb*th;
    d := d...dir(co*(th+alfa))*(pitchb*th+mar);
  endfor;
  for i=1 upto numarquid:
    drawoptions( rotated (i*alfa*co) scaled 1000 withpen pencircle scaled 0pt);
    draw a;
    draw b;
    draw c;
    draw d;
  endfor;
  %%%%%%%%%%%%%%%%%%%%%%%% A very nice case of knitting Bezier splines...
endfig;

beginfig(2);
  numeric u, numb, faca, numa, numc, rmin, rmax, i, j, pitch, phi;
  numeric petch, startang, angstep, ampi, ampe, arcstep, paperfactor;
  path b[], base, a[], alto;
  numeric valmax, valmin, valstep, valhmin, valhmax;
  u = 15mm;
  numb = 2; % start at one
  faca = 7;
  numc = 9;
  rmin = 0.3;
  rmax = 1.5;
  paperfactor = 0.25*u/rmin;
  angstep = 1/30;
  valstep = 0.017;
  numa = numb*faca;
  startang = 1/numb;
  arcstep = 1/numa;
  ampi = numc/faca; 
  ampe = ampi-startang;
  pitch = (rmax-rmin)/ampi; %show pitch;
  petch = (rmax-rmin)/ampe; %show petch;
  valmin = 1/oper(rmax,pitch);
  valmax = 1/oper(rmin,pitch);
  valhmin = 1/oper(rmax,petch);
  valhmax = 1/oper(rmin,petch);
  for i=1 step 2 until 2*numb:
    phi := 0.5*i/numb;
    b[i-1] = dir(phi*360)*rmin for j=angstep step angstep until ampi:
      ...dir((phi+j)*360)*(rmin+pitch*j)
    endfor;
    b[i] = dir((phi+startang)*360)*rmin for j=angstep step angstep until ampe:
      ...dir((phi+startang+j)*360)*(rmin+petch*j)
    endfor;
  endfor;
  base = for i=1 step 2 until 2*numb: b[i-1]...reverse b[i]-- endfor cycle;
  fill base scaled paperfactor;
  for i=1 step 2 until 2*numa:
    phi := 0.5*(i-1)*arcstep*360;
    a[i-1] = hiperspiral(valmin,phi,pitch)
    for j=valmin+valstep step valstep until valmax:
      ...hiperspiral(j,phi,pitch)
    endfor;
    draw a[i-1] scaled paperfactor withcolor red+green;
    phi := 0.5*i*arcstep*360;
    a[i] = hiperspiral(valhmin,phi,petch)
    for j=valhmin+valstep step valstep until valhmax:
      ...hiperspiral(j,phi,petch)
    endfor;
    draw a[i] scaled paperfactor withcolor red;
  endfor;
endfig;

beginfig(1);
  numeric u, numb, faca, numa, numc, rmin, rmax, i, j, pitch, phi;
  numeric petch, startang, angstep, ampi, ampe, arcstep;
  path b[], base, a[], alto;
  u = 15mm;
  numb = 3; % start at one
  faca = 3;
  numc = 5;
  rmin = 2u;
  rmax = 10u;
  angstep = 1/90; % If too big, the figure may be broken...
  numa = numb*faca;
  ampi = numc/faca; 
  pitch = (rmax-rmin)/ampi;
  startang = 1/numb;
  arcstep = 1/numa;
  ampe = numc/faca-startang;
  petch = (rmax-rmin)/ampe;
  for i=1 step 2 until 2*numb:
    phi := 0.5*i/numb;
    b[i-1] = dir(phi*360)*rmin for j=angstep step angstep until ampi:
      ...dir((phi+j)*360)*(rmin+pitch*j)
    endfor;
    b[i] = dir((phi+startang)*360)*rmin for j=angstep step angstep until ampe:
      ...dir((phi+startang+j)*360)*(rmin+petch*j)
    endfor;
  endfor;
  base = for i=1 step 2 until 2*numb: b[i-1]...reverse b[i]-- endfor cycle;
  fill base withcolor red;
  for i=1 step 2 until 2*numa:
    phi := 0.5*i/numa;
    a[i-1] = dir(phi*360)*rmin for j=angstep step angstep until 0.5:
      ...dir(-(phi+j)*360)*(rmin+2*(rmax-rmin)*j)
    endfor;
    a[i] = dir((phi+arcstep)*360)*rmin
    for j=angstep step angstep until 0.5-arcstep:
      ...dir(-(phi+arcstep+j)*360)*(rmin+(rmax-rmin)/(0.5-arcstep)*j)
    endfor;
  endfor;
  alto = for i=1 step 2 until 2*numa: a[i-1]...reverse a[i]-- endfor cycle;
  fill alto withcolor green+blue;
endfig;

end.