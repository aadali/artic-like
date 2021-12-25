% fis3prex.mp
% L. Nobre G. 
% 2004

input featpost3Dplus2D;

verbatimtex
\documentclass{article}
\usepackage{beton}
\usepackage{concmath}
\usepackage{ccfonts}
\begin{document}
etex

def platecap( expr Xpos, Shade, Angle, Oray, Iray, NumPoints ) =
  begingroup
    path perim, op, ip;
    numeric i;
    color axey, axmz, tp;
    axey = cosd( Angle )*green + sind( Angle )*blue;
    axmz = -cosd( Angle )*blue + sind( Angle )*green;
    tp = Xpos*red;
    op = rp( tp+Oray*axey )
      for i=1 upto NumPoints:
      ..rp(tp+Oray*(cosd(i*180/NumPoints)*axey+sind(i*180/NumPoints)*axmz))
      endfor;
    ip = rp( tp+Iray*axey )
      for i=1 upto NumPoints:
      ..rp(tp+Iray*(cosd(i*180/NumPoints)*axey+sind(i*180/NumPoints)*axmz))
      endfor;
    perim = op--reverse ip--cycle;
    fill perim withcolor Shade;
    draw perim
  endgroup
enddef;
    
beginfig(1); 
  f := ( 4, 1, 2 );
  numeric i, j, xsiz, ysiz, gstep, angl, vlen, fact, xa, xb, yb;
  color vdir, bpos, ap, bp, cp, dp;
  path base;
  xsiz = 2;
  ysiz = 3;
  gstep = 0.23;
  angl = 33;
  vlen = 0.45;
  fact = 1.3;
  vdir = vlen*( 0, cosd(angl), sind(angl) );
  xa = 0.5*gstep;
  xb = (xsiz+0.5)*gstep;
  yb = (ysiz+0.5)*gstep;
  ap = (-xa,-xa,0);
  bp = (-xa,yb,0);
  cp = (xb,yb,0);
  dp = (xb,-xa,0);
  for i=0 upto xsiz:
    for j=0 upto ysiz:
      bpos := gstep*(i,j,0);
      draw rp(bpos-fact*vdir)--rp(bpos);
      drawarrow rp(bpos-fact*vdir)--rp(bpos-vdir);
    endfor;
  endfor;
  base = rp(ap)--rp(bp)--rp(cp)--rp(dp)--cycle;
  fill base withcolor 0.7background;
  draw base;
  for i=0 upto xsiz:
    for j=0 upto ysiz:
      bpos := gstep*(i,j,0);
      draw rp(bpos+fact*vdir)--rp(bpos);
      drawarrow rp(bpos)--rp(bpos+vdir);
    endfor;
  endfor;
  draw rp(black)--rp( (0,0,vlen) ) dashed evenly;
  angline( fact*vdir, (0,0,vlen), black, 0.7*vlen, btex 57$^\circ$ etex, urt );
  label.top( btex $\vec{B}$ etex, rp(vdir) );
endfig;

beginfig(2); 						
  f := ( 4, 3, 2 );
  numeric i, numplates, xdelta, iangle, ang, oray, iray, numpois;
  numeric aang, imray, omray, mmray, axray, numouts;
  color shade, shadea, shadeb, p[];
  path fucir, outerarrow, theknob;
  numplates = 11;
  numouts = 2;
  xdelta = 0.1;
  iangle = 53;
  oray = 0.5;
  iray = 0.07;
  imray = 0.52;
  omray = 0.67;
  mmray = 0.64;
  axray = 0.04;
  numpois = 10;
  shadea = 0.8background;
  shadeb = 0.6background;
  fucir = rigorouscircle( (-(numplates-1)*xdelta,0,0), red, imray );
  outerarrow = subpath (1.7,2.3) of fucir;
  p1 = ( -(numplates-1)*xdelta, imray, 0 );
  p2 = ( -(numplates-1)*xdelta, omray, 0 );
  p3 = ( -(numplates-2)*xdelta, imray, 0 );
  p4 = ( -(numplates-2)*xdelta, omray, 0 );
  p5 = ( -(numplates-4)*xdelta, mmray, 0 );
  p6 = ( -(numplates-2)*xdelta, mmray, 0 );
  p7 = ( -(numplates+1)*xdelta, mmray, 0 );
  p8 = ( -(numplates-1)*xdelta, mmray, 0 );
  p9 = ( -xdelta, cosd( 0.5*iangle )*0.5*oray, sind( 0.5*iangle )*0.5*oray );
  p0 = ( 0, cosd( 0.5*iangle )*0.5*oray, -sind( 0.5*iangle )*0.5*oray );
  for i=1-numplates upto 0:
    if odd(i):
      shade := shadea;
      ang := iangle;
    else:
      shade := shadeb;
      ang := 0;
    fi;
    platecap( i*xdelta, shade, ang, oray, iray, numpois );
    background := shadea;
    rigorousdisc( 0, true, (i*xdelta,0,0), axray, xdelta*red );
    background := shadeb;
    rigorousdisc( 0, true, (i*xdelta,-oray,0), axray, xdelta*red );
    background := background;
  endfor;
  for i=1 upto numouts:
    background := shadea;
    rigorousdisc( 0, true, (i*xdelta,0,0), axray, xdelta*red );
%    background := shadeb;
%    rigorousdisc( 0, true, (i*xdelta,-oray,0), axray, xdelta*red );
  endfor;
  theknob = rigorousfearpath( ((numouts+1)*xdelta,0,0), 2*axray );
  fill theknob;
  draw theknob;
  background := background;
  draw rp(p1)--rp(p2);
  draw rp(p3)--rp(p4);
  drawarrow rp(p5)--rp(p6);
  drawarrow rp(p7)--rp(p8);
  drawarrow outerarrow;
  label.ulft( btex $d$ etex, rp(p7) );
  label( btex $A$ etex, rp(p9) );  
  label( btex $A$ etex, rp(p0) );  
endfig;

beginfig(3); 						
  f := ( 3, 2, 1 );
  numeric oray, iray;
  path opath, ipath;
  color lRpos, Rpos, Ppos, zpos, rpos, spos;
  oray = 0.5;
  iray = 0.2oray;
  opath = rigorouscircle( black, blue, oray );
  ipath = rigorouscircle( black, blue, 0.2oray );
  fill opath withcolor 0.7background;
  fill ipath withcolor background;
  draw opath;
  draw ipath;
  lRpos = ( oray, 0, 0 );
  Rpos = 0.5lRpos;
  Ppos = ( 0, 0, 2oray );
  zpos = 0.5Ppos;
  rpos = ( 0, iray, 0 );
  spos = -Rpos;
  drawarrow rp(black)--rp(lRpos);
  drawarrow rp(black)--rp(Ppos);
  drawarrow rp(black)--rp(rpos);
  dotlabel.lft( btex $P$ etex, rp(Ppos) );
  label.rt( btex $z$ etex, rp(zpos) );
  label.ulft( btex $R$ etex, rp(Rpos) );
  label.lrt( btex $r$ etex, rp(rpos) );
  label( btex $\sigma$ etex, rp(spos) );
endfig;

beginfig(4); 						
  f := ( 3, 4, 2 );
  numeric raya, rayb, rayc, lena, lenb;
  color baseac, baseav, basebc, basebv, rayabbeg, rayabend, lambdapos;
  raya = 0.017;
  rayb = 0.1;
  rayc = 0.2;
  lena = 1.4;
  lenb = 0.6;
  baseac = 0.5*lena*red;
  baseav = -lena*red;
  basebc = 0.5*lenb*red;
  basebv = -lenb*red;
  rayabbeg = ( 0.5*lenb, raya, 0 );
  rayabend = ( 0.5*lenb, rayc, 0 );
  lambdapos = ( 0, rayc, 0 );
  rigorousdisc( 0, true, basebc, rayc, basebv );
  rigorousdisc( 0, true, baseac, raya, baseav );
  rigorousdisc( 0, false, basebc, rayc, basebv );
  drawarrow rp( rayabbeg )--rp( rayabend );
  label.bot( btex $R$ etex, rp( 0.5[rayabbeg,rayabend] ) );
  label( btex $\lambda_t$ etex, rp ( lambdapos ) );
  label.lrt( btex $\lambda_f$ etex, rp( 0.3[baseac,basebc] ) );
endfig;

beginfig(5); 						
  f := ( 5, 5, 3 );
  numeric i, numpoints, ray, aux, blen, bang, fa, fb, fc, fd;
  path wirea, wireb, wirec, wired, wirepath;
  color wa[], wb[], btip, xpos;
  pen bigpen, smallpen;
  numpoints = 36;
  ray = 0.5;
  blen = 0.8;
  fa = -1.1;
  fb = -0.8;
  fc = -fb;
  fd = -fa;
  bang = -45;
  btip = blen*(0,cosd( bang ), sind( bang ) );
  for i=1 upto numpoints:
    aux := i*180/numpoints;
    wa[i] = ray*(red*cosd(aux)+green*sind(aux));
%    drawarrow rp( wa[i] )--rp( wa[i]+btip );
    wb[i] = ray*(-red*cosd(aux)-blue*sind(aux));
%    drawarrow rp( wb[i] )--rp( wb[i]+btip );
  endfor;
  wirepath = for i=1 upto numpoints:
    rp( wa[i] )-- endfor for i=1 upto numpoints:
    rp( wb[i] )-- endfor cycle;
  bigpen = pencircle scaled 3pt;
  smallpen = pencircle scaled 2.2pt;
  draw wirepath withpen bigpen;
  wirea = rp( wb[numpoints] ) for i=1 upto numpoints:
    --rp( wa[i] ) endfor;
  wireb = rp( wa[numpoints] ) for i=1 upto numpoints:
    --rp( wb[i] ) endfor;
  wirec = rp( wb[numpoints-1] )--wirea--rp( wb[1] );
%  wired = rp( wa[numpoints-1] )--wireb--rp( wa[1] );
  draw wireb withpen smallpen withcolor background;
%  cartaxes( 2ray, 2ray, 1.2ray );
  draw rp( 2*ray*red )--rp( -2*ray*red ) dashed evenly;
  draw rp( black )--rp( -ray*blue );
  draw rp( black )--rp( ray*green );
  for i=-1 upto 1:
    xpos := (i*ray,0,0);
    drawarrow rp( xpos+fa*btip )--rp( xpos+fb*btip );
    drawarrow rp( xpos+fb*btip )--rp( xpos+fc*btip );
    draw rp( xpos+fc*btip )--rp( xpos+fd*btip );
  endfor;
  squareangline( -ray*blue, ray*green, black, 0.1*ray);
  draw wirea withpen bigpen;
  draw wirec withpen smallpen withcolor background;
  label.urt( btex $\vec{B}$ etex, rp( btip ) );
endfig;

beginfig(6);
  path baset, sphere;
  color va, vb, vc, qpos;
  va = (0,0,0);
  vb = (0.5,0.5*sqrt(3),0);
  vc = (-0.5,0.5*sqrt(3),0);
  baset = rp( va )--rp( vb )--rp( vc )--cycle;
  draw baset;
  qpos = (0,0.333*sqrt(3),0.333*sqrt(1-0.333*sqrt(3)) );
  sphere = rigorousfearpath( qpos, 0.1 );
  fill sphere withcolor background;
  draw sphere;
  label( btex +$q$ etex, rp( qpos ) );
endfig;

beginfig(7);
  Spread := 50;
  f := (4,6,3.5);
  path fiodecima, fioabaixo, malhacub, rama, ramb, ramc;
  path malhatet, ramd, rame;
  pen grossa, nuclear, negat;
  color p[];
  numeric altura, margem, largneg, largpos, afast, lengax, lscal;
  pair topo, base, decima, abaixo, centerfont, centerneg, centerpos;
  pair a, b, c, d, tx, ty, tz;
  altura = 2cm;
  margem = 3mm;
  largneg = 9mm;
  largpos = 18mm;
  afast = 5cm;
  lengax = 0.7;
  lscal = 2;
  grossa = pencircle scaled 2.7pt;
  nuclear = pencircle scaled 2.1pt;
  negat = pencircle scaled 1.8pt;
  p1 = (1,1,-1);
  p2 = (-1,1,-1);
  p3 = (-1,-1,-1);
  p4 = (1,-1,-1);
  p5 = (1,1,1);
  p6 = (-1,1,1);
  p7 = (-1,-1,1);
  p8 = (1,-1,1);
  p0 = (0,0,1);
  p9 = (0,0,-1);
  z1 = rp(p1);
  z2 = rp(p2);
  z3 = rp(p3);
  z4 = rp(p4);
  z5 = rp(p5);
  z6 = rp(p6);
  z7 = rp(p7);
  z8 = rp(p8);
  z0 = rp(p0);
  z9 = rp(p9);
  tx = rp( p1 - (lengax,0,0) ); 
  ty = rp( p1 - (0,lengax,0) ); 
  tz = rp( p1 + (0,0,lengax) ); 
  topo = (xpart z0, altura + ypart z0);
  base = (xpart z9, -altura + ypart z9);
  decima = (afast + xpart topo, ypart topo);
  abaixo = (afast + xpart base, ypart base);
  centerfont = 0.5*(decima+abaixo);
  centerpos = (xpart centerfont, 0.5*margem + ypart centerfont);
  centerneg = (xpart centerfont, -0.5*margem + ypart centerfont);
  fiodecima = z0--topo--decima--centerpos;
  fioabaixo = z9--base--abaixo--centerneg;
  a = ( -0.5*largpos + xpart centerfont, ypart centerpos);
  b = ( 0.5*largpos + xpart centerfont, ypart centerpos);
  c = ( -0.5*largneg + xpart centerfont, ypart centerneg);
  d = ( 0.5*largneg + xpart centerfont, ypart centerneg);
  malhacub = z1--z2--z3--z4--z1--z5--z6--z7--z8--z5;
  malhatet = z2--z7--z4--z5--cycle;
  rama = z2--z6;
  ramb = z3--z7;
  ramc = z4--z8;
  ramd = z5--z7;
  rame = z2--z4;
  drawoptions( dashed evenly );
  draw malhacub;
  draw rama;
  draw ramb;
  draw ramc;
  drawoptions( withpen grossa );
  draw malhatet;
  draw ramd;
  draw rame;
%  drawoptions( withpen nuclear withcolor background );
%  draw malhatet;
%  draw ramd;
%  draw rame;
  drawoptions();
  draw fiodecima;
  draw fioabaixo;
  draw a--b;
  drawarrow z1--tx withpen negat;
  drawarrow z1--ty withpen negat;
  drawarrow z1--tz withpen negat;
  label.lrt(btex $x$ etex scaled lscal, tx);
  label.llft(btex $y$ etex scaled lscal, ty);
  label.lft(btex $z$ etex scaled lscal, tz);
  linecap := butt;
  draw c--d withpen negat;
endfig;

beginfig(8);
  Spread := 50;
  f := (6,6,3.5);
  numeric squareside, sray, cdiameter, cray, distmargin, linemargin, eloing;
  squareside = 1;
  cdiameter = squareside/0.88623;
  cray = cdiameter/2.0;
  sray = squareside/2.0;
  distmargin = 0.5;
  linemargin = 0.15;
  eloing = 0.4;
  pen pena, penb, penc;
  pena = pencircle scaled 1.3mm;
  penb = pencircle scaled 0.95mm;
  penc = pencircle scaled 0.4mm;
  color tracecol, outcol, incol;
  tracecol = 0.5white;
  outcol = black;
  incol = 0.876white;
  color va, vb, vc, vd, cen, bs, bc, fs, fc, ds, dc;
  va = (0,-sray,-sray);
  vb = (0,sray,-sray);
  vc = (0,sray,sray);
  vd = (0,-sray,sray);
  cen = (-distmargin,0,0);
  bs = (0,-sray-linemargin,0);
  bc = bs+cen;
  fs = (0,sray+eloing,0);
  ds = (0,sray+eloing-linemargin,0);
  fc = fs+cen;
  dc = ds+cen;
  path rigcirc, carrpath, square, circ, sline, cline, dline;
  square = rp(va)--rp(vb)--rp(vc)--rp(vd)--cycle;
  circ = goodcirclepath( cen, red, cray );
  rigcirc = rigorouscircle( cen, red, cray+linemargin );
  carrpath = subpath (1,2) of rigcirc;
  sline = rp(bs)--rp(fs);
  cline = rp(bc)--rp(fc);
  dline = rp(ds)--rp(dc);
  
  draw cline withpen penc withcolor tracecol;
  draw circ withpen pena withcolor outcol;
  draw circ withpen penb withcolor incol;
  draw sline withpen penc withcolor tracecol;
  draw square withpen pena withcolor outcol;
  draw square withpen penb withcolor incol;
  drawdblarrow dline withpen penc withcolor tracecol;
  
  pair begsarr, endsarr;
  begsarr = rp((0,-sray-linemargin,linemargin));
  endsarr = rp((0,-sray-linemargin,sray));
  drawarrow begsarr--endsarr withpen penc withcolor tracecol;
  drawarrow carrpath withpen penc withcolor tracecol;

  label.lft( btex $I$ etex, rp((0,-sray-linemargin,sray-linemargin)) );
  label.top( btex $+$ etex, point 1.5 of rigcirc );
  label.top( btex $d$ etex, rp(dc) );
endfig;

beginfig(9);
  Spread := 50;
  f := (6,6,3.5);
  numeric squareside, sray, cdiameter, cray, distmargin, linemargin, eloing;
  squareside = 1;
  cdiameter = squareside/0.88623;
  cray = cdiameter/2.0;
  sray = squareside/2.0;
  distmargin = 0.5;
  linemargin = 0.15;
  eloing = 0.4;
  pen pena, penb, penc;
  pena = pencircle scaled 1.3mm;
  penb = pencircle scaled 0.95mm;
  penc = pencircle scaled 0.4mm;
  color tracecol, outcol, incol;
  tracecol = 0.5white;
  outcol = black;
  incol = 0.876white;
  color va, vb, vc, vd, cen, bs, bc, fs, fc, ds, dc, os, oc, ncc;
  va = (0,-sray,-sray);
  vb = (0,sray,-sray);
  vc = (0,sray,sray);
  vd = (0,-sray,sray);
  cen = (-distmargin,0,0);
  bs = (0,-sray-linemargin,0);
  bc = bs+cen;
  fs = (0,sray+eloing,0);
  ds = (0,sray+eloing-linemargin,0);
  fc = fs+cen;
  dc = ds+cen;
  os = (linemargin,sray,0);
  oc = (-distmargin-cdiameter-linemargin,sray,0);
  ncc = (-distmargin-cray,sray,0);
  path rigcirc, carrpath, square, circ, sline, cline, dline, oline;
  square = rp(va)--rp(vb)--rp(vc)--rp(vd)--cycle;
  circ = goodcirclepath( ncc, green, cray );
  rigcirc = rigorouscircle( ncc, green, cray+linemargin );
  carrpath = subpath (3,4) of rigcirc;
  sline = rp(bs)--rp(fs);
  cline = rp(bc)--rp(fc);
  dline = rp(ds)--rp(dc);
  oline = rp(os)--rp(oc);
  
  draw oline withpen penc withcolor tracecol;
  draw cline withpen penc withcolor tracecol;
  draw circ withpen pena withcolor outcol;
  draw circ withpen penb withcolor incol;
  draw sline withpen penc withcolor tracecol;
  draw square withpen pena withcolor outcol;
  draw square withpen penb withcolor incol;
  drawdblarrow dline withpen penc withcolor tracecol;
  
  pair begsarr, endsarr;
  begsarr = rp((0,-sray-linemargin,linemargin));
  endsarr = rp((0,-sray-linemargin,sray));
  drawarrow begsarr--endsarr withpen penc withcolor tracecol;
  drawarrow carrpath withpen penc withcolor tracecol;

  label.lft( btex $I$ etex, rp((0,-sray-linemargin,sray-linemargin)) );
  label.top( btex $+$ etex, point 3.5 of rigcirc );
  label.bot( btex $d$ etex, rp(dc) );
endfig;

verbatimtex
\end{document}
etex

end;
