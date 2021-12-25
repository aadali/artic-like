% aledlogo.mp
% L. Nobre G.
% 2007

input featpost3Dplus2D;

prologues := 1;

verbatimtex
\documentclass{article}
\usepackage{marvosym}
\begin{document}
etex

vardef calxtep( expr Beglim, Endlim, Instep ) =
  ((Endlim-Beglim)/floor((Endlim-Beglim)/Instep))
enddef;

def drawfila( expr H, Ray, Astart, Afinis, Astep, Sray ) =
  begingroup
    save a, c, p, n, i, s, l, r, u, d, inc, rad, ba, bb, vertex;
    numeric a, n, i, s, inc, sa, sb;
    color c, p, l, r, u, d, rad, vertex;
    path shape;
    boolean ba, bb;
    picture aux;
    inc = -5;
    sa = 2;
    sb = 0.27;
    vertex = (-0.7,0.9,0.15);
    n = floor((Afinis-Astart)/Astep);
    s = (Afinis-Astart)/n;
    for i=0 upto n:
      a := Astart+i*s;
      p := (Ray*cosd(a),Ray*sind(a),H);
      rad := (cosd(a),sind(a),0);
      u := Sray*(blue*cosd(inc)+rad*sind(inc));
      r := Sray*ncrossprod( p-f, u );
      l := -r;
      d := -u;
      ba := (ypart rp(p)) < (ypart rp(vertex)) 
                            + sa*((xpart rp(p))-(xpart rp(vertex)));
      bb := (ypart rp(p)) > (ypart rp(vertex))
			    + sb*((xpart rp(p))-(xpart rp(vertex)));
      shape := rp(p+r+d)--rp(p+r+u)--rp(p+l+u)--rp(p+l+d)--cycle;
      if ( uniformdeviate(1) < 0.5 ):                                        
%%    if ( uniformdeviate(1) < 0.7 ): %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    if ( ba and bb ): 
	c := black;
	unfill shape;
	draw shape withcolor c;
	aux := thelabel( btex \EyesDollar etex scaled (ps( p, 7 )), rp(p) );
	draw aux withcolor c;
%%%%	fill rp(p+r+d)--rp(p+r+u)--rp(p+l)--cycle withcolor c;
      else:
%%%%	c := 0.5white;
%%	c := 0.5red;
	c := black;
	unfill shape;
	draw shape withcolor c;
	aux := thelabel( btex \EURhv etex scaled (ps( p, 9 )), rp(p) );
	draw aux withcolor c;
%%%%	fill rp(p+l+d)--rp(p+l+u)--rp(p+r)--cycle withcolor c;
      fi;
%%    fill rp(p+r+u)--rp(p+l+u)--rp(p+d)--cycle withcolor c;
%%%%  fill rp(p+r)--rp(p+u)--rp(p+l)--rp(p+d)--cycle withcolor c;
%%%%  fill shape withcolor c;
%%%%  draw shape withcolor white withpen pensquare scaled 3pt;
%%%%  fill rigorousfearpath( p, Sray ) withcolor c;
    endfor
  endgroup
enddef;

%SphericalDistortion := true;

beginfig(1);
  numeric heig, bray, tray, vang, sray, ldis, fdis, cwid, cnum, ylen;
  heig = 0.35;
  bray = 0.95;
  tray = 2.25;
  vang = 95;
  sray = 0.029;
  ldis = 0.02;
  fdis = 0.075;
  cwid = 0.15;
  cnum = 2;
  ylen = -0.16;
  f := (0,ylen,heig);
  viewcentr := (0,tray,heig);
  color thered;
  thered = 0.66red;
  numeric degperrad, asector, s;
  degperrad = 180.0/3.14159;
  numeric a, h, astep, ray, raystep, hstep, abeg, aend, astart, afinis, awid;
  raystep = calxtep( bray, tray, fdis );
  hstep = raystep/(tray-bray);
  h = 1;
  abeg = 0.5*(180-vang);
  aend = 180-abeg;
  for ray=tray step -raystep until bray:
    awid := degperrad*cwid/ray;
    asector := (vang-cnum*awid)/(cnum+1);
    astep := degperrad*(ldis+2*sray)/ray;
    for s=0 upto cnum:
      astart := abeg+s*(asector+awid);
      afinis := astart+asector;
      drawfila( h, ray, astart, afinis, astep, sray );
    endfor;
    h := h - hstep;
  endfor;
  picture stor;
  stor = currentpicture;
  bboxmargin := 2cm;
  fill bbox stor withcolor thered;
  draw stor;  
endfig;

beginfig(2);
  numeric angout, anglow, rayout, distra, lifactor, lidist, madist;
  numeric maleng, numstp, madstp, malstp, mafrac, i, j, begy, endy, xcoord;
  color cenpos, vechor, vecver, yela, yelb, dblu, mgre;
  pair ulc, llc, urc;
  picture block;
  angout=30;
  anglow=-6;
  lifactor=2;
  lidist=3.5;
  mafrac=0.6;
  numstp=25;
  f:=(7,6,5);
  viewcentr:=(0,-1.4,0);
  Spread:=35;
  yela=(1,1,0.7);
  yelb=(1,1,0.45);
  background:=yelb;
  dblu=0.4blue;
  mgre=0.6green;
  rayout=sind(angout);
  distra=cosd(angout);
  llc = point 0 of VGAborder;
  ulc = point 3 of VGAborder;
  urc = point 2 of VGAborder;
  madist = abs( ulc-llc );
  maleng = abs( urc-ulc );
  malstp = maleng/numstp;
  madstp = 5*malstp/3;
  z2 = mafrac*(malstp,madstp);
  z1 = mafrac*(0,madstp);
  z3 = mafrac*(malstp,0);
  fill origin--z1--z2--z3--cycle shifted (0.5*(1-mafrac)*(malstp,madstp))
                                 withcolor mgre;
  block = currentpicture;
  currentpicture := nullpicture;
  fill VGAborder withcolor dblu;
  for i=0 upto numstp-1:
    begy:= (ypart ulc) + madstp*(1+normaldeviate);
    endy:= (ypart 0.5[llc,ulc]) + 4*madstp*normaldeviate;
    xcoord:= (xpart llc) + i*malstp;
    for j=begy step -madstp until endy:
      draw block shifted (xcoord,j);
    endfor;
  endfor;
  clip currentpicture to VGAborder;
  cenpos:=distra*(-cosd(30)*cosd(anglow),-sind(30)*cosd(anglow),sind(anglow));
  vechor:=rayout*(sind(30),-cosd(30),0);
  vecver:=rayout*(cosd(30)*sind(anglow),sind(30)*sind(anglow),cosd(anglow));
  drawoptions( withcolor yelb );
  whatisthis(cenpos,vechor,vecver,3lidist,1+3*(lifactor-1));
  fill rigorousfearpath( black, 1 ) withcolor yela;
  cenpos:=distra*(cosd(30)*cosd(anglow),-sind(30)*cosd(anglow),sind(anglow));
  vechor:=rayout*(sind(30),cosd(30),0);
  vecver:=rayout*(-cosd(30)*sind(anglow),sind(30)*sind(anglow),cosd(anglow));
  drawoptions( withcolor yelb );
  whatisthis(cenpos,vechor,vecver,lidist,lifactor);
  cenpos:=distra*(0,cosd(anglow),sind(anglow));
  vechor:=rayout*(-1,0,0);
  vecver:=rayout*(0,-sind(anglow),cosd(anglow));
  whatisthis(cenpos,vechor,vecver,lidist,lifactor);
  drawoptions();
  produce_auto_scale;
endfig;

end.