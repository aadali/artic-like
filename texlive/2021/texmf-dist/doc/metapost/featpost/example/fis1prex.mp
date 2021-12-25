% fis1prex.mp
% L. Nobre G. 
% 2003

input featpost3Dplus2D;

verbatimtex
\documentclass{article}
\usepackage{beton}
\usepackage{concmath}
\usepackage{ccfonts}
\begin{document}
etex

beginfig(1); 						
	numeric largtecto, molasfrac, altutecto, compcorda, espetecto;
	numeric largmolas, largnum, cordafrac, espesfrac, largmolfr;
	numeric u, alpha, i, numturns, alphafrac, alpharadi, rcf, rcd;
	pen mainp, auxip;
	path leftspr, rightsp, corda, boxtecto, traco, linhalpha, rcp;
	pair esqutecto, diretecto, esquemola, direimola, espessvec;
	largnum = 100;
	u = 1mm;
	molasfrac = 0.8;
	espesfrac = 0.1;
	cordafrac = 0.3;
	alphafrac = 0.3;
	largmolfr = 0.12;
	rcf = 0.1;
	numturns = 11;
	alpha = 120;
	mainp = pencircle scaled 0.4u;
	auxip = pencircle scaled 0.2u;
	largtecto = largnum*u;
	rcd = rcf*largtecto;
	largmolas = largmolfr*molasfrac*largtecto;
	compcorda = cordafrac*largtecto;
	espetecto = espesfrac*largtecto;
	altutecto = 0.5*largtecto*molasfrac*cosd(alpha/2)/sind(alpha/2);
	alpharadi = 2*alphafrac*altutecto;
	esqutecto = (-largtecto/2,altutecto);
	diretecto = (largtecto/2,altutecto);
	esquemola = (-largtecto*molasfrac/2,altutecto);
	direimola = (largtecto*molasfrac/2,altutecto);
	espessvec = (0,espetecto);
	boxtecto = esqutecto--
		   diretecto--
		  (diretecto+espessvec)--
		  (esqutecto+espessvec)--cycle;
	traco = (origin--(dir(45))) scaled 6in;
	for i=-largnum upto largnum:
		draw traco shifted (2*i*u,0) withpen auxip;
	endfor;
	clip currentpicture to boxtecto;
	pickup mainp;
	draw esqutecto--diretecto;
	leftspr = springpath(direimola,origin,numturns,largmolas,0.55);
	rightsp = springpath(origin,esquemola,numturns,largmolas,0.55);
	draw leftspr;
	draw rightsp;
	rcp = randomcirc( rcd, 0.1*rcd, 180 ) shifted (0,-compcorda);
	draw rcp;
	draw origin--(0,-compcorda) cutafter rcp;
	linhalpha = fullcircle scaled alpharadi;
	linhalpha := linhalpha cutafter rightsp;
	linhalpha := linhalpha cutbefore leftspr;
	draw linhalpha withpen auxip;
	label(btex m etex scaled 3, (0,-compcorda) );
	label.urt(btex $\alpha$ etex scaled 3, (0,alpharadi/2) );
endfig;

beginfig(3); 						
  numeric u, largtecto, altutecto, espetecto, largmolas, rcd, largnum;
  numeric spfraca, spfracb, angex, lenfactor, lenspringb, largfac;
  numeric lenfel, lenwei;
  pen mainp, auxip;
  path molapatha, molapathb, rcpa, rcpb, boxtecto, traco;
  pair esqutecto, diretecto, espessvec, fixpoi, springbextreme, fel, wei, res;
  color graycol;
  mainp = pencircle scaled 0.4mm;
  auxip = pencircle scaled 0.2mm;
  angex = 36;
  lenfactor = 1.35;
  largtecto = 100mm;
  espetecto = 10mm;
  altutecto = 50mm;
  largmolas = 7mm;
  rcd = 6mm;
  lenfel = 34.5mm;
  lenwei = 20mm;
  largnum = 100;
  spfraca = 0.65;
  graycol = 0.67white;
  largfac = 0.5*(1+lenfactor);
  fel = lenfel*dir(90+angex);
  wei = lenwei*down;
  res = fel+wei;
  lenspringb = altutecto*lenfactor;
  spfracb = 1-altutecto*(1-spfraca)/lenspringb;
  randomseed := 3976;
  esqutecto = (-largtecto/2,altutecto);
  diretecto = (largtecto/2,altutecto);
  espessvec = (0,espetecto);
  fixpoi = (0,altutecto);
  springbextreme = fixpoi+lenspringb*dir(angex-90);
  fel := springbextreme+fel;
  wei := springbextreme+wei;
  res := springbextreme+res;
  boxtecto = esqutecto--
    diretecto--
    (diretecto+espessvec)--
    (esqutecto+espessvec)--cycle;
  traco = (origin--(dir(45))) scaled 6in;
  for i=-largnum upto largnum:
    draw traco shifted (2*i*mm,0) withpen auxip;
  endfor;
  clip currentpicture to boxtecto;
  pickup mainp;
  draw esqutecto--diretecto;
  molapatha = springpath(fixpoi,origin,11,largmolas,spfraca);
  draw molapatha withcolor graycol;
  rcpa = randomcirc( rcd, 0.1*rcd, 24 );
  unfill rcpa;
  draw rcpa withcolor graycol;
  molapathb = springpath(fixpoi,springbextreme,11,largmolas/largfac,spfracb);
  draw molapathb;
  rcpb = rcpa rotated angex shifted springbextreme;
  unfill rcpb;
  draw rcpb;
  drawoptions( withcolor blue );
  drawarrow springbextreme--fel;
  drawarrow springbextreme--wei;
  drawarrow springbextreme--res;
  label.urt( btex $\vec{F}_{el}$ etex, fel );
  label.llft( btex $\vec{F}_{R}$ etex, res );
  label.ulft( btex $\vec{F}_{g}$ etex, wei );
endfig;

beginfig(2);
  numeric u, basesize, baseheight, alphang, thetang, sfearsize;
  numeric armsize, marginraise, totalraise, axeray, cardansize;
  color basepos, tdorigin, oneaxe, otheraxe, rotaxbascenter, axedir;
  color gluaxbascenter, rotrajcenter, cardanpos, refpoint, anotheraxe;
  color armdir, armbasonepos, armbastwopos, sfearpos;
  path glueaxebase, rotataxebase, rotrajector, axepath, sfearpath;
  numeric armray;
  u = 0.0125mm;
  basesize = 30u;
  baseheight = 3u;
  alphang = 20;
  thetang = -18;
  armsize = 25u;
  marginraise = 6u;
  axeray = 3u;
  cardansize = 3u;
  sfearsize = 4u;
  armray = 1.5u;
  totalraise = armsize*sind(alphang)+marginraise+cardansize*cosd(alphang);
  basepos = (0,0,-0.5baseheight);
  f := -2*(-2basesize,-basesize,-5baseheight);
  kindofcube(false, false, basepos, 0, 0, 0, 2basesize, 2basesize, baseheight );
  tdorigin = (0,0,0);
  oneaxe = (0,axeray/cosd(alphang),0);
  otheraxe = (1,0,0);
  anotheraxe = (0,cosd(alphang),-sind(alphang));
  gluaxbascenter = tdorigin - (0,totalraise*sind(alphang),0);
  glueaxebase = ellipticpath( gluaxbascenter, oneaxe, axeray*otheraxe );
  axedir = (0,sind(alphang),cosd(alphang));
  rotrajcenter = (0,0,totalraise);
  rotaxbascenter = rotrajcenter-axedir*cardansize;
  refpoint = rotrajcenter+axedir*cardansize;
  rotataxebase = rigorouscircle( rotaxbascenter, axedir, axeray );
  rotrajector = rigorouscircle( rotrajcenter, axedir, armsize );
  axepath = twocyclestogether( glueaxebase, rotataxebase );
  fill axepath withcolor background;
  draw axepath;
  draw rotataxebase;
  drawarrow subpath (3,10.5) of rotrajector dashed evenly;
%  dotlabel.top( btex 1 etex, point 1 of rotrajector );
%  dotlabel.top( btex 5 etex, point 5 of rotrajector );
%  dotlabel.top( btex 3 etex, point 3 of rotrajector );
  kindofcube(false,false,rotrajcenter,90,90-alphang,thetang,2cardansize,2cardansize,2cardansize);
  draw rp(refpoint);
  armdir := otheraxe*cosd(thetang)+anotheraxe*sind(thetang);
  armbasonepos = rotrajcenter + cardansize*armdir;
  armbastwopos = rotrajcenter + (armsize-sfearsize)*armdir; % NOT RIGOROUS
  sfearpos = rotrajcenter + armsize*armdir;
  rigorousdisc( 0, true, armbasonepos, armray, armbastwopos-armbasonepos );
  sfearpath = rigorousfearpath( sfearpos, sfearsize );
  fill sfearpath withcolor background;
  draw sfearpath;
endfig;

verbatimtex
\end{document}
etex

end;
