% ledlogo.mp
% L. Nobre G.
% 2007

beginfig(1);
  numeric u, penangle, peneccen, penminor, letterheig, letterwid, letterang;
  u=2mm; 
  penangle=112;
  peneccen=2.1;
  penminor=2u;
  letterheig=15u;
  letterwid=12u;
  letterang=-25;
  numeric hdist, vdist, letterslantang, margindist, humanang, humanfactor;
  hdist=7u;
  vdist=7u;
  letterslantang=20;
  margindist=4.7u;
  humanang=7;
  humanfactor=0.2;
  pen paintpen;
  paintpen = pencircle scaled penminor xscaled peneccen rotated penangle;
  pickup paintpen;
  pair letterup, letterig, dirabove, dirbelow;
  letterup=dir(90+letterang-letterslantang);
  letterig=dir(letterang);
  dirabove=dir(90+letterang-letterslantang+humanang);
  dirbelow=dir(90+letterang-letterslantang-humanang);
  z0=origin;
  z1=z0-letterheig*letterup;
  z2=z1+letterwid*letterig;
  z3=z0+hdist*letterig+0.5*letterheig*letterup;
  z4=z3-letterheig*letterup;
  z5=z4+letterwid*letterig;
  z6=z3+hdist*letterig+0.5*letterheig*letterup;
  z7=z6-letterheig*letterup;
  z8=z6-(letterheig-margindist)*letterup;
  z9=z3+letterwid*letterig;
  z10=z9+margindist*letterig;
  z11=z0-(letterheig-humanfactor*margindist)*letterup;
  z12=z1+humanfactor*margindist*letterig;
  z13=z3-(letterheig-humanfactor*margindist)*letterup;
  z14=z4+humanfactor*margindist*letterig;
  z15=z7-2*humanfactor*hdist*letterig;
  drawoptions( withcolor 0.876*red+0.7*green );
%   draw z0--z1--z2;
%   draw z9--z3--z4--z5;
%   draw z8---z6{letterig}..{unitvector(z2-z5)}z10..z7;
  draw z0{-dirabove}..z11...z12..{letterig}z2;
  draw z9--z3{-dirabove}..z13...z14..{letterig}z5;
  draw z8{dirbelow}..z6{letterig}..{unitvector(z2-z5)}z10..z15;
endfig;

verbatimtex
\documentclass{article}
\usepackage[T1]{fontenc}
\renewcommand*\rmdefault{cmdh}
\begin{document}
etex

beginfig(2);
  picture led;
  led = thelabel( btex led etex, origin );
  pair ledshift;
  ledshift = (0.3pt,-4.4pt);
  numeric dc, db, da;
  dc = 24pt;
  db = 6.9pt;
  da = 4.98pt;
  fill fullcircle scaled dc withcolor red;
%  fill fullcircle scaled db withcolor (red+green);
  fill fullcircle scaled db withcolor 0.5green;
  draw led shifted ledshift withcolor (red+green);
endfig;

prologues := 1;

beginfig(3);
  numeric u, altur, diam, espessum, margem, distdisp, separ, frac;
  color vermel, amarel;
  vermel = 0.75red;
  amarel = 0.876*red+0.7*green;
  u = 5cm;
  diam = 0.23u;
  margem = 0.85diam;
  altur = u-margem;
  frac = 0.7;
  numeric auxdist, lencirc, orientang, limitshift, i;
  limitshift = 1.6u;
  auxdist = (1u-2margem)*0.5;
  distdisp = auxdist*sqrt(2) - diam;
  3*separ+2*espessum = distdisp;
  separ = frac*espessum;
  lencirc = length fullcircle;
  path uma, umb;
  uma = (0.5separ,-0.5diam)--(0.5separ+espessum,-0.5diam)
        --(0.5separ+espessum,0.5diam)--(0.5separ,0.5diam)--cycle;
  umb = uma rotated 180;
  fill uma withcolor amarel;
  fill umb withcolor amarel;
  picture parduns;
  parduns = currentpicture;
  currentpicture := nullpicture;
  path botcirc, caminhescudo, circulbase;
  circulbase = fullcircle scaled diam;
  botcirc = subpath (0.5*lencirc,lencirc) of fullcircle scaled u;
  caminhescudo = botcirc--(0.5u,altur)--(-0.5u,altur)--cycle;
  fill caminhescudo withcolor vermel;
  z1 = ( auxdist, altur-margem );
  z2 = ( -auxdist, altur-margem );
  z3 = ( auxdist, 0 );
  z4 = ( -auxdist, 0 );
  z5 = ( 0, auxdist );
  fill circulbase shifted z1 withcolor amarel;
  fill circulbase shifted z2 withcolor amarel;
  fill circulbase shifted z3 withcolor amarel;
  fill circulbase shifted z4 withcolor amarel;
  fill circulbase shifted z5 withcolor amarel;
  z6 = 0.5[z5,z1];
  z7 = 0.5[z5,z2];
  z8 = 0.5[z5,z3];
  z9 = 0.5[z5,z4];
  draw parduns rotated 45 shifted z6;
  draw parduns rotated -45 shifted z7;
  draw parduns rotated -45 shifted z8;
  draw parduns rotated 45 shifted z9;
  numeric actshift;
  for orientang=45 step 90 until 360:
    actshift := auxdist*sqrt(2) + 0.5diam + separ;
    forever:
      if uniformdeviate(1)<0.5:
	actshift := actshift + espessum + 0.5separ;
	draw parduns shifted (actshift,0) rotated orientang shifted z5;
        exitif actshift>limitshift;
	actshift := actshift + espessum + 0.5separ;
      else:
	actshift := actshift + 0.5diam;
        fill circulbase shifted (actshift,0) rotated orientang shifted z5
	                                                 withcolor amarel;
        exitif actshift>limitshift;
        actshift := actshift + 0.5diam;
      fi;
      actshift := actshift + separ;
    endfor;
  endfor;  
endfig;

beginfig(4);
  numeric u, hthick, vthick, lspace, espace, dspace, h, lbase, dbase;
  numeric dinray, halfh, diffe;
  u = 0.5mm;
  hthick = 20u;
  vthick = 10u;
  lspace = 15u;
  espace = 8u;
  dspace = 8u;
  h = 3*vthick+2*espace;
  lbase = 38u;
  dbase = 16u;
  halfh = h/2;
  dinray = halfh-vthick;
  diffe = 2*vthick+espace-halfh;
  path pathl, pathe, pathd, pathu;
  z0 = (0,0);
  z1 = (lbase,0);
  z2 = (lbase,vthick);
  z3 = (hthick,vthick);
  z4 = (hthick,h);
  z5 = (0,h);
  pathl = z0--z1--z2--z3--z4--z5--cycle;
  z6 = (lbase+lspace,0);
  z7 = z6 + (h,0);
  z8 = z7 + (0,vthick);
  z9 = z6 + (vthick,vthick);
  z10= z9 + (0,espace);
  z11= z10+ (vthick,0);
  z12= z11+ (-vthick,vthick);
  z13= z6 + (-diffe,2*vthick+espace);
  z14=z13+(diffe,-diffe);
  pathe = z6--z7--z8--z9--z10--z11--z12--z13--z14--cycle;
  z15= z7 + (dspace,0);
  z16= z15+ (halfh,halfh);
  z17= z15+ (0,h);
  z18= z17- (halfh,0);
  z19= z18+ (halfh,-halfh);
  z20= z15+ (0,h-vthick);
  z21= z20+ (dinray,-dinray);
  z22= z15+ (0,vthick); 
  pathd = z15..z16..z17---z18--z19--z20..z21..z22--cycle;
  z24= z6 + (-halfh,h);
  z25= z24+ (h,0);
  z26= z25+ (vthick,-vthick);
  z27= z26- (h,0);
  pathu = z24--z25--z26--z27--cycle;
  fill pathl withcolor 0.78red;
  fill pathe;
  fill pathd;
  fill pathu;
endfig;

end.
  
