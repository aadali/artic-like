% acoplanv.mp
% L. Nobre G.
% 2003

verbatimtex 
\documentclass{article} 
\usepackage{mathpazo} 
\begin{document} 
etex

beginfig(1);
  numeric overallang, torsionang, u, curvfact;
  u = 1mm;
  overallang = -22;
  torsionang = 40;
  curvfact = 0.6;
  numeric armsize, roundray, direcstep;
  armsize = 25u;
  roundray = 10u;
  direcstep = 2.8u;
  pair restaxedir, restdir, torsaxedir, torsdir;
  restaxedir = dir( 90+overallang );
  restdir = dir( overallang );
  torsaxedir = dir( 90+overallang+torsionang );
  torsdir = dir( overallang+torsionang );
  numeric leng;  
  pair midpos;
  pickup pencircle scaled 0.5u;
  
  drawoptions( withcolor 0.52white );
  for i=1.5*direcstep step direcstep until roundray: 
    leng := roundray +-+ i;
    midpos := i*restaxedir;
    draw (midpos-leng*restdir)--(midpos+leng*restdir);
    midpos := -i*restaxedir;
    draw (midpos-leng*restdir)--(midpos+leng*restdir);
  endfor;
  path minuslin, plusline;
  plusline = (0.5*direcstep*restaxedir-armsize*restdir)--
             (0.5*direcstep*restaxedir+armsize*restdir);
  minuslin = (-0.5*direcstep*restaxedir-armsize*restdir)--
             (-0.5*direcstep*restaxedir+armsize*restdir);
  draw plusline;
  draw minuslin;
  
  drawoptions( withcolor black );
  for i=1.5*direcstep step direcstep until roundray: 
    leng := roundray +-+ i;
    midpos := i*torsaxedir;
    draw (midpos-leng*torsdir)--(midpos+leng*torsdir);
    midpos := -i*torsaxedir;
    draw (midpos-leng*torsdir)--(midpos+leng*torsdir);
  endfor;
  path minuscur, pluscurv;
  numeric curvheig;
  curvheig = curvfact*armsize*sind( torsionang )/cosd( torsionang );
  leng := roundray +-+ (0.5*direcstep);
  minuscur = ((-curvheig-0.5*direcstep)*restaxedir-armsize*restdir){restdir}..
             (-0.5*direcstep*torsaxedir-leng*torsdir)--
             (-0.5*direcstep*torsaxedir+leng*torsdir)..
             {restdir}((curvheig-0.5*direcstep)*restaxedir+armsize*restdir);
  pluscurv = ((-curvheig+0.5*direcstep)*restaxedir-armsize*restdir){restdir}..
             (0.5*direcstep*torsaxedir-leng*torsdir)--
             (0.5*direcstep*torsaxedir+leng*torsdir)..
             {restdir}((curvheig+0.5*direcstep)*restaxedir+armsize*restdir);
  draw minuscur;
  draw pluscurv;
  
  drawoptions( withcolor red );
  path velbase;
  for i=roundray+direcstep step direcstep until armsize:
    velbase := (i*restdir)--(i*restdir+curvheig*restaxedir);
    drawarrow velbase cutbefore plusline cutafter minuscur;
  endfor;
endfig;

verbatimtex 
\end{document} 
etex

bye;
