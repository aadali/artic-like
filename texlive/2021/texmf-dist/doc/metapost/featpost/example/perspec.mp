% perspec.mp
% L. Nobre G. 
% 2004

input featpost3Dplus2D;

% bash: export TEX=latex
% tcsh: setenv TEX latex

f := (4,2,2.5);
Spread := 45;

verbatimtex \documentclass{article} \begin{document} etex

color psdf, psdu, psdv, rll, rlr, rur, rul, pll, plr, pur, pul;
color pc;
numeric projecplanesize, sizefraction, positfraction, aspectratio, rectsize;
pen realpen, lightpen, projpen, dotpen;
path projecplanepath, realpath;
psdf = 3.1*(0,1,0.37);
projecplanesize = 1.7;
sizefraction = 0.19;
positfraction = 0.65;
aspectratio = 1.0;
rectsize = 0.15;
realpen = pencircle scaled 2.6pt;
lightpen = pencircle scaled 0.55pt;
projpen = pencircle scaled 1.8pt;
dotpen = pencircle scaled 3.8pt;
psdv = N( (-Y(psdf), X(psdf), 0) );
psdu = ncrossprod( psdf, psdv );
pll = projecplanesize*(-aspectratio*psdv-psdu);
plr = projecplanesize*(+aspectratio*psdv-psdu);
pur = projecplanesize*(+aspectratio*psdv+psdu);
pul = projecplanesize*(-aspectratio*psdv+psdu);
projecplanepath = rp(pll)--rp(plr)--rp(pur)--rp(pul)--cycle;
pc = positfraction*psdf;
rll = pc+sizefraction*projecplanesize*(-psdv-psdu);
rlr = pc+sizefraction*projecplanesize*(+psdv-psdu);
rur = pc+sizefraction*projecplanesize*(+psdv+psdu);
rul = pc+sizefraction*projecplanesize*(-psdv+psdu);
realpath = rp(rll)--rp(rlr)--rp(rur)--rp(rul)--cycle;

beginfig(1);
  draw projecplanepath withpen realpen;
  draw realpath withpen realpen;
  draw rp(psdf) withpen lightpen;
  drawoptions( withpen lightpen );
  squareangline( psdf, psdv, black, rectsize );
  squareangline( psdf, psdu, black, rectsize );
  drawoptions();
  pair prll, prlr, prur, prul;
  LightSource:=psdf;
  prll = rp(projectpoint(black,rll));
  prlr = rp(projectpoint(black,rlr));
  prur = rp(projectpoint(black,rur));
  prul = rp(projectpoint(black,rul));
  path projpath;
  projpath = prll--prlr--prur--prul--cycle;
  draw projpath withpen projpen;	     
  draw prll--rp(psdf) withpen lightpen;
  draw prlr--rp(psdf) withpen lightpen;
  draw prur--rp(psdf) withpen lightpen;
  draw prul--rp(psdf) withpen lightpen;
  draw rp(black)--rp(psdf) withpen lightpen;
  draw rp(black) withpen dotpen;
  draw rp(psdf) withpen dotpen;
  label.rt( btex \texttt{f} etex, rp(psdf) );
  label.bot( btex \texttt{viewcentr} etex, rp(black) );
endfig;

verbatimtex 
\end{document} 
etex

end;
