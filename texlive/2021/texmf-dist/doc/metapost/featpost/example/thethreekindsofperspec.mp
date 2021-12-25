% thethreekindsofperspec.mp
% L. Nobre G. 
% 2003

input featpost3Dplus2D;

% bash: export TEX=latex
% tcsh: setenv TEX latex

f := (4,2,2.5);
Spread := 45;

verbatimtex
\documentclass{article} 
\usepackage{beton,concmath,ccfonts} 
\begin{document} 
etex

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
  prll = rp(rll-pc);
  prlr = rp(rlr-pc);
  prur = rp(rur-pc);
  prul = rp(rul-pc);
  path projpath;
  projpath = prll--prlr--prur--prul--cycle;
  draw projpath withpen projpen;	     
  draw prll--rp(rll+(1-positfraction)*psdf) withpen lightpen;
  draw prlr--rp(rlr+(1-positfraction)*psdf) withpen lightpen;
  draw prur--rp(rur+(1-positfraction)*psdf) withpen lightpen;
  draw prul--rp(rul+(1-positfraction)*psdf) withpen lightpen;
  draw rp(black)--rp(psdf) withpen lightpen;
  draw rp(black) withpen dotpen;
  draw rp(psdf) withpen dotpen;
  label.rt( btex f etex, rp(psdf) );
%  label.bot( btex viewcentr etex, rp(black) );
endfig;

beginfig(2);
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
  label.rt( btex f etex, rp(psdf) );
  label.bot( btex viewcentr etex, rp(black) );
endfig;

beginfig(3);
  draw projecplanepath withpen realpen;
  draw realpath withpen realpen;
  draw rp(psdf) withpen lightpen;
  drawoptions( withpen lightpen );
  squareangline( psdf, psdv, black, rectsize );
  squareangline( psdf, psdu, black, rectsize );
  numeric dist;
  dist = conorm( psdf );
  color ml, mt, mr, mb;
  mb = 0.5[rll,rlr];
  ml = 0.5[rll,rul];
  mr = 0.5[rur,rlr];
  mt = 0.5[rur,rul];
  drawoptions( withpen projpen );
  angline( rll, mb, psdf, dist, "", top );
  angline( rll, ml, psdf, dist, "", top );
  angline( rur, mt, psdf, dist, "", top );
  angline( rur, mr, psdf, dist, "", top );  
  angline( mb, rlr, psdf, dist, "", top );
  angline( ml, rul, psdf, dist, "", top );
  angline( mr, rlr, psdf, dist, "", top );
  angline( mt, rul, psdf, dist, "", top );
  angline( black, mb, psdf, dist, "", top );
  angline( black, ml, psdf, dist, "", top );
  angline( black, mr, psdf, dist, "", top );
  angline( black, mt, psdf, dist, "", top );
  drawoptions();
  numeric marc, darc, dmarc, ddarc;
  marc = 0.01745*getangle(-psdf,mt-psdf)*dist;
  darc = 0.01745*getangle(-psdf,rll-psdf)*dist*0.707;
  color pro[], esf[], proc[], esfc[];
  pro1 = darc*(-psdv-psdu);
  pro2 = marc*(-psdu);
  pro3 = darc*(+psdv-psdu);
  pro4 = marc*(+psdv);
  pro5 = darc*(+psdv+psdu);
  pro6 = marc*(+psdu);
  pro7 = darc*(-psdv+psdu);
  pro8 = marc*(-psdv);
  esf1 = psdf+dist*N(rll-psdf);
  esf2 = psdf+dist*N(mb-psdf);
  esf3 = psdf+dist*N(rlr-psdf);
  esf4 = psdf+dist*N(mr-psdf);
  esf5 = psdf+dist*N(rur-psdf);
  esf6 = psdf+dist*N(mt-psdf);
  esf7 = psdf+dist*N(rul-psdf);
  esf8 = psdf+dist*N(ml-psdf);
  ddarc = conorm(pro1-esf1)/3;
  dmarc = conorm(pro2-esf2)/3;
  proc1 = pro1+ddarc*N(psdf);
  proc2 = pro2+dmarc*N(psdf);
  proc3 = pro3+ddarc*N(psdf);
  proc4 = pro4+dmarc*N(psdf);
  proc5 = pro5+ddarc*N(psdf);
  proc6 = pro6+dmarc*N(psdf);
  proc7 = pro7+ddarc*N(psdf);
  proc8 = pro8+dmarc*N(psdf);
  esfc1 = esf1+ddarc*N(rll-psdf);
  esfc2 = esf2+dmarc*N(mb-psdf);
  esfc3 = esf3+ddarc*N(rlr-psdf);
  esfc4 = esf4+dmarc*N(mr-psdf);
  esfc5 = esf5+ddarc*N(rur-psdf);
  esfc6 = esf6+dmarc*N(mt-psdf);
  esfc7 = esf7+ddarc*N(rul-psdf);
  esfc8 = esf8+dmarc*N(ml-psdf);
  numeric i;
  drawoptions( withpen lightpen );
  for i=1 upto 8:
    draw rp(esf[i])..controls rp(esfc[i]) and rp(proc[i])..rp(pro[i]);
    draw rp(psdf)--rp(esf[i]);
  endfor;
  draw rp(pro2)--rp(pro6);
  draw rp(pro4)--rp(pro8);
  draw rp(mb)--rp(mt);
  draw rp(ml)--rp(mr);
  drawoptions( withpen projpen );
  draw rp(pro1)...rp(pro2)..rp(pro3);
  draw rp(pro3)...rp(pro4)..rp(pro5);
  draw rp(pro5)...rp(pro6)..rp(pro7);
  draw rp(pro7)...rp(pro8)..rp(pro1);  
  draw rp(black)--rp(psdf) withpen lightpen;
  draw rp(black) withpen dotpen;
  draw rp(psdf) withpen dotpen;
  label.rt( btex f etex, rp(psdf) );
%  label.ulft( btex viewcentr etex, rp(black) );
endfig;

verbatimtex 
\end{document} 
etex

end;