% f1aula05.mp
% L. Nobre G.
% 2007

verbatimtex
\documentclass{beamer}
\usepackage{amsmath,amssymb}
\usepackage[squaren]{SIunits}
\usepackage[latin1]{inputenc}
\usepackage[portuges]{babel}
\newcommand{\veri}{\ensuremath{\hat{\text{i}}}}
\newcommand{\verj}{\ensuremath{\hat{\text{j}}}}
\newcommand{\verk}{\ensuremath{\hat{\text{k}}}}
\newcommand{\veru}{\ensuremath{\hat{\text{u}}}}
\newcommand{\emf}{\ensuremath{\mathcal{E}}}
\newcommand{\vecnorm}[1]{\left|\left|\vec{#1}\right|\right|}
\renewcommand{\arraystretch}{1.45}
\begin{document}
\frame{
etex

ahangle := 29;
ahlength := 4.7mm;
labeloffset := 3.5mm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% {} []
beginfig(1);
  pen basepen, vecpen;
  basepen = pencircle scaled 0.5mm;
  vecpen = pencircle scaled 1.15mm;
  numeric ray, lscale, ang, stepang, rang, actang;
  ray = 5cm;
  lscale = 2.5;
  rang = -9;
  stepang = 18;
  path thearc;
  thearc = halfcircle scaled 2ray rotated rang;
  for ang=0 step stepang until 180+rang:
    drawarrow origin--(ray*dir(ang)) withpen vecpen;
    actang := ang;
  endfor;
  draw thearc withpen basepen withcolor red;
  label.llft( btex $\vec{r}$ etex scaled lscale, 0.5*ray*dir(actang) );
  for ang=0 step stepang until 180+rang:
    drawarrow (ray,0)--(ray*dir(ang)) withpen vecpen withcolor blue;
  endfor;
  label.bot( btex $\vec{r}_0$ etex scaled lscale, (0.5ray,0) );
  picture colabel;
  colabel = thelabel.urt(
    btex $\Delta\vec{r}\rightarrow d\vec{r}\parallel\vec{v}\perp\vec{r}_0$ etex
         scaled lscale, (ray,0) );
  draw colabel withcolor blue;
endfig;

beginfig(2);
  pen basepen, vecpen;
  basepen = pencircle scaled 0.5mm;
  vecpen = pencircle scaled 1.15mm;
  numeric ray, lscale, ang, stepang, rang, actang, fac;
  ray = 5cm;
  lscale = 2.5;
  rang = -9;
  stepang = 18;
  fac = 0.7;     
  path thearc, actv;
  thearc = halfcircle scaled 2ray rotated rang;
  for ang=0 step stepang until 180+rang:
    drawarrow (ray,0)--((ray,0)+fac*ray*dir(ang+90)) withpen vecpen
    withcolor green;
    actv := ((ray*dir(ang))--(ray*(dir(ang)+fac*dir(ang+90))));
    drawarrow actv withpen basepen;
    actang := ang;
  endfor;
  draw thearc withpen basepen withcolor red;
  label.lrt( btex $\vec{v}$ etex scaled lscale, point 0.5 of actv );
  for ang=0 step stepang until 180+rang:
    drawarrow (ray,fac*ray)--((ray,0)+fac*ray*dir(ang+90)) withpen vecpen
                                                           withcolor blue;
  endfor;
  label.rt( btex $\vec{v}_0$ etex scaled lscale, (ray,0.5*fac*ray) );
  picture colabel, relabel;
  colabel = thelabel.rt(
    btex $\Delta\vec{v}\rightarrow d\vec{v}\parallel\vec{a}\perp\vec{v}_0$ etex
         scaled lscale, (ray,fac*ray) );
  draw colabel withcolor blue;
%   drawarrow ((0,fac*ray)--(fac*ray*dir(stepang+90))) shifted (ray,-fac*ray)
%                                         withpen vecpen withcolor blue+green;
%   relabel = thelabel.rt(
%     btex $\Delta\vec{v}\parallel\vec{a}$ etex scaled lscale, (ray,0) );
%   draw relabel withcolor blue+green;
endfig;


verbatimtex
}
\end{document}
etex

end.